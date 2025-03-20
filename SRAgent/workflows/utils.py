import os
import asyncio
import xml.etree.ElementTree as ET
import aiohttp

async def entrez_id_to_srx_async(ids, email, api_key=None, max_concurrent=5):
    """
    Asynchronously convert SRA internal IDs or Entrez IDs to SRX/ERX accessions
    
    Parameters:
    -----------
    ids : str or list
        Single ID or list of IDs to convert
    email : str
        Email address (required by NCBI for Entrez queries)
    api_key : str, optional
        NCBI API key for higher request limits
    max_concurrent : int, optional
        Maximum number of concurrent requests to NCBI
        
    Returns:
    --------
    dict
        Dictionary mapping input IDs to their corresponding SRX/ERX accessions
        If no accession is found, the value will be None
    """
    # Ensure ids is a list
    if isinstance(ids, str):
        ids = [ids]
    
    # Base URLs for Entrez API
    elink_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    efetch_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    esummary_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    
    # Common parameters for all requests
    base_params = {
        'tool': 'python_async_script',
        'email': email
    }
    
    if api_key:
        base_params['api_key'] = api_key
    
    # Create a semaphore to limit concurrent connections
    semaphore = asyncio.Semaphore(max_concurrent)
    
    async def fetch_url(session, url, params):
        """Helper function to fetch and parse from a URL with parameters"""
        async with semaphore:
            merged_params = {**base_params, **params}
            async with session.get(url, params=merged_params) as response:
                if response.status != 200:
                    text = await response.text()
                    raise Exception(f"HTTP error {response.status}: {text}")
                return await response.text()
    
    async def direct_sra_fetch(session, sra_id):
        """
        Get SRX/ERX accession directly from SRA using efetch
        This approach treats the ID as an SRA internal ID
        """
        params = {
            'db': 'sra',
            'id': sra_id,
            'retmode': 'xml'
        }
        
        try:
            xml_text = await fetch_url(session, efetch_base_url, params)
            root = ET.fromstring(xml_text)
            
            accessions = []
            # Look for experiment accessions in various places
            for exp in root.findall(".//EXPERIMENT"):
                accession = exp.get("accession")
                if accession and (accession.startswith("SRX") or accession.startswith("ERX")):
                    accessions.append(accession)
            
            # Look for experiment references
            for exp_ref in root.findall(".//EXPERIMENT_REF"):
                accession = exp_ref.get("accession")
                if accession and (accession.startswith("SRX") or accession.startswith("ERX")):
                    accessions.append(accession)
            
            return accessions
        except Exception as e:
            print(f"Error in direct SRA fetch for ID {sra_id}: {str(e)}")
            return []
    
    async def try_sra_summary(session, sra_id):
        """
        Try using esummary to get the SRX/ERX accession
        Sometimes works when efetch doesn't
        """
        params = {
            'db': 'sra',
            'id': sra_id,
            'retmode': 'xml'
        }
        
        try:
            xml_text = await fetch_url(session, esummary_base_url, params)
            root = ET.fromstring(xml_text)
            
            accessions = []
            
            # Process DocSum elements which contain experiment info
            for doc_sum in root.findall(".//DocSum"):
                # Look for ExpAcc or Experiment field
                for item in doc_sum.findall("./Item"):
                    if item.get("Name") == "ExpAcc" and item.text:
                        acc = item.text
                        if acc.startswith("SRX") or acc.startswith("ERX"):
                            accessions.append(acc)
                    
                    # Sometimes the experiment info is in ExpXml
                    elif item.get("Name") == "ExpXml" and item.text:
                        # Try to extract SRX/ERX from the XML text
                        import re
                        exp_matches = re.findall(r'accession="(SRX\d+|ERX\d+)"', item.text)
                        accessions.extend(exp_matches)
            
            return accessions
        except Exception as e:
            print(f"Error in SRA summary for ID {sra_id}: {str(e)}")
            return []
    
    async def get_sra_links(session, entrez_id, from_db="nucleotide"):
        """Get SRA links for a given Entrez ID (original method)"""
        params = {
            'dbfrom': from_db,
            'db': 'sra',
            'id': entrez_id
        }
        
        xml_text = await fetch_url(session, elink_base_url, params)
        root = ET.fromstring(xml_text)
        
        sra_links = []
        for linkset in root.findall(".//LinkSetDb"):
            if linkset.find("DbTo") is not None and linkset.find("DbTo").text == "sra":
                for link in linkset.findall("Link/Id"):
                    sra_links.append(link.text)
        
        return sra_links
    
    async def process_id(session, input_id):
        """Process a single ID to get SRX/ERX accessions"""
        try:
            # First, try direct SRA fetch (treating it as an SRA internal ID)
            accessions = await direct_sra_fetch(session, input_id)
            
            if accessions:
                return input_id, accessions
            
            # If that fails, try esummary
            accessions = await try_sra_summary(session, input_id)
            
            if accessions:
                return input_id, accessions
            
            # If direct methods fail, try the original approach (treating as Entrez ID)
            # Try different source databases
            databases = ["sra", "nucleotide", "gene", "bioproject", "biosample"]
            
            for db in databases:
                sra_links = await get_sra_links(session, input_id, db)
                
                if sra_links:
                    # For each SRA link, try to get SRX/ERX accessions
                    all_accessions = []
                    for sra_id in sra_links:
                        direct_results = await direct_sra_fetch(session, sra_id)
                        if direct_results:
                            all_accessions.extend(direct_results)
                    
                    if all_accessions:
                        return input_id, all_accessions
            
            # If all methods fail, return None
            return input_id, None
            
        except Exception as e:
            print(f"Error processing ID {input_id}: {str(e)}")
            return input_id, None
    
    async with aiohttp.ClientSession() as session:
        # Create tasks for all IDs
        tasks = [process_id(session, input_id) for input_id in ids]
        results = await asyncio.gather(*tasks)
        
        # Convert results to dictionary
        return dict(results)

# Example usage
async def main():
    from dotenv import load_dotenv
    load_dotenv(override=True)
    
    email = os.getenv("EMAIL")
    ids_to_convert = ["30604662", "21712488", "28707082"]
    
    result = await entrez_id_to_srx_async(ids_to_convert, email)
    
    for input_id, accessions in result.items():
        if accessions:
            print(f"ID {input_id} corresponds to: {', '.join(accessions)}")
        else:
            print(f"No SRX/ERX accession found for ID {input_id}")

if __name__ == "__main__":
    asyncio.run(main())