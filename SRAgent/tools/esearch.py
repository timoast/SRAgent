# import
## batteries
import os
import sys
import time
from datetime import datetime, timedelta
from typing import Annotated, List, Dict, Optional
from urllib.error import HTTPError
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_core.runnables.config import RunnableConfig
## package
from SRAgent.tools.utils import set_entrez_access
from SRAgent.db.connect import db_connect
from SRAgent.db.get import db_get_entrez_ids
from SRAgent.organisms import OrganismEnum

# functions

def to_sci_name(organism: str) -> str:
    """
    Convert organism name to scientific name using OrganismEnum.
    """
    organism_str = organism.replace(" ", "_").upper()
    
    try:
        enum_name = OrganismEnum[organism_str].value
        return f'"{enum_name}"'
    except KeyError:
        raise ValueError(f"Organism '{organism}' not found in OrganismEnum.")

## default time span limits
#MAX_DATE = (datetime.now()).strftime('%Y/%m/%d')
#MIN_DATE = (datetime.now() - timedelta(days=7 * 365)).strftime('%Y/%m/%d')

@tool 
def esearch_scrna(
    query_terms: Annotated[List[str], "Entrez query terms"]=["10X Genomics", "single cell RNA sequencing", "single cell RNA-seq"],
    database: Annotated[str, "Database name ('sra' or 'gds')"]="sra",
    organisms: Annotated[List[str], "List of organisms to search."]=["human", "mouse"],
    max_ids: Annotated[Optional[int], "Maximum number of IDs to return."]=10,
    config: RunnableConfig=None,
    )-> Annotated[List[str], "Entrez IDs of database records"]:
    """
    Find scRNA-seq datasets in the SRA or GEO databases.
    """
    esearch_query = ""    
   
    # check if query terms are provided
    query_terms = " OR ".join([f'"{x}"' for x in query_terms])
    esearch_query += f"({query_terms})"
    
    # add date range, if provided in the config
    min_date = config.get("configurable", {}).get("min_date")
    max_date = config.get("configurable", {}).get("max_date")
    if min_date and max_date:
        date_range = f"{min_date}:{max_date}[PDAT]"
        esearch_query += f" AND ({date_range})"
    
    # add organism
    ## override organisms, if provided in the config
    if config.get("configurable", {}).get("organisms"):
        organisms = config["configurable"]["organisms"]
    ## create organism query
    organisms = " OR ".join([f"{to_sci_name(x)}[Organism]" for x in organisms])
    if organisms:
        esearch_query += f" AND ({organisms})"

    # other filters
    esearch_query += ' AND "transcriptomic single cell"[Source]'
    esearch_query += ' AND "public"[Access]'
    esearch_query += ' AND "has data"[Properties]'
    esearch_query += ' AND "library layout paired"[Filter]'
    esearch_query += ' AND "platform illumina"[Filter]'
    esearch_query += ' AND "sra bioproject"[Filter]'
    ## library prep methods to exclude
    esearch_query += ' NOT ("Smart-seq" OR "Smart-seq2" OR "Smart-seq3" OR "MARS-seq")'

    # override max_ids parameter with config
    if config.get("configurable", {}).get("max_datasets"):
        max_ids = config["configurable"]["max_datasets"]

    # use database to filter existing?
    filter_existing = config.get("configurable", {}).get("use_database", False)

    # return entrez IDs 
    return esearch_batch(esearch_query, database, max_ids=max_ids, filter_existing=filter_existing)

def esearch_batch(
    esearch_query: str, 
    database: str, 
    max_ids: Optional[int]=None,
    verbose: bool=False, 
    filter_existing: bool=False,
    max_retries: int=3, 
    base_delay: float=3.0
    ) -> List[str]:
    # get existing Entrez IDs
    existing_ids = set()
    if filter_existing:
        with db_connect() as conn:
            existing_ids = {str(x) for x in db_get_entrez_ids(conn=conn, database=database)}
    
    # search for novel Entrez IDs
    ids = []
    retstart = 0
    retmax = min(10000, max_ids) if max_ids else 10000  # NCBI limit is 10000 per request
    total_retmax = max_ids if max_ids else float('inf')
    
    while True:
        for attempt in range(max_retries):
            set_entrez_access()
            try:
                search_handle = Entrez.esearch(
                    db=database, 
                    term=esearch_query, 
                    retstart=retstart,
                    retmax=retmax,
                    sort='pub+date'
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                # Add new IDs that aren't in existing_ids
                new_ids = [x for x in search_results['IdList'] if x not in existing_ids]
                ids.extend(list(set(new_ids)))
                
                # Update retstart for next batch
                retstart += retmax
                
                # Sleep to respect NCBI rate limits
                time.sleep(0.34)
                break
                
            except HTTPError as e:
                if e.code == 429 and attempt < max_retries - 1:
                    wait_time = base_delay * 2 ** attempt
                    if verbose:
                        print(f"Got HTTP 429; retrying in {wait_time} s...", file=sys.stderr)
                    time.sleep(wait_time)
                else:
                    print(f"Error searching {database} with query: {esearch_query}: {str(e)}")
                    return list(set(ids))
            except Exception as e:
               print(f"Error searching {database} with query: {esearch_query}: {str(e)}")
               return list(set(ids))
        else:
            break
            
        # Check if we've reached our target number of IDs
        if len(ids) >= total_retmax:
            break
            
        # Check if we've retrieved all available results
        if retstart >= int(search_results['Count']):
            break
            
    # Remove duplicates and limit to max_ids if specified
    ids = list(set(ids))
    if max_ids:
        ids = ids[:max_ids]
    return ids

@tool 
def esearch(
    esearch_query: Annotated[str, "Entrez query string."],
    database: Annotated[str, "Database name (e.g., sra, gds, or pubmed)"],
    config: RunnableConfig,
    )-> Annotated[str, "Entrez IDs of database records"]:
    """
    Run an Entrez search query and return the Entrez IDs of the results.
    Example query for single cell RNA-seq:
        `("single cell"[Title] OR "single-cell"[Title] OR "scRNA-seq"[Title])`
    Example query for an ENA accession number (database = sra):
        `ERX13336121`
    Example query for a GEO accession number (database = gds):
        `GSE51372`
    """
    # debug model
    max_records = 3 if os.getenv("DYNACONF", "").lower() == "test" else None

    # check input
    if esearch_query == "":
        return "Please provide a valid query."
    for x in ["SRR", "ERR", "GSE", "GSM", "GDS", "ERX", "DRR", "PRJ", "SAM", "SRP", "SRX"]:
        if esearch_query == x:
            return f"Invalid query: {esearch_query}"
    
    # query
    records = []
    retstart = 0
    retmax = 50
    max_retries = 3
    base_delay = 3.0
    search_results = None
    while True:
        for attempt in range(max_retries):
            set_entrez_access()
            try:
                search_handle = Entrez.esearch(db=database, term=esearch_query, retstart=retstart, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                for k in ["RetMax","RetStart"]:
                    if k in search_results: del search_results[k]
                records.append(str(search_results))
                retstart += retmax
                time.sleep(0.5)
                break
            except HTTPError as e:
                if e.code == 429 and attempt < max_retries - 1:
                    time.sleep(base_delay * 2 ** attempt)
                else:
                    msg = f"HTTP Error searching {database} with query {esearch_query}: {e}"
                    print(msg, file=sys.stderr)
                    return msg
            except Exception as e:
                msg = f"Error searching {database} with query {esearch_query}: {str(e)}"
                print(msg, file=sys.stderr)
                return msg
        else:
            break
        if max_records and len(records) >= max_records:
            break
        if search_results is None or retstart >= int(search_results['Count']):
            break
        
    # return records
    if len(records) == 0:
        return f"No records found for query: {esearch_query}"
    if max_records:
        records = records[:max_records] 
    return ", ".join([str(x) for x in records])


if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv(override=True)
    Entrez.email = os.getenv("EMAIL")

    # query for scRNA-seq 
    config = {"configurable": {
        "organisms": ["human", "mouse", "rat", "dog"],
        "min_date": "2021/01/01",
        "max_date": "2025/12/31",
    }}
    query = '("single cell RNA sequencing" OR "single cell RNA-seq")'
    #query = '("bulk RNA sequencing")'
    input = {"esearch_query" : query, "database" : "sra"}
    #input = {"esearch_query" : query, "database" : "gds"}
    #input = {"organisms" : ["Homo sapien", "Mus musculus"]} 
    #input = {"organisms" : ["yeast"], "max_ids" : 10000}
    print(esearch_scrna.invoke(input, config=config))

    # esearch accession
    #input = {"esearch_query" : "GSE51372", "database" : "sra"}
    #input = {"esearch_query" : "GSE121737", "database" : "gds"}
    #input = {"esearch_query" : "35447314", "database" : "sra"}
    #input = {"esearch_query" : "35447314", "database" : "gds"}
    #print(esearch.invoke(input))


