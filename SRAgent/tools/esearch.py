# import
## batteries
import os
import sys
import time
from datetime import datetime, timedelta
from typing import Annotated, List, Dict, Optional
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
## package
from SRAgent.tools.utils import set_entrez_access
from SRAgent.db.connect import db_connect
from SRAgent.db.get import db_get_unprocessed_records

# functions
ORGANISMS = {
    'human': 'Homo sapiens',
    'mouse': 'Mus musculus',
    'rat': 'Rattus norvegicus',
    'monkey': 'Simiiformes', 
    'macaque': 'Macaca mulatta',
    'marmoset': 'Callithrix jacchus',
    'horse': 'Equus caballus',
    'dog': 'Canis lupus',
    'bovine': 'Bos taurus',
    'chicken': 'Gallus gallus',
    'sheep': 'Ovis aries',
    'pig': 'Sus scrofa',
    'fruit_fly': 'Drosophila melanogaster',
    'roundworm': 'Caenorhabditis elegans',
    'zebrafish': 'Danio rerio'
}

def to_sci_name(organism: str) -> str:
    """
    Convert organism name to scientific name.
    """
    try:
        return f'"{ORGANISMS[organism]}"'
    except KeyError:
        raise ValueError(f"Organism '{organism}' not found in list.")

## default time span limits
MAX_DATE = (datetime.now()).strftime('%Y/%m/%d')
MIN_DATE = (datetime.now() - timedelta(days=5 * 365)).strftime('%Y/%m/%d')

@tool 
def esearch_scrna(
    query_terms: Annotated[List[str], "Entrez query terms"]=["single cell RNA sequencing", "single cell RNA-seq"],
    database: Annotated[str, "Database name ('sra' or 'gds')"]="sra",
    organisms: Annotated[List[str], "List of organisms to search."]=["human", "mouse"],
    min_date: Annotated[str, "Minimum date to search back (%Y/%m/%d)."]=MIN_DATE,
    max_date: Annotated[str, "Maximum date to search back (%Y/%m/%d)."]=MAX_DATE,
    max_ids: Annotated[Optional[int], "Maximum number of IDs to return."]=10,
    )-> Annotated[List[str], "Entrez IDs of database records"]:
    """
    Find single cell RNA-seq datasets in the SRA or GEO databases.
    """
    set_entrez_access()
    esearch_query = ""

    # check if query terms are provided
    query_terms = " OR ".join([f'"{x}"' for x in query_terms])
    esearch_query += f"({query_terms})"
    
    # add date range
    date_range = f"{min_date}:{max_date}[PDAT]"
    esearch_query += f" AND ({date_range})"
    
    # add organism
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

    # debug model
    max_ids = 2 if os.getenv("DEBUG_MODE") == "true" else max_ids

    # return entrez IDs 
    return esearch_batch(esearch_query, database, max_ids, filter_existing=True)

def esearch_batch(
    esearch_query: str, database: str, max_ids: Optional[int], 
    verbose: bool=False,
    filter_existing: bool=False
    ) -> List[str]:
    # if filter existing, connect to database
    existing_ids = []
    if filter_existing:
        with db_connect() as conn:
            existing_ids = db_get_unprocessed_records(conn)["entrez_id"].tolist()
    # query 
    ids = []
    retstart = 0
    retmax = 10000
    while True:
        try:
            # search
            search_handle = Entrez.esearch(
                db=database, 
                term=esearch_query, 
                retstart=retstart, 
                retmax=retmax
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            # add IDs 
            ids.extend(
                [x for x in search_results['IdList'] if x not in existing_ids]
            )
            retstart += retmax
            time.sleep(0.34)
            if verbose:
                print(f"No. of IDs found: {len(ids)}", file=sys.stderr)
            if max_ids and len(ids) >= max_ids:
                break
            if retstart >= int(search_results['Count']):
                break
        except Exception as e:
            print(f"Error searching {database} with query: {esearch_query}: {str(e)}")
            break 
        
    # just unique IDs
    ids = list(set(ids))
        
    # return IDs
    if max_ids:
        ids = ids[:max_ids]
    return ids

@tool 
def esearch(
    esearch_query: Annotated[str, "Entrez query string."],
    database: Annotated[str, "Database name (e.g., sra, gds, or pubmed)"],
    )-> Annotated[List[str], "Entrez IDs of database records"]:
    """
    Run an Entrez search query and return the Entrez IDs of the results.
    Example query for single cell RNA-seq:
        `("single cell"[Title] OR "single-cell"[Title] OR "scRNA-seq"[Title])`
    Example query for an ENA accession number (database = sra):
        `ERX13336121`
    Example query for a GEO accession number (database = gds):
        `GSE51372`
    """
    set_entrez_access()

    # debug model
    max_records = 2 if os.getenv("DEBUG_MODE") == "TRUE" else None

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
    while True:
        try:
            search_handle = Entrez.esearch(
                db=database, 
                term=esearch_query, 
                retstart=retstart, 
                retmax=retmax
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            # delete unneeded keys
            to_rm = ["RetMax", "RetStart"]
            for key in to_rm:
                if key in search_results.keys():
                    del search_results[key]
            # add to records
            records.append(str(search_results))
            # update retstart
            retstart += retmax
            time.sleep(0.34)
            if max_records and len(records) >= max_records:
                break
            if retstart >= int(search_results['Count']):
                break
        except Exception as e:
            print(f"Error searching {database} with query: {esearch_query}: {str(e)}", file=sys.stderr)
            break 
        
    # return records
    if len(records) == 0:
        return f"No records found for query: {esearch_query}"
    if max_records:
        records = records[:max_records]  # debug
    return str(records)


if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # scRNA-seq
    #query = '("single cell RNA sequencing" OR "single cell RNA-seq")'
    #query = '("bulk RNA sequencing")'
    #input = {"esearch_query" : query, "database" : "sra", "previous_days" : 90}
    #input = {"esearch_query" : query, "database" : "gds", "previous_days" : 60}
    input = {}
    print(esearch_scrna.invoke(input))

    # esearch accession
    input = {"esearch_query" : "GSE51372", "database" : "sra"}
    input = {"esearch_query" : "GSE121737", "database" : "gds"}
    input = {"esearch_query" : "35447314", "database" : "sra"}
    input = {"esearch_query" : "35447314", "database" : "gds"}
    #print(esearch.invoke(input))

