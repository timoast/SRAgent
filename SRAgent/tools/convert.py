# langchain custom tools
import os
import sys
import time
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any
from dotenv import load_dotenv
from Bio import Entrez
from langchain_core.tools import tool

# functions
@tool 
def geo2sra(
    entrez_ids: Annotated[List[str], "List of Entrez IDs for the GEO database"],
    )-> Annotated[List[str], "List linked Entrez IDs for SRA records"]:
    """
    Convert a GEO Entrez ID to SRA Entrez IDs.
    """
    sra_ids = []
    for entrez_id in entrez_ids:
        # Fetch detailed GEO record to get links to SRA
        handle = Entrez.elink(dbfrom="gds", db="sra", id=entrez_id)
        links = Entrez.read(handle)
        handle.close()
        
        if links[0]['LinkSetDb']:
            sra_ids += [link['Id'] for link in links[0]['LinkSetDb'][0]['Link']]
        time.sleep(0.34) 

    # debug mode
    if os.getenv("DEBUG_MODE") == "TRUE":
        sra_ids = sra_ids[:2]

    # return SRA IDs
    return sra_ids



# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # test
    ## GEO IDs
    entrez_ids = {"entrez_ids" : ["200274955", "200120926"]}
    #print(geo2sra(entrez_ids))