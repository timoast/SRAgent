#!/usr/bin/env python
# import
import os
import argparse
from datetime import datetime, timedelta
from typing import List, Dict, Any, Tuple, Annotated
from dotenv import load_dotenv
import pandas as pd
from Bio import Entrez
from SRAgent.tools.esearch import esearch_batch
from SRAgent.cli.utils import CustomFormatter
from SRAgent.workflows.SRX_info import create_SRX_info_graph
from SRAgent.record_db import db_connect, db_get_srx_records


# global variables
## organism names
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
## default time span limits
MAX_DATE = datetime.now()
MIN_DATE = datetime.now() - timedelta(days=5 * 365) 

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def parse_args():
    desc = 'Obtain accessions and call metadata agent workflow'
    epi = """DESCRIPTION:
Run Entrez esearch to find dataset Entrez IDs and provide the IDs
to the Metadata Agent Workflow.
    """
    parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                     formatter_class=CustomFormatter)
    parser.add_argument('--query-terms', type=str, nargs='+',
                        default=['single cell RNA sequencing', 'single cell RNA-seq'],
                        help='esearch query terms')
    parser.add_argument('--database', type=str, default='sra',
                        choices=['sra', 'gds'], 
                        help='Database name (sra or gds).')
    parser.add_argument('--organisms', type=str, nargs='*', 
                        default=['human', 'mouse'],
                        choices=ORGANISMS.keys(),
                        help='Database name (sra or gds).')
    parser.add_argument('--min-date', type=str, default=MIN_DATE.strftime('%Y/%m/%d'),
                        help='Minimum date to search back.')
    parser.add_argument('--max-date', type=str, default=MAX_DATE.strftime('%Y/%m/%d'),
                        help='Maximum date to search back.')
    parser.add_argument('--email', type=str, default=os.getenv('EMAIL'),
                        help='Email address for Entrez')
    parser.add_argument('--max-concurrency', type=int, default=3, 
                        help='Maximum number of concurrent processes')
    parser.add_argument('--recursion-limit', type=int, default=100,
                        help='Maximum recursion limit')
    return parser.parse_args()

def to_sci_name(organism: str) -> str:
    """
    Convert organism name to scientific name.
    """
    try:
        return f'"{ORGANISMS[organism]}"'
    except KeyError:
        raise ValueError(f"Organism '{organism}' not found in list.")

def esearch(
    query_terms: Annotated[List[str], "Entrez query terms"],
    database: Annotated[str, "Database name ('sra' or 'gds')"],
    organisms: Annotated[List[str], "List of organisms to search."],
    min_date: Annotated[str, "Minimum date to search back."],
    max_date: Annotated[str, "Maximum date to search back."],
    to_exclude: Annotated[List[str], "Entrez IDs to exclude."] = [],
    )-> Annotated[List[str], "Entrez IDs of database records"]:
    """
    Perform esearch to find database records.
    """
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
    #esearch_query += ' AND "study"'

    # add exclusions
    exclusions = " OR ".join([f"{x}[UID]" for x in to_exclude])
    if exclusions:
        esearch_query += f" NOT ({exclusions})"

    # debug model
    max_ids = 2 if os.getenv("DEBUG_MODE") == "true" else None

    # return entrez IDs
    return esearch_batch(esearch_query, database, max_ids, verbose=True)

def filter_existing_entrez_ids(entrez_ids: List[str]) -> List[str]:
    """
    Filter entrez IDs by checking if they are already in the database.
    """
    # get exclusions
    existing_entrez_ids = set()
    with db_connect() as conn:
        existing_entrez_ids = set(
            db_get_srx_records(conn, column="entrez_id", database=args.database)
        )

    # filter
    return [x for x in entrez_ids if x not in existing_entrez_ids]

def main(args):
    # set email
    if args.email:
        Entrez.email = args.email
    # set API key
    if 'NCBI_API_KEY' in os.environ:
        Entrez.api_key = os.environ['NCBI_API_KEY']

    # get entrez IDs for novel datasets
    entrez_ids = esearch(
        args.query_terms, 
        args.database, 
        args.organisms, 
        args.min_date, 
        args.max_date,
    )
    print(f"No. of entrez IDs (all): {len(entrez_ids)}")

    # filter 
    entrez_ids = filter_existing_entrez_ids(entrez_ids)
    print(f"No. of entrez IDs (novel): {len(entrez_ids)}")
    
    # create supervisor agent
    graph = create_SRX_info_graph()

    # invoke agent
    config = {
        "max_concurrency" : args.max_concurrency,
        "recursion_limit": args.recursion_limit
    }
    for entrez_id in entrez_ids:
        print(f"#-- Entrez ID: {entrez_id} (database={args.database}) --#")
        input = {"entrez_id": entrez_id, "database": args.database}
        # stream invoke graph
        final_state = None
        for i,step in enumerate(graph.stream(input, config=config)):
            final_state = step
            print(f"Step {i+1}: {step}")
        # print final state
        if final_state:
            print("")
            try:
                print(final_state["final_state_node"]["messages"][-1].content)
            except KeyError:
                print("No final state message.")
        print("")

## script main
if __name__ == '__main__':
    load_dotenv()
    args = parse_args()
    main(args)