# import
## batteries
import os
import sys
import argparse
from typing import List
from Bio import Entrez
from SRAgent.cli.utils import CustomFormatter
from SRAgent.workflows.metadata import create_metadata_workflow
from SRAgent.agents.utils import create_step_summary_chain
from SRAgent.utils import filter_by_db

# functions
def metadata_agent_parser(subparsers):
    help = 'Metadata Agent: Obtain metadata for SRA experiments.'
    desc = """
    """
    sub_parser = subparsers.add_parser(
        'metadata', help=help, description=desc, formatter_class=CustomFormatter
    )
    sub_parser.set_defaults(func=metadata_agent_main)
    sub_parser.add_argument('entrez_ids', type=str, nargs='+',
                            help='>=1 Entrez IDs to query')    
    sub_parser.add_argument('--database', type=str, default='sra',
                            choices=['gds', 'sra'],
                            help='Entrez database origin of the Entrez IDs')
    sub_parser.add_argument('--no-filter', action='store_true', default=False,
                            help='Do not filter Entrez IDs already in the metadata database')
    sub_parser.add_argument('--no-summaries', action='store_true', default=False,
                            help='No LLM summaries')
    sub_parser.add_argument('--max-concurrency', type=int, default=3, 
                            help='Maximum number of concurrent processes')
    sub_parser.add_argument('--recursion-limit', type=int, default=200,
                            help='Maximum recursion limit')

def metadata_agent_main(args):
    """
    Main function for invoking the entrez agent
    """
    # filter entrez_ids
    if not args.no_filter:
        args.entrez_ids = filter_by_db(args.entrez_ids)
        if len(args.entrez_ids) == 0:
            print("All Entrez IDs are already in the metadata database.", file=sys.stderr)
            return 0

    # set email and api key
    Entrez.email = os.getenv("EMAIL")
    Entrez.api_key = os.getenv("NCBI_API_KEY")

    # create supervisor agent
    graph = create_metadata_workflow()
    step_summary_chain = create_step_summary_chain()

    # invoke agent
    config = {
        "max_concurrency" : args.max_concurrency,
        "recursion_limit": args.recursion_limit
    }
    for entrez_id in args.entrez_ids:
        print(f"#-- Entrez ID: {entrez_id} (database={args.database}) --#")
        input = {"entrez_id": entrez_id, "database": args.database}
        # stream invoke graph
        final_state = None
        for i,step in enumerate(graph.stream(input, config={"max_concurrency" : 3, "recursion_limit": 200})):
            final_state = step
            if args.no_summaries:
                print(f"Step {i+1}: {step}")
            else:
                msg = step_summary_chain.invoke({"step": step})
                print(f"Step {i+1}: {msg.content}")
        # print final state
        if final_state:
            try:
                print(final_state["final_state_node"]["messages"][-1].content)
            except KeyError:
                print("No final state message.")
        print("")

# main
if __name__ == '__main__':
    pass