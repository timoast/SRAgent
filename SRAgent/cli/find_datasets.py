# import
## batteries
import os
import sys
import asyncio
import argparse
from typing import List
## 3rd party
from Bio import Entrez
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
## package
from SRAgent.cli.utils import CustomFormatter
from SRAgent.workflows.find_datasets import create_find_datasets_graph
from SRAgent.agents.utils import create_step_summary_chain
from SRAgent.tools.utils import set_entrez_access

# functions
def find_datasets_parser(subparsers):
    help = 'Obtain datasets and process each via the srx-info workflow'
    desc = """
    Example message: "Obtain recent single cell RNA-seq datasets in the SRA database"
    """
    sub_parser = subparsers.add_parser(
        'find-datasets', help=help, description=desc, formatter_class=CustomFormatter
    )
    sub_parser.set_defaults(func=find_datasets_main)
    sub_parser.add_argument(
        'message', type=str, help='Message to instruct the agent. See the Description'
    ) 
    sub_parser.add_argument(
        '--max-datasets', type=int, default=10, help='Maximum number of datasets to analyze'
    )
    sub_parser.add_argument(
        '--no-summaries', action='store_true', default=False, help='No LLM summaries'
    )
    sub_parser.add_argument(
        '--max-concurrency', type=int, default=6, help='Maximum number of concurrent processes'
    )
    sub_parser.add_argument(
        '--recursion-limit', type=int, default=200, help='Maximum recursion limit'
    )
    sub_parser.add_argument(
        '-o', '--organisms', type=str, nargs='+', default=["human", "mouse"],
        choices=[
            "human", "mouse", "rat", "macaque", "marmoset", "horse", "dog", "bovine", "sheep", "pig", 
            "rabbit", "naked_mole_rat", "chimpanzee", 
            "chicken", "frog", "zebrafish", 
            "fruit_fly", "blood_fluke", "roundworm", "mosquito", 
            "thale_cress", "rice", "tomato", "corn"
        ],
        help='Organisms to search for'
    )
    sub_parser.add_argument(
        '--use-database', action='store_true', default=False, help='Use the SQL database'
    )

async def _find_datasets_main(args):
    """
    Main function for invoking the find-datasets workflow
    """
    # set email and api key
    set_entrez_access()
    
    # create supervisor agent
    graph = create_find_datasets_graph()
    if not args.no_summaries:
        step_summary_chain = create_step_summary_chain()

    # set graph inpu
    config = {
        "max_concurrency": args.max_concurrency,
        "recursion_limit": args.recursion_limit,
        "configurable": {
            "organisms": args.organisms,
            "max_datasets": args.max_datasets,
            "use_database": args.use_database
        }
    }
    input = {"messages" : [HumanMessage(content=args.message)]}

    # call the graph
    final_state = None
    i = 0
    async for step in graph.astream(input, config=config):
        i += 1
        final_state = step
        if args.no_summaries:
            nodes = ",".join(list(step.keys()))
            print(f"Step {i}: {nodes}")
        else:
            msg = await step_summary_chain.ainvoke({"step": step})
            print(f"Step {i}: {msg.content}")

    # print final results
    if final_state:
        print(f"\n#-- Final results --#")
        try:
            for msg in final_state["final_state_node"]["messages"]:
                try:
                    print(msg.content)
                except AttributeError:
                    if isinstance(msg, list):
                        for x in msg:
                            print(x)
                    else:
                        print(msg)
        except KeyError:
            print("Processing skipped")

def find_datasets_main(args):
    asyncio.run(_find_datasets_main(args))

# main
if __name__ == '__main__':
    pass