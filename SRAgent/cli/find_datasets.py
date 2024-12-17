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
    sub_parser.add_argument('message', type=str,
                            help='Message to instruct the agent. See the Description') 
    sub_parser.add_argument('--no-summaries', action='store_true', default=False,
                            help='No LLM summaries')
    sub_parser.add_argument('--max-concurrency', type=int, default=6, 
                            help='Maximum number of concurrent processes')
    sub_parser.add_argument('--recursion-limit', type=int, default=200,
                            help='Maximum recursion limit')


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

    # invoke agent
    config = {
        "max_concurrency": args.max_concurrency,
        "recursion_limit": args.recursion_limit
    }

    #
    input = {"messages" : [HumanMessage(content=args.message)]}
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

    if final_state:
        print(f"\n#-- Final results --#")
        try:
            for msg in final_state["final_state_node"]["messages"]:
                print(msg); 
        except KeyError:
            print("Processing skipped")

def find_datasets_main(args):
    asyncio.run(_find_datasets_main(args))

# main
if __name__ == '__main__':
    pass