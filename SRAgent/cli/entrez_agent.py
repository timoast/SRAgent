# import
## batteries
import os
import sys
import argparse
from Bio import Entrez
from SRAgent.cli.utils import CustomFormatter
from SRAgent.agents.supervisors import create_supervisor_agent, create_step_summary_chain, invoke_entrez_agent

# functions
def entrez_agent_parser(subparsers):
    help = 'Entrez Agent: general agent for working with Entrez databases.'
    desc = """
    # Example prompts:
    1. "Convert GSE121737 to SRX accessions"
    """
    sub_parser = subparsers.add_parser(
        'entrez-agent', help=help, description=desc, formatter_class=CustomFormatter
    )
    sub_parser.set_defaults(func=entrez_agent_main)
    sub_parser.add_argument('prompt', type=str, help='Prompt for the agent')    


def entrez_agent_main(args):
    """
    Main function for invoking the entrez agent
    """
    Entrez.email = os.getenv("EMAIL")

    # create supervisor agent
    agent = create_supervisor_agent()
    step_summary_chain = create_step_summary_chain()

    # invoke agent
    input = {"messages": [("user", args.prompt)]}
    invoke_entrez_agent(input, agent, step_summary_chain)

# main
if __name__ == '__main__':
    pass