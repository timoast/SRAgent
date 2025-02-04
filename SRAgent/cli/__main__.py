#!/usr/bin/env python
# libraries
## batteries
import os
import sys
import argparse
## 3rd party
from dotenv import load_dotenv
## package
from SRAgent.cli.utils import CustomFormatter
from SRAgent.cli.entrez import entrez_agent_parser, entrez_agent_main
from SRAgent.cli.sragent import sragent_parser, sragent_main
from SRAgent.cli.srx_info import SRX_info_agent_parser, SRX_info_agent_main
from SRAgent.cli.find_datasets import find_datasets_parser, find_datasets_main


# functions
def arg_parse(args=None) -> dict:
    """
    Parse command line arguments.
    """
    desc = "SRAgent: A multi-agent tool for working with the SRA"
    epi = """DESCRIPTION:
    SRAgent is a multi-agent tool for working with the Sequence Read Archive (SRA) database
    and other Entriz databases. It is designed to be a flexible and easy-to-use tool for
    interacting with the SRA and Entrez.
    """
    # check for OP
    if os.getenv("OPENAI_API_KEY") is None:
        raise ValueError("OPENAI_API_KEY not found in environment")
    
    # main parser
    parser = argparse.ArgumentParser(
        description=desc,
        epilog=epi,
        formatter_class=CustomFormatter
    )
    parser.add_argument(
        '--tenant', type=str, default="prod", choices = ["test", "prod"],
        help='Database tenant to connect to. Defaults to DYNACONF env variable'
    )

    # subparsers
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")
    ## Entrez agent
    entrez_agent_parser(subparsers)
    ## SR agent
    sragent_parser(subparsers)
    ## SRX info agent
    SRX_info_agent_parser(subparsers)
    ## Find datasets
    find_datasets_parser(subparsers)
    # parsing args
    return parser.parse_args()

def main():
    # load environment variables
    load_dotenv(override=True)
    # parsing args
    args = arg_parse()
    # set database tenant
    if args.tenant:
        os.environ["DYNACONF"] = args.tenant
    
    # which subcommand
    if not args.command:
        print("Provide a subcommand or use -h/--help for help")
        sys.exit(0)
    elif args.command == "entrez":
        entrez_agent_main(args)
    elif args.command == "sragent":
        sragent_main(args)
    elif args.command.lower() == "srx-info":
        SRX_info_agent_main(args)
    elif args.command.lower() == "find-datasets":
        find_datasets_main(args)
    else:
        print("No command specified. Exiting ...")
        sys.exit(0)

    
if __name__ == "__main__":
    main()