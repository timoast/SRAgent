# import
## batteries
import os
import asyncio
from Bio import Entrez
from langchain_core.messages import HumanMessage
from SRAgent.cli.utils import CustomFormatter
from SRAgent.workflows.tissue_ontology import create_tissue_ontology_workflow
from SRAgent.agents.utils import create_agent_stream

# functions
def tissue_ontology_parser(subparsers):
    help = 'Tissue Ontology: categorize tissue descriptions using the Uberon ontology.'
    desc = """
    # Example prompts:
    1. "Categorize the following tissue: brain"
    2. "What is the Uberon ID for hippocampus?"
    3. "Tissues: lung, heart, liver"
    4. "Find the ontology term for the thin layer of epithelial cells lining the alveoli in lungs"
    5. "What is the Uberon classification for skeletal muscle tissue?"
    """
    sub_parser = subparsers.add_parser(
        'tissue-ontology', help=help, description=desc, formatter_class=CustomFormatter
    )
    sub_parser.set_defaults(func=tissue_ontology_main)
    sub_parser.add_argument('prompt', type=str, help='Tissue description(s) to categorize')    
    sub_parser.add_argument('--no-summaries', action='store_true', default=False,
                            help='No LLM summaries')
    sub_parser.add_argument('--max-concurrency', type=int, default=3, 
                            help='Maximum number of concurrent processes')
    sub_parser.add_argument('--recursion-limit', type=int, default=40,
                            help='Maximum recursion limit')
    
def tissue_ontology_main(args):
    """
    Main function for invoking the tissue ontology workflow
    """
    # set email and api key
    Entrez.email = os.getenv("EMAIL")
    Entrez.api_key = os.getenv("NCBI_API_KEY")

    # create workflow
    workflow = create_tissue_ontology_workflow()

    # create config
    config = {
        "max_concurrency": args.max_concurrency,
        "recursion_limit": args.recursion_limit
    }

    # prepare input
    input = {"messages": [HumanMessage(content=args.prompt)]}
    
    # invoke workflow
    results = asyncio.run(workflow.ainvoke(input, config=config))
    
    # print results
    if results:
        print("\nUberon IDs for the tissue descriptions:")
        for uberon_id in results:
            print(f"- {uberon_id}")
    else:
        print("No suitable Uberon ontology terms found for the provided tissue descriptions.")

# main
if __name__ == '__main__':
    pass
