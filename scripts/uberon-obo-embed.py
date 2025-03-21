import os
import sys
from typing import List

import obonet
import chromadb
import networkx as nx
from langchain_core.documents import Document
from langchain_openai import OpenAIEmbeddings
from langchain_chroma import Chroma


def parse_cli_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments containing GCP path, feature type, and number of workers.
    """
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description='Convert OBO file to Chroma vector store',
        epilog="""Example:
    ./uberon-obo-embed.py /path/to/uberon.obo
    """,
        formatter_class=CustomFormatter
    )
    parser.add_argument(
        'obo_path', type=str, help='Path to the OBO file'
    )
    parser.add_argument(
        '--output-db-path', type=str, default='uberon_chroma',
        help='Path to the output Chroma database'
    )
    parser.add_argument(
        '--collection-name', type=str, default='uberon',
        help='Name of the collection in the Chroma database'
    )
    parser.add_argument(
        '--model', type=str, default='text-embedding-3-small',
        help='Model to use for embeddings'
    )
    parser.add_argument(
        '--max-embeddings', type=int, default=None,
        help='Maximum number of embeddings to generate'
    )
    return parser.parse_args()

def extract_definitions(
    graph: nx.MultiDiGraph,
    max_embeddings: int | None = None
) -> List[Document]:
    """
    Extracts definition texts, metadata, and node IDs from the ontology graph.
    
    Returns:
        A list of definition texts, a list of metadata dictionaries, and a list of node IDs.
    """
    documents = []
    for node_id, data in graph.nodes(data=True):
        definition = data.get("def")
        if not definition:
            continue
        documents.append(
            Document(
                page_content=definition,
                metadata={
                    "id": node_id,
                    "name": data.get("name", ""),
                }
            )
        )
        if max_embeddings and len(documents) >= max_embeddings:
            break
    return documents

def main(args: argparse.Namespace) -> None:
    # Read the OBO file into a graph.
    graph = obonet.read_obo(args.obo_path)
    documents = extract_definitions(graph, args.max_embeddings)

    if not documents:
        print("No definitions found in the OBO file.")
        return

    # Use LangChain's OpenAIEmbeddings.
    embeddings = OpenAIEmbeddings(model=args.model)  

    # Set output based on output-db-path
    os.makedirs(args.output_db_path, exist_ok=True) 

    # Create (or load) a Chroma vector store from the extracted texts.
    persistent_client = chromadb.PersistentClient(path=args.output_db_path)
    vector_store = Chroma(
        client=persistent_client,
        collection_name=args.collection_name,
        embedding_function=embeddings,  
    )

    # Add the documents to the vector store.
    vector_store.add_documents(documents=documents)

    # Verify documents were added
    collection = persistent_client.get_collection(args.collection_name)
    count = collection.count()
    print(f"Stored {len(documents)} embeddings in collection '{args.collection_name}' at '{args.output_db_path}'.")
    print(f"Collection now contains {count} total documents.")

if __name__ == "__main__":
    args = parse_cli_args()
    main(args)