# import
## batteries
import os
import re
import requests
from typing import Annotated
## 3rd party
from langchain_core.tools import tool
from langchain_openai import OpenAIEmbeddings
from langchain_chroma import Chroma
import chromadb
import obonet
## package


# functions
def verify_collection(persistent_client: chromadb.PersistentClient, collection_name: str) -> None:
    """
    Verify that the collection exists and has documents
    """
    try:
        collection = persistent_client.get_collection(collection_name)
        count = collection.count()
        print(f"Found {count} documents in collection '{collection_name}'")
    except Exception as e:
        msg = f"Error accessing collection: {e}"
        msg += f"\nAvailable collections: {persistent_client.list_collections()}"
        raise Exception(msg)

def load_vector_store(chroma_path: str) -> Chroma:
    """Load a Chroma vector store from the specified path."""
    # Ensure the path exists
    if not os.path.exists(chroma_path):
        raise FileNotFoundError(f"Chroma DB directory not found: {chroma_path}")
    
    # Use the base name of the path as the collection name
    collection_name = os.path.basename(chroma_path.rstrip("/"))

    # Initialize embeddings
    embeddings = OpenAIEmbeddings(model="text-embedding-3-small")
    
    # Load the persistent Chroma client
    persistent_client = chromadb.PersistentClient(path=chroma_path)
    
    # Load the existing vector store
    vector_store = Chroma(
        client=persistent_client,
        collection_name=collection_name,
        embedding_function=embeddings,
    )
    return vector_store

@tool
def query_vector_db(
    query: Annotated[str, "The semantic search query"],
    k: Annotated[int, "The number of results to return"]=3,
    #db_path: Annotated[str, "The path to the Chroma DB"] = 
    ) -> str: 
    """
    Perform a semantic search by querying a vector store 
    """
    # DEBUG: 
    db_path = "/home/nickyoungblut/dev/python/SRAgent/tmp/tissue_ontology/uberon-simple_chroma"
    # load the vector store
    vector_store = load_vector_store(db_path)
    # query the vector store
    message = ""
    try:
        results = vector_store.similarity_search(query, k=k)
        message += f"# Results for query: \"{query}\"\n"
        for i, res in enumerate(results, 1):
            id = res.metadata.get("id", "No ID available")
            if not id:
                continue
            message += f"{i}. {id}\n"
            name = res.metadata.get("name", "No name available")
            message += f"  Ontology name: {name}\n"
            message += f"  Description: {res.page_content}\n"
    except Exception as e:
        return f"Error performing search: {e}"
    if not message:
        message = f"No results found for query: \"{query}\". Consider refining your query."
    return message

@tool 
def get_neighbors(
    uberon_id: Annotated[str, "The Uberon ID (UBERON:XXXXXXX)"],
    ) -> str: 
    """
    Get the neighbors of a given Uberon ID in the tissue ontology.
    """
    # check the ID format
    if not re.match(r"UBERON:\d{7}", uberon_id):
        return f"Invalid Uberon ID format: \"{uberon_id}\". The format must be \"UBERON:XXXXXXX\"."

    # DEBUG:
    obo_path = "/home/nickyoungblut/dev/python/SRAgent/tmp/tissue_ontology/uberon-simple.obo"

    # read the ontology graph
    g = obonet.read_obo(obo_path)

    # get neighbors
    message = ""
    try:
        message += f"# Neighbors in the ontology for: \"{uberon_id}\"\n"
        for i,node_id in enumerate(g.neighbors(uberon_id), 1):
            if not g.nodes[node_id]:
                continue
            node_name = g.nodes[node_id]["name"]
            node_def = g.nodes[node_id]["def"]
            message += f"{i}. {node_id}\n"
            message += f"  Ontology name: {node_name}\n"
            message += f"  Description: {node_def}\n"
    except Exception as e:
        return f"Error getting neighbors: {e}"

    if not message:
        message = f"No neighbors found for ID: \"{uberon_id}\"."
    return message

@tool
def query_uberon_ols(
    search_term: Annotated[str, "The term to search for in the Uberon ontology"]
) -> str:
    """
    Query the Ontology Lookup Service (OLS) for Uberon terms matching the search term.
    """
    url = f"https://www.ebi.ac.uk/ols/api/search?q={search_term}&ontology=uberon"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
    except Exception as e:
        return f"Error querying OLS API: {e}"

    results = data.get("response", {}).get("docs", [])
    if not results:
        return f"No results found for search term: '{search_term}'."

    message = f"# Results from OLS for '{search_term}':\n"
    for i, doc in enumerate(results, 1):
        # Each doc should have an 'obo_id', a 'label', and possibly a 'description'
        obo_id = doc.get("obo_id", "No ID")
        label = doc.get("label", "No label")
        description = doc.get("description", ["No description"])[0]  # description is usually a list
        message += f"{i}. {obo_id} - {label}\n   Description: {description}\n"
    return message


if __name__ == "__main__":
    from dotenv import load_dotenv
    load_dotenv(override=True)

    # semantic search
    #query = "brain"
    #results = query_vector_db.invoke({"query" : query})
    #print(results)

    # get neighbors
    # input = {'uberon_id': 'UBERON:0000010'}
    # #input = {'uberon_id': 'UBERON:0000013'}
    # neighbors = get_neighbors.invoke(input)
    # print(neighbors)

    # query OLS
    input = {'search_term': 'brain'}
    results = query_uberon_ols.invoke(input)
    print(results)
    

    
