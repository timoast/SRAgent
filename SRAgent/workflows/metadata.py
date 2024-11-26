# import
## batteries
import os
import operator
from functools import partial
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, TypedDict, Sequence
## 3rd party
from pydantic import BaseModel, Field
from Bio import Entrez
import pandas as pd
from dotenv import load_dotenv
from langgraph.types import Send
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from langgraph.graph import START, END, StateGraph
## package
from SRAgent.agents.convert import create_convert_graph, invoke_convert_graph
from SRAgent.agents.metadata import create_metadata_graph, invoke_metadata_graph, get_metadata_items

# classes
class GraphState(TypedDict):
    """
    Shared state of the agents in the graph
    """
    # messages
    messages: Annotated[Sequence[BaseMessage], operator.add]
    # database
    database: str
    # dataset entrez ID
    entrez_id: str
    # accessions
    SRX: Annotated[List[str], operator.add]
    # is_illumina
    is_illumina: Annotated[List[str], operator.add]
    # is_single_cell
    is_single_cell: Annotated[List[str], operator.add]
    # is_paired_end
    is_paired_end: Annotated[List[str], operator.add]
    # is_10x
    is_10x: Annotated[List[str], operator.add]
    # organism
    organism: Annotated[List[str], operator.add]


# functions
def create_convert_graph_node():
    """
    Create a convert graph node
    """
    graph = create_convert_graph()
    def invoke_convert_graph_node(state: GraphState) -> Dict[str, Any]:
        entrez_id = state["entrez_id"]
        database = state["database"]
        message = "\n".join([
            f"Convert Entrez ID {entrez_id} to SRX (or ERX) accessions.",
            f"The Entrez ID is associated with the {database} database."
        ])
        input = {"messages": [HumanMessage(message)]}
        return graph.invoke(input)
    return invoke_convert_graph_node

def continue_to_metadata(state: GraphState) -> Dict[str, Any]:
    """
    Parallel invoke of the metadata graph
    """
    # format the prompt for the metadata graph
    prompt = "\n".join([
        "For the SRA accession {SRX_accession}, find the following information:",
        ] + list(get_metadata_items().values())
    )
    
    # submit each accession to the metadata graph
    responses = []
    for SRX_accession in state["SRX"]:
        input = {
            "database": state["database"],
            "entrez_id": state["entrez_id"],
            "SRX": SRX_accession,
            "messages": [HumanMessage(prompt.format(SRX_accession=SRX_accession))]
        }
        responses.append(Send("metadata_graph_node", input))
    return responses

def final_state(state: GraphState) -> Dict[str, Any]:
    """
    Return the final state of the graph
    """
    # filter to messages that 
    messages = []
    for msg in state["messages"]:
        try:
            msg = [msg.content]
        except AttributeError:
            msg = [x.content for x in msg]
        for x in msg:
            if x.startswith("# SRX accession: "):
                messages.append(x)
    # final message
    message = "\n".join(messages)
    return {
        "messages": [AIMessage(content=message)]
    }

def create_metadata_workflow(db_add:bool = True):
    # metadata subgraph
    invoke_metadata_graph_p = partial(
        invoke_metadata_graph,
        graph=create_metadata_graph(db_add=db_add),
        to_return=["messages"]
    )

    #-- graph --#
    workflow = StateGraph(GraphState)

    # nodes
    workflow.add_node("convert_graph_node", create_convert_graph_node())
    workflow.add_node("metadata_graph_node", invoke_metadata_graph_p)
    workflow.add_node("final_state_node", final_state)

    # edges
    workflow.add_edge(START, "convert_graph_node")
    workflow.add_conditional_edges("convert_graph_node", continue_to_metadata, ["metadata_graph_node"])
    workflow.add_edge("metadata_graph_node", "final_state_node")
    workflow.add_edge("final_state_node", END)

    # compile the graph
    graph = workflow.compile()
    return graph