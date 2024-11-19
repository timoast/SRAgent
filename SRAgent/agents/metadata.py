# import 
import os
import operator
from functools import partial
from enum import Enum
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Sequence, TypedDict
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from pydantic import BaseModel, Field
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langgraph.graph import START, END, StateGraph, MessagesState
## package
from SRAgent.agents.entrez import create_entrez_agent

# classes
class YesNo(Enum):
    """
    Yes or no
    """
    YES = "yes"
    NO = "no"
    UNSURE = "unsure"

class Organism(Enum):
    """
    Organism
    """
    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"
    MONKEY = "monkey"
    HORSE = "horse"
    DOG = "dog"
    BOVINE = "bovine"
    CHICKEN = "chicken"
    SHEEP = "sheep"
    PIG = "pig"
    FRUIT_FLY = "fruit_fly"
    ROUNDWORM = "roundworm"
    ZEBRAFISH = "zebrafish"
    OTHER = "other"

class Metadata(BaseModel):
    is_illumina: List[YesNo] 
    is_single_cell: List[YesNo]
    is_paired_end: List[YesNo]
    is_10x: List[YesNo]
    organism: List[Organism]

class Choices(Enum):
    CONTINUE = "Continue"
    STOP = "Stop"

class Choice(BaseModel):
    Choice: Choices

class GraphState(TypedDict):
    """
    Shared state of the agents in the graph
    """
    messages: Annotated[Sequence[BaseMessage], operator.add, "Messages"]
    is_illumina: Annotated[List[YesNo], operator.add, "Is Illumina sequence data?"]
    is_single_cell: Annotated[List[YesNo], operator.add, "Is single cell RNA-seq data?"]
    is_paired_end: Annotated[List[YesNo], operator.add, "Is paired-end sequencing data?"]
    is_10x: Annotated[List[YesNo], operator.add, "Is 10X Genomics data?"]
    organism: Annotated[List[Organism], operator.add, "Organism sequenced"]
    route: Annotated[str, "Route"]
    round: Annotated[int, "Round"]


# functions
def get_metadata_items():
    metadata_items = [
        " - Is the study Illumina sequence data?",
        " - Is the study Illumina sequence data?",
        " - Is the study single cell RNA-seq data?",
        " - Is the study paired-end sequencing data?",
        " - Is the study 10X Genomics data?",
        " - Which organism was sequenced?"
    ]
    return metadata_items

def invoke_entrez_agent_node(state: GraphState):
    entrez_agent = create_entrez_agent()
    response = entrez_agent.invoke({"messages" : state["messages"]})
    return {"messages" : [response["messages"][-1]]}

def create_get_metadata_node():
    model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)

    def invoke_get_metadata_node(state: GraphState):
        """
        Structured data extraction
        """
        message = state["messages"][-1].content
        prompt = "\n".join([
            "Your job is to extract metadata from the provided text on a Sequence Read Archive (SRA) experiment.",
            "If there is not enough information to determine the metadata, please respond with 'unsure'.",
            "The specific metadata to extract:"] + get_metadata_items() + [
            "\n",
            "The provided text:",
            message
        ])
        model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)
        response = model.with_structured_output(Metadata, strict=True).invoke(prompt)
        return {
            "is_illumina" : response.is_illumina,
            "is_single_cell" : response.is_single_cell,
            "is_paired_end" : response.is_paired_end,
            "is_10x" : response.is_10x,
            "organism" : response.organism
        }
    
    return invoke_get_metadata_node

def create_router_node():
    model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)

    def invoke_router_node(state: GraphState):
        """
        Router for the graph
        """
        # create prompt
        prompt = "\n".join([
            "You determine whether all of appropriate metadata has been extracted.",
            "The metadata of interest:"
            ] + get_metadata_items() + [
            "\n",
            "If most or all of the metadata is \"unsure\", then the task is incomplete."
            ]
        )

        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", prompt),
            ("system", "Here are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
            # Add the final question/instruction
            ("human", "Based on the messages above, select STOP if the task is complete or CONTINUE if more information is needed."),
        ])
        formatted_prompt = prompt.format_messages(history=state["messages"])
        # call the model
        response = model.with_structured_output(Choice, strict=True).invoke(formatted_prompt)
        # format the response
        if response.Choice.value == Choices.CONTINUE.value:
            message = "At least some of the meadata is still uncertain. Please try to provide more information."
        else:
            message = "Enough of the metadata has been extracted."
        return {"route": response.Choice.value, "rounds": 1, "messages": [AIMessage(content=message)]}
    
    return invoke_router_node

def route_interpret(state: GraphState) -> str:
    """
    Determine the route based on the current state of the conversation.
    """
    if state["rounds"] >= 2:
        return END
    return "entrez_agent_node" if state["route"] == "Continue" else END

def create_metadata_graph():
    #-- graph --#
    workflow = StateGraph(GraphState)

    # nodes
    workflow.add_node("entrez_agent_node", invoke_entrez_agent_node)
    workflow.add_node("get_metadata_node", create_get_metadata_node())
    workflow.add_node("router_node", create_router_node())

    # edges
    workflow.add_edge(START, "entrez_agent_node")
    workflow.add_edge("entrez_agent_node", "get_metadata_node")
    workflow.add_edge("get_metadata_node", "router_node")
    workflow.add_conditional_edges("router_node", route_interpret)

    # compile the graph
    graph = workflow.compile()
    return graph


# main
if __name__ == "__main__":
    from functools import partial
    from Bio import Entrez

    #-- setup --#
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    #-- graph --#
    SRX_accession = "SRX11740066"
    prompt = "\n".join([
        f"For the SRA accession {SRX_accession}, find the following information:",
        ] + get_metadata_items()
    )
    input = {"messages" : [HumanMessage(content=prompt)]}
    graph = create_metadata_graph()
    for step in graph.stream(input, config={"max_concurrency" : 3, "recursion_limit": 30}):
        print(step)

    ## invoke with graph object directly provided
    #invoke_metadata_graph = partial(invoke_metadata_graph, graph=graph)
    #print(invoke_metadata_graph(input))