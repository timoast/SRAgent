# import 
import os
import operator
from functools import partial
from enum import Enum
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Sequence, TypedDict
import gspread
from gspread_dataframe import set_with_dataframe
import pandas as pd
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
    Yes, no, or unsure
    """
    YES = "yes"
    NO = "no"
    UNSURE = "unsure"

class OrganismEnum(Enum):
    """
    Organism sequenced
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

class MetadataEnum(BaseModel):
    """
    Metadata to extract
    """
    is_illumina: YesNo
    is_single_cell: YesNo
    is_paired_end: YesNo
    is_10x: YesNo
    organism: OrganismEnum

class ChoicesEnum(Enum):
    CONTINUE = "Continue"
    STOP = "Stop"

class Choice(BaseModel):
    Choice: ChoicesEnum

class GraphState(TypedDict):
    """
    Shared state of the agents in the graph
    """
    database: Annotated[str, "Database"]
    entrez_id: Annotated[str, "Entrez ID"]
    SRX: Annotated[str, "SRX accession to query"]
    is_illumina: Annotated[str, "Is Illumina sequence data?"]
    is_single_cell: Annotated[str, "Is single cell RNA-seq data?"]
    is_paired_end: Annotated[str, "Is paired-end sequencing data?"]
    is_10x: Annotated[str, "Is 10X Genomics data?"]
    organism: Annotated[str, "Organism sequenced"]
    messages: Annotated[Sequence[BaseMessage], operator.add]
    route: Annotated[str, "Route"]
    rounds: Annotated[int, operator.add]
    

# functions
def get_metadata_items():
    metadata_items = [
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
        # format prompt
        prompt = "\n".join([
            "Your job is to extract metadata from the provided text on a Sequence Read Archive (SRA) experiment.",
            "If there is not enough information to determine the metadata, please respond with 'unsure'.",
            "The specific metadata to extract:"] + get_metadata_items())
        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt),
            ("system", "\nHere are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
        ])
        prompt = prompt.format_messages(history=state["messages"])

        # call the model
        model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)
        response = model.with_structured_output(MetadataEnum, strict=True).invoke(prompt)
        return {
            "is_illumina" : response.is_illumina.value,
            "is_single_cell" : response.is_single_cell.value,
            "is_paired_end" : response.is_paired_end.value,
            "is_10x" : response.is_10x.value,
            "organism" : response.organism.value
        }
    
    return invoke_get_metadata_node

def create_router_node():
    model = ChatOpenAI(model_name="gpt-4o-mini", temperature=0)

    def invoke_router_node(state: GraphState):
        """
        Router for the graph
        """
        # create prompt
        prompt1 = "\n".join([
            "You determine whether all of appropriate metadata has been extracted.",
            "If most or all of the metadata is \"unsure\", then the task is incomplete."
            ]
        )
        prompt2 = "\n".join([
            "\nThe extracted metadata:",
            " - Is the study Illumina sequence data?",
            "    - " + state['is_illumina'],
            " - Is the study single cell RNA-seq data?",
            "    - " + state['is_single_cell'],
            " - Is the study paired-end sequencing data?",
            "    - " + state['is_paired_end'],
            " - Is the study 10X Genomics data?",
            "    - " + state['is_10x'],
            " - Which organism was sequenced?",
            "    - " + state['organism']
        ])

        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", prompt1),
            ("system", "\nHere are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
            ("system", prompt2),
            # Add the final question/instruction
            ("human", "Based on the information above, select STOP if the task is complete or CONTINUE if more information is needed."),
        ])
        formatted_prompt = prompt.format_messages(history=state["messages"])
        # call the model
        response = model.with_structured_output(Choice, strict=True).invoke(formatted_prompt)
        # format the response
        if response.Choice.value == ChoicesEnum.CONTINUE.value:
            message = "At least some of the metadata is still uncertain. Please try to provide more information by using a different approach (e.g., different tool calls)."
        else:
            message = "Enough of the metadata has been extracted."
        return {"route": response.Choice.value, "rounds": 1, "messages": [AIMessage(content=message)]}
    
    return invoke_router_node

def route_interpret(state: GraphState) -> str:
    """
    Determine the route based on the current state of the conversation.
    """
    if state["rounds"] >= 2:
        return "add2db_node"
    return "entrez_agent_node" if state["route"] == "Continue" else "add2db_node"

def fmt(x):
    if type(x) != list:
        return x
    return ";".join([str(y) for y in x])

def add2db(state: GraphState):
    """
    Add results to the database
    """
    # create dataframe from state
    results = pd.DataFrame([{
        "database": state["database"],
        "entrez_id": state["entrez_id"],
        "SRX": fmt(state["SRX"]),
        "is_illumina": state["is_illumina"],
        "is_single_cell": state["is_single_cell"],
        "is_paired_end": state["is_paired_end"],
        "is_10x": state["is_10x"],
        "organism": state["organism"]
    }])

    # Authenticate and open the Google Sheet
    db_name = "SRAgent_database"
    gc = gspread.service_account(filename=os.getenv("GOOGLE_APPLICATION_CREDENTIALS"))
    sheet = gc.open(db_name)
    worksheet = sheet.worksheet("database")
    # Read existing data to find where to append
    existing_data = pd.DataFrame(worksheet.get_all_values())
    next_row = len(existing_data) + 1
    # Append the data starting from the next row
    set_with_dataframe(worksheet, results, row=next_row, col=1, include_index=False, include_column_header=False)
    # Return
    return {"messages": [AIMessage(content=f"The results have been added to the database: {db_name}")]}

def create_metadata_graph():
    #-- graph --#
    workflow = StateGraph(GraphState)

    # nodes
    workflow.add_node("entrez_agent_node", invoke_entrez_agent_node)
    workflow.add_node("get_metadata_node", create_get_metadata_node())
    workflow.add_node("router_node", create_router_node())
    workflow.add_node("add2db_node", add2db)

    # edges
    workflow.add_edge(START, "entrez_agent_node")
    workflow.add_edge("entrez_agent_node", "get_metadata_node")
    workflow.add_edge("get_metadata_node", "router_node")
    workflow.add_conditional_edges("router_node", route_interpret)

    # compile the graph
    graph = workflow.compile()
    return graph

def invoke_metadata_graph(
    state: GraphState, 
    graph: StateGraph,
    to_return: List[str] = MetadataEnum.__fields__.keys()
) -> Annotated[dict, "Response from the graph"]:
    """
    Invoke the graph to convert Entrez IDs & non-SRA accessions to SRA accessions
    """
    response = graph.invoke(state)
    filtered_response = {key: [response[key]] for key in to_return}
    return filtered_response

# main
if __name__ == "__main__":
    from functools import partial
    from Bio import Entrez

    #-- setup --#
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    #-- graph --#
    #SRX_accession = "SRX25716878"
    #SRX_accession = "SRX20554853"
    #SRX_accession = "SRX20554856"
    SRX_accession = "SRX25994842"
    prompt = "\n".join([
        f"For the SRA accession {SRX_accession}, find the following information:",
        ] + get_metadata_items()
    )
    input = {
        "SRX": SRX_accession,
        "database": "sra",
        "entrez_id": "",
        "messages": [HumanMessage(content=prompt)],
    }
    graph = create_metadata_graph()
    for step in graph.stream(input, subgraphs=True, config={"max_concurrency" : 3, "recursion_limit": 40}):
        print(step)

    ## invoke with graph object directly provided
    invoke_metadata_graph = partial(invoke_metadata_graph, graph=graph)
    #print(invoke_metadata_graph(input))