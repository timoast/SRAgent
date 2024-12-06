# import 
import os
import re
import operator
from functools import partial
from enum import Enum
from typing import Annotated, List, Dict, Any, Sequence, TypedDict, Callable, get_args, get_origin
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
    MACAQUE = "macaque"
    MARMOSET = "marmoset"
    HORSE = "horse"
    DOG = "dog"
    BOVINE = "bovine"
    CHICKEN = "chicken"
    SHEEP = "sheep"
    PIG = "pig"
    FRUIT_FLY = "fruit_fly"
    ROUNDWORM = "roundworm"
    ZEBRAFISH = "zebrafish"
    METAGENOME = "metagenome"
    OTHER = "other"

class Tech10XEnum(Enum):
    THREE_PRIME_GEX = "3_prime_gex"
    FIVE_PRIME_GEX = "5_prime_gex"
    ATAC = "atac"
    MULTIOME = "multiome"
    FLEX = "flex"
    VDJ = "vdj"
    FIXED_RNA = "fixed_rna"
    CELLPLEX = "cellplex"
    CNV = "cnv"
    FEATURE_BARCODING = "feature_barcoding"
    OTHER = "other"

class MetadataEnum(BaseModel):
    """
    Metadata to extract
    """
    is_illumina: YesNo
    is_single_cell: YesNo
    is_paired_end: YesNo
    is_10x: YesNo
    tech_10x: Tech10XEnum
    organism: OrganismEnum

class ChoicesEnum(Enum):
    CONTINUE = "Continue"
    STOP = "Stop"

class Choice(BaseModel):
    Choice: ChoicesEnum

class SRR(BaseModel):
    """
    SRR accessions
    """
    SRR: List[str]

class GraphState(TypedDict):
    """
    Shared state of the agents in the graph
    """
    messages: Annotated[Sequence[BaseMessage], operator.add]
    database: Annotated[str, "Database"]
    entrez_id: Annotated[str, "Entrez ID"]
    SRX: Annotated[str, "SRX accession to query"]
    SRR: Annotated[List[str], "SRR accessions for the SRX"]
    is_illumina: Annotated[str, "Is the dataset Illumina sequence data?"]
    is_single_cell: Annotated[str, "Is the dataset single cell RNA-seq data?"]
    is_paired_end: Annotated[str, "Is the dataset paired-end sequencing data?"]
    is_10x: Annotated[str, "Is the dataset 10X Genomics data?"]
    tech_10x: Annotated[str, "10X Genomics technology"]
    organism: Annotated[str, "Which organism was sequenced"]
    route: Annotated[str, "Route"]
    rounds: Annotated[int, operator.add]
    

# functions
def get_metadata_items() -> Dict[str, str]:
    """
    Set metadata items based on graph state annotationes
    """
    to_include = ["is_illumina", "is_single_cell", "is_paired_end", "is_10x", "tech_10x", "organism"]
    metadata_items = {}
    for key, value in GraphState.__annotations__.items():
        if key not in to_include or not get_origin(value) is Annotated:
            continue
        metadata_items[key] = get_args(value)[1]
    return metadata_items

def invoke_entrez_agent_node(state: GraphState) -> Dict[str, Any]:
    entrez_agent = create_entrez_agent()
    response = entrez_agent.invoke({"messages" : state["messages"]})
    return {"messages" : [response["messages"][-1]]}

def create_get_metadata_node() -> Callable:
    model = ChatOpenAI(model_name="gpt-4o", temperature=0)

    def invoke_get_metadata_node(state: GraphState):
        """
        Structured data extraction
        """
        # format prompt
        prompt = "\n".join([
            "Your job is to extract metadata from the provided text on a Sequence Read Archive (SRA) experiment.",
            "If there is not enough information to determine the metadata, please respond with 'unsure'.",
            "The specific metadata to extract:"] + list(get_metadata_items().values()))
        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt),
            ("system", "\nHere are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
        ])
        prompt = prompt.format_messages(history=state["messages"])

        # call the model
        response = model.with_structured_output(MetadataEnum, strict=True).invoke(prompt)
        return {
            "is_illumina" : response.is_illumina.value,
            "is_single_cell" : response.is_single_cell.value,
            "is_paired_end" : response.is_paired_end.value,
            "is_10x" : response.is_10x.value,
            "tech_10x" : response.tech_10x.value,
            "organism" : response.organism.value
        }
    
    return invoke_get_metadata_node

def create_router_node():
    model = ChatOpenAI(model_name="gpt-4o", temperature=0)

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
            " - Which 10X Genomics technology was used?",
            "    - " + state['tech_10x'],
            " - Which organism was sequenced?",
            "    - " + state['organism']
        ])

        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", prompt1),
            ("system", "\nHere are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
            ("system", prompt2),
            # Add the final instruction
            ("human", "Based on the information above, select STOP if the task is complete or CONTINUE if more information is needed."),
        ])
        formatted_prompt = prompt.format_messages(history=state["messages"])
        # call the model
        response = model.with_structured_output(Choice, strict=True).invoke(formatted_prompt)
        # format the response
        if response.Choice.value == ChoicesEnum.CONTINUE.value:
            message = "\n".join([
                f"At least some of the metadata for {state['SRX']} is still uncertain." 
                "Please try to provide more information by using a different approach (e.g., different tool calls)."
            ])
        else:
            message = f"Enough of the metadata has been extracted for {state['SRX']}."
        return {"route": response.Choice.value, "rounds": 1, "messages": [AIMessage(content=message)]}
    
    return invoke_router_node

def route_retry_metadata(state: GraphState) -> str:
    """
    Determine the route based on the current state of the conversation.
    """
    #continue_node = "add2db_node" if db_add else END
    if state["rounds"] >= 2:
        return "SRX2SRR_node"
    return "entrez_agent_node" if state["route"] == "Continue" else "SRX2SRR_node"

def invoke_SRX2SRR_entrez_agent_node(state: GraphState) -> Dict[str, Any]:
    # format the message
    if state["SRX"].startswith("SRX"):
        message = f"Find the SRR accessions for {state['SRX']}. Provide a list of SRR accessions."
    elif state["SRX"].startswith("ERX"):
        message = f"Find the ERR accessions for {state['SRX']}. Provide a list of ERR accessions."
    else:
        message = f"The wrong accession was provided: \"{state['SRX']}\". The accession must start with \"SRX\" or \"ERR\"."
    # call the agent
    entrez_agent = create_entrez_agent()
    response = entrez_agent.invoke({"messages" : [HumanMessage(content=message)]})
    # extract all SRR/ERR accessions in the message
    regex = re.compile(r"(?:SRR|ERR)\d{4,}")
    SRR_acc = regex.findall(response["messages"][-1].content)
    return {"SRR" : list(set(SRR_acc))}

def fmt(x):
    if type(x) != list:
        return x
    return ";".join([str(y) for y in x])

def add2db(state: GraphState):
    """
    Add results to the database
    """
    # create dataframe from state
    SRX = fmt(state["SRX"])
    results = pd.DataFrame([{
        "database": state["database"],
        "entrez_id": state["entrez_id"],
        "SRX": SRX,
        "SRR": fmt(state["SRR"]),
        "is_illumina": state["is_illumina"],
        "is_single_cell": state["is_single_cell"],
        "is_paired_end": state["is_paired_end"],
        "is_10x": state["is_10x"],
        "tech_10x": state["tech_10x"],
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
    return {"messages": [AIMessage(content=f"{SRX} added to database \"{db_name}\"")]}

def final_state(state: GraphState):
    """
    Final state
    """
    message = "\n".join([
        "# SRX accession: " + state["SRX"],
        " - SRR accessions: " + fmt(state["SRR"]),
        " - Is the study Illumina sequence data? " + state["is_illumina"],
        " - Is the study single cell RNA-seq data? " + state["is_single_cell"],
        " - Is the study paired-end sequencing data? " + state["is_paired_end"],
        " - Is the study 10X Genomics data? " + state["is_10x"],
        " - Which 10X Genomics technology was used? " + state["tech_10x"],
        " - Which organism was sequenced? " + state["organism"]
    ])
    return {"messages": [AIMessage(content=message)]}

def create_metadata_graph(db_add: bool=True):
    #-- graph --#
    workflow = StateGraph(GraphState)

    # nodes
    workflow.add_node("entrez_agent_node", invoke_entrez_agent_node)
    workflow.add_node("get_metadata_node", create_get_metadata_node())
    workflow.add_node("router_node", create_router_node())
    workflow.add_node("SRX2SRR_node", invoke_SRX2SRR_entrez_agent_node)
    if db_add:
        workflow.add_node("add2db_node", add2db)
    workflow.add_node("final_state_node", final_state)

    # edges
    workflow.add_edge(START, "entrez_agent_node")
    workflow.add_edge("entrez_agent_node", "get_metadata_node")
    workflow.add_edge("get_metadata_node", "router_node")
    workflow.add_conditional_edges("router_node", route_retry_metadata)
    if db_add:
        workflow.add_edge("SRX2SRR_node", "add2db_node")
        workflow.add_edge("add2db_node", "final_state_node")
    else:
        workflow.add_edge("SRX2SRR_node", "final_state_node")

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
    SRX_accession = "SRX25716878"
    #SRX_accession = "SRX20554856"
    prompt = "\n".join([
        f"For the SRA accession {SRX_accession}, find the following information:",
        ] + list(get_metadata_items().values())
    )
    input = {
        "SRX": SRX_accession,
        "database": "sra",
        "entrez_id": "",
        "messages": [HumanMessage(content=prompt)],
    }
    graph = create_metadata_graph(db_add=False)
    for step in graph.stream(input, subgraphs=True, config={"max_concurrency" : 3, "recursion_limit": 40}):
        print(step)

    ## invoke with graph object directly provided
    invoke_metadata_graph = partial(invoke_metadata_graph, graph=graph)
    #print(invoke_metadata_graph(input))

