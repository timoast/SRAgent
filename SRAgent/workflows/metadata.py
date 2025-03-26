# import 
import os
import re
import asyncio
import operator
from enum import Enum
from typing import Annotated, List, Dict, Any, Sequence, TypedDict, Callable, Union, get_args, get_origin
import pandas as pd
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from pydantic import BaseModel, Field
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langgraph.graph import START, END, StateGraph, MessagesState
from langchain_core.runnables.config import RunnableConfig
## package
from SRAgent.agents.utils import set_model
from SRAgent.db.connect import db_connect 
from SRAgent.db.upsert import db_upsert
from SRAgent.agents.sragent import create_sragent_agent
from SRAgent.workflows.tissue_ontology import create_tissue_ontology_workflow

# classes
class YesNo(Enum):
    """Choices: yes, no, or unsure"""
    YES = "yes"
    NO = "no"
    UNSURE = "unsure"

class OrganismEnum(Enum):
    """Organism sequenced"""
    # mammals
    HUMAN = "Homo sapiens"
    MOUSE = "Mus musculus"
    RAT = "Rattus norvegicus"
    MACAQUE = "Macaca mulatta"
    MARMOSET = "Callithrix jacchus"
    HORSE = "Equus caballus"
    DOG = "Canis lupus"
    BOVINE = "Bos taurus"
    SHEEP = "Ovis aries"
    PIG = "Sus scrofa"
    RABBIT = "Oryctolagus cuniculus"
    NAKED_MOLE_RAT = "Heterocephalus glaber"
    CHIMPANZEE = "Pan troglodytes"
    GORILLA = "Gorilla gorilla"
    CAT = "Felis catus"     # NEW
    BONOBO = "Pan paniscus" # NEW
    GREEN_MONKEY = "Chlorocebus aethiops" # NEW
    GRAY_SHORT_TAILED_OPPOSUM = "Monodelphis domestica" # NEW
    VERVET_MONKEY = "Chlorocebus pygerythrus" # NEW
    GOAT = "Capra aegagrus" # NEW
    ALPACA = "Vicugna pacos" # NEW
    CHINCHILLA = "Chinchilla lanigera" # NEW
    DOMESTIC_GUINEA_PIG = "Cavia porcellus" # NEW
    GOLDEN_HAMSTER = "Mesocricetus auratus" # NEW
    EURASIAN_HEDGEHOG = "Erinaceus europaeus" # NEW
    AMERICAN_MINK = "Neovison vison" # NEW
    REDNECKED_WALLABY = "Macropus rufogriseus" # NEW
    SUNDA_PANGOLIN = "Manis javanica" # NEW
    PLATYPUS = "Ornithorhynchus anatinus" # NEW
    FERRET = "Mustela putorius" # NEW
    NORTHERN_TREE_SHREW = "Tupaia belangeri" # NEW
    # birds
    CHICKEN = "Gallus gallus"
    ZEBRAFINCH = "Taeniopygia guttata"   # NEW
    GOOSE = "Anser cygnoides"        # NEW
    DUCK = "Anas platyrhynchos"      # NEW
    # reptiles
    TURTLE = "Trachemys scripta"   # NEW
    # amphibians
    FROG = "Xenopus tropicalis"
    AXOLOTL = "Ambystoma mexicanum"   # NEW
    # fish
    ZEBRAFISH = "Danio rerio"
    SALMON = "Salmo salar"                  # NEW
    STICKLEBACK = "Gasterosteus aculeatus"   # NEW
    # invertebrates
    FRUIT_FLY = "Drosophila melanogaster"
    ROUNDWORM = "Caenorhabditis elegans"
    MOSQUITO = "Anopheles gambiae"
    BLOOD_FLUKE = "Schistosoma mansoni"
    # plants
    THALE_CRESS = "Arabidopsis thaliana"
    RICE = "Oryza sativa"
    TOMATO = "Solanum lycopersicum"
    CORN = "Zea mays" 
    # microorganisms
    METAGENOME = "metagenome"
    # other
    OTHER = "other"

class Tech10XEnum(Enum):
    """10X Genomics library preparation technology"""
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
    NA = "not_applicable"

class LibPrepEnum(Enum):
    """scRNA-seq library preparation technology"""
    TENX = "10x_Genomics"
    SMART_SEQ = "Smart-seq"
    SMART_SEQ2 = "Smart-seq2"
    SMART_SEQ3 = "Smart-seq3"
    CEL_SEQ = "CEL-seq"
    CEL_SEQ2 = "CEL-seq2"
    DROP_SEQ = "Drop-seq"
    IN_DROPS = "indrops"
    SCALE_BIO = "Scale Bio"
    PARSE = "Parse"
    PARSE_EVERCODE = "Parse_evercode"
    PARSE_SPLIT_SEQ = "Parse_split-seq"
    FLUENT = "Fluent"
    PLEXWELL = "plexWell"
    MARS_SEQ = "MARS-seq"
    BD_RHAPSODY = "BD_Rhapsody"
    OTHER = "other"
    NA = "not_applicable"

class CellPrepEnum(Enum):
    """Distinguishes between single nucleus and single cell RNA sequencing methods"""
    SINGLE_NUCLEUS = "single_nucleus"
    SINGLE_CELL = "single_cell" 
    UNSURE = "unsure"   
    NA = "not_applicable"

class PrimaryMetadataEnum(BaseModel):
    """Metadata to extract"""
    is_illumina: YesNo = Field(description="Is the dataset Illumina sequence data?")
    is_single_cell: YesNo = Field(description="Is the dataset single cell?")
    is_paired_end: YesNo = Field(description="Is the dataset paired-end?")
    lib_prep: LibPrepEnum = Field(description="The library preparation technology")
    tech_10x: Tech10XEnum = Field(description="The 10x Genomics technology")
    cell_prep: CellPrepEnum = Field(description="The cell preparation technology")

class SecondaryMetadataEnum(BaseModel):
    organism: OrganismEnum = Field(description="The organism sequenced")
    tissue: str = Field(description="The tissues where sequenced")
    disease: str = Field(description="The diseases of interest")
    perturbation: str = Field(description="The perturbations of interest")
    cell_line: str = Field(description="The cell lines of interest")

class TertiaryMetadataEnum(BaseModel):
    tissue_ontology_term_id: List[str] = Field(description="A list of tissue ontology terms")

class ChoicesEnum(Enum):
    """Choices for the router"""
    CONTINUE = "CONTINUE"
    STOP = "STOP"

class Choice(BaseModel):
    """Choice to continue or stop"""
    Choice: ChoicesEnum

class SRR(BaseModel):
    """SRR accessions"""
    SRR: List[str] = Field(description="A list of SRR accessions")

class GraphState(TypedDict):
    """Shared state of the agents in the graph"""
    messages: Annotated[Sequence[BaseMessage], operator.add]
    route: Annotated[str, "Route"]
    attempts: Annotated[int, "Number of attempts to extract metadata"]
    metadata_level: Annotated[str, "Metadata level"]
    # metadata
    ## IDs
    database: Annotated[str, "Database"]
    entrez_id: Annotated[str, "Entrez ID"]
    SRX: Annotated[str, "SRX accession to query"]
    SRR: Annotated[List[str], "SRR accessions for the SRX"]
    ## primary metadata
    is_illumina: Annotated[str, "Is the dataset Illumina sequence data?"]
    is_single_cell: Annotated[str, "Is the dataset single cell RNA-seq data?"]
    is_paired_end: Annotated[str, "Is the dataset paired-end sequencing data?"]
    lib_prep: Annotated[str, "Which scRNA-seq library preparation technology?"]
    tech_10x: Annotated[str, "If 10X Genomics, which particular 10X technology?"]
    organism: Annotated[str, "Which organism was sequenced?"]
    ## secondary metadata
    cell_prep: Annotated[str, "Single nucleus or single cell RNA sequencing?"]
    tissue: Annotated[str, "Which tissues were sequenced?"]
    disease: Annotated[str, "Any disease information?"]
    perturbation: Annotated[str, "Any treatment/perturbation information?"]
    cell_line: Annotated[str, "Any cell line information?"]
    ## tertiary metadata
    tissue_ontology_term_id: Annotated[List[str], "The ontology terms corresponding to the sequenced tissues"]
    

# functions
def get_metadata_items(metadata_level: str="primary") -> Dict[str, str]:
    """
    Get primary metadata items based on graph state annotations
    Return:
        A dictionary of metadata items
    """
    # which metadata items to include?
    if metadata_level == "primary":
        to_include = PrimaryMetadataEnum.model_fields.keys()
    elif metadata_level == "secondary":
        to_include = SecondaryMetadataEnum.model_fields.keys()
    elif metadata_level == "tertiary":
        to_include = TertiaryMetadataEnum.model_fields.keys()
    else:
        raise ValueError("The metadata_level must be 'primary', 'secondary', or 'tertiary'.")

    # get the metadata items
    metadata_items = {}
    for key, value in GraphState.__annotations__.items():
        if key not in to_include or not get_origin(value) is Annotated:
            continue
        metadata_items[key] = get_args(value)[1]
    return metadata_items

def create_sragent_agent_node():
    # create the agent
    agent = create_sragent_agent()

    # create the node function
    async def invoke_sragent_agent_node(state: GraphState) -> Dict[str, Any]:
        """Invoke the SRAgent to get the initial messages"""
    
        # create message prompt
        metadata_level = state.get("metadata_level", "primary")
        metadata_items = get_metadata_items(metadata_level).values()
        prompt = "\n".join([
            "# Instructions",
            f"For the SRA experiment accession {state['SRX']}, find the following dataset metadata:",
            "\n".join([f" - {x}" for x in metadata_items]),
            "# IMPORTANT NOTES",
            " - If the dataset is not single cell, then some of the other metadata fields may not be applicable",
            " - Try to confirm all metadata values with two data sources",
            " - Do NOT make assumptions about the metadata values; find explicit evidence",
        ])
        # call the agent
        response = await agent.ainvoke({"messages" : [HumanMessage(content=prompt)]})
        # return the last message in the response
        return {
            "messages" : [response["messages"][-1]],
            "metadata_level" : metadata_level
        }
    return invoke_sragent_agent_node

def max_str_len(x: str, max_len: int=300) -> str:
    """
    Find the maximum length string in a list. If the string is longer than max_len, truncate it with "..."
    Args:
        x: The string or list of strings
        max_len: The maximum length of the string
    Returns:
        The maximum length string
    """
    if isinstance(x, list):
        x = ",".join(x)
    if not isinstance(x, str):
        return x
    return x if len(x) <= max_len else x[:max_len-3] + "..."
    
def get_extracted_fields(response) -> Dict[str, str]:
    """
    Dynamically extract fields from the response model
    Args:
        response: The response object
    Return:
        A dictionary of the extracted fields: {field_name: field_value}
    """
    # get the extracted metadata fields
    fields = {}
    for field_name in response.model_fields.keys():
        # set the max string length
        if field_name == "organism":
            max_len = 100
        else:
            max_len = 300
        # get the field value
        field_value = getattr(response, field_name)
        # add to fields dict
        if hasattr(field_value, 'value'):
            fields[field_name] = max_str_len(field_value.value, max_len=max_len)
        else:
            fields[field_name] = max_str_len(field_value, max_len=max_len)
    return fields

def get_annot(key: str, state: dict) -> str:
    """If the key matches a graph state field, return the field annotation"""
    try:
        return get_args(GraphState.__annotations__[key])[1]
    except KeyError:
        return key

def create_get_metadata_node() -> Callable:
    """Create a node to extract metadata"""
    model = set_model(agent_name="metadata")

    async def invoke_get_metadata_node(state: GraphState, config: RunnableConfig) -> Dict[str, Any]:
        """Structured data extraction"""
        metadata_items = "\n".join([f" - {x}" for x in get_metadata_items(state["metadata_level"]).values()])
        # format prompt
        prompt = "\n".join([
            "# Instructions",
            " - Your job is to extract metadata from the provided text on a Sequence Read Archive (SRA) experiment.",
            " - The provided text is from 1 or more attempts to find the metadata, so you many need to combine information from multiple sources.",
            " - If there are multiple sources, use majority rules to determine the metadata values, but weigh ambiguous values less (e.g., \"unknown\", \"likely\", or \"assumed\").",
            " - If there is not enough information to determine the metadata, respond with \"unsure\" or \"other\", depending on the metadata field.",
            " - If the selected \"lib_prep\" field is NOT \"10X_Genomics\", the \"tech_10x\" field should be \"not_applicable\".",
            " - \"single cell\" typically refers to whole-cell sequencing; \"nucleus\" is usually stated if single nucleus sequencing.",
            " - Keep free text responses short; use less than 300 characters.",
            "# The specific metadata to extract",
            metadata_items
        ])
        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt),
            ("system", "\nHere are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
        ])
        prompt = prompt.format_messages(history=state["messages"]) 
        # call the model with a certain enum
        if state["metadata_level"] == "primary":
            selected_enum = PrimaryMetadataEnum
        elif state["metadata_level"] == "secondary":
            selected_enum = SecondaryMetadataEnum
        else:
            raise ValueError("The metadata_level must be 'primary' or 'secondary'.")
        response = await model.with_structured_output(selected_enum, strict=True).ainvoke(prompt)
        extracted_fields = get_extracted_fields(response)
        # check logic
        try:
            if extracted_fields["is_single_cell"] != "yes":
                extracted_fields["tech_10x"] = "not_applicable"
            if extracted_fields["is_single_cell"] == "yes" and extracted_fields["lib_prep"] == "not_applicable":
                extracted_fields["lib_prep"] = "other"
            if extracted_fields["lib_prep"] != "10x_Genomics":
                extracted_fields["tech_10x"] = "not_applicable"
            if extracted_fields["lib_prep"] == "10x_Genomics" and extracted_fields["tech_10x"] == "not_applicable":
                extracted_fields["tech_10x"] = "other"
        except KeyError:
            pass
        # create the natural language response message   
        message = "\n".join(
            ["# The extracted metadata:"] + 
            [f" - {get_annot(k, GraphState)}: {fmt(v)}" for k,v in extracted_fields.items()]
        )
        return {
            "messages" : [HumanMessage(content=message)],
            **extracted_fields
        }
    return invoke_get_metadata_node

def create_router_node() -> Callable:
    """Routing based on percieved completion of metadata extraction"""
    model = set_model(agent_name="metadata_router")

    async def invoke_router_node(state: GraphState):
        """
        Router for the graph
        """
        # no need to evaluate secondary metadata
        if state["metadata_level"] != "primary":
            return {
                "route": "STOP", 
                "attempts": state.get("attempts", 0) + 1, 
                "messages": [AIMessage(content="No evaluation needed for secondary metadata")]
            }

        # create prompt
        prompt = "\n".join([
            "# Instructions",
            " - You are a helpful bioinformatican who is evaluating the metadata extracted from the SRA experiment.",
            " - You will be provided with the extracted metadata and will determine if the metadata is complete.",
            " - Metadata values of \"unsure\" or \"other\" are considered incomplete.",
            " - \"not_applicable\" is considered complete if the metadata field is not applicable, given other metadata values.",
            " - If the metadata is incomplete, you will respond to let the system know if more information is needed.",
            "# Notes",
            " - The organism may be \"other\" if the organism is not a common model organism.",
            " - If the library preparation method is not 10X Genomics, then there is no need to provide a 10X technology.",
            "# Response",
            " - Based on your evaluation of the metadata, select \"STOP\" if the task is complete or \"CONTINUE\" if more information is needed.",
        ])
        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt),
            MessagesPlaceholder(variable_name="history")
        ])
        prompt = prompt.format_messages(history=[state["messages"][-1]]) 
        # call the model
        response = await model.with_structured_output(Choice, strict=True).ainvoke(prompt)
        # format the response
        if response.Choice.value == ChoicesEnum.CONTINUE.value:
            message = "\n".join([
                f"At least some of the metadata for {state['SRX']} is incomplete (uncertain).",
                "Please try to provide more information by using a different approach (e.g., different agent calls)."
            ])
        else:
            message = f"Enough of the metadata has been extracted for {state['SRX']}."
        return {
            "route": response.Choice.value, 
            "attempts": state.get("attempts", 0) + 1, 
            "messages": [HumanMessage(content=message)]
        }
    return invoke_router_node

def route_retry_metadata(state: GraphState) -> str:
    """Determine the route based on the current state of the conversation."""
    # move on to which node?
    next_node = "bump_metadata_level_node" if state["metadata_level"] == "primary" else "tissue_ontology_node"
    max_attempts = 2 if state["metadata_level"] == "primary" else 1
    # limit the number of attempts
    if state["attempts"] >= max_attempts:
        return next_node
    # return route choice
    return "sragent_agent_node" if state["route"] == "CONTINUE" else next_node 

def bump_metadata_level(state: GraphState) -> str:
    """Bump the metadata level"""
    return {
        "metadata_level": "secondary",
        "attempts" : 0
    }

def create_tissue_ontology_node() -> Callable:
    """Create a node to extract tissue ontology"""
    agent = create_tissue_ontology_workflow()
    
    # create the node function
    async def invoke_tissue_ontology_node(state: GraphState) -> Dict[str, Any]:
        """Invoke the tissue ontology workflow"""
        # create message prompt
        tissues = state.get("tissue")
        if tissues:
            tissues = fmt(tissues)
        else:
            return {"tissue_ontology_term_id" : []}
        organism = state.get("organism", "No organism provided")
        disease= state.get("disease", "No disease provided")
        perturbation = state.get("perturbation", "No perturbation provided")
        cell_line = state.get("cell_line", "No cell line provided")
        message = "\n".join([
            "# Primary information",
            f"The tissues: {tissues}",
            "# Secondary information",
            f"The organism: {organism}",
            f"The disease: {disease}",
            f"The perturbation: {perturbation}",
            f"The cell line: {cell_line}",
        ])
        # call the agent
        response = await agent.ainvoke({"messages" : [HumanMessage(content=message)]})
        # return the structured response (term IDs)
        return {"tissue_ontology_term_id" : response}

    return invoke_tissue_ontology_node

async def invoke_SRX2SRR_sragent_agent_node(state: GraphState) -> Dict[str, Any]:
    """Invoke the SRAgent to get the SRR accessions for the SRX accession"""
    # format the message
    suffix = "Generally, the bigquery agent can handle this task."
    if state["SRX"].startswith("SRX"):
        message = f"Find the SRR accessions for {state['SRX']}. Provide a list of SRR accessions. {suffix}"
    elif state["SRX"].startswith("ERX"):
        message = f"Find the ERR accessions for {state['SRX']}. Provide a list of ERR accessions. {suffix}"
    else:
        message = f"The wrong accession was provided: \"{state['SRX']}\". The accession must start with \"SRX\" or \"ERR\"."
    # call the agent
    agent = create_sragent_agent()
    response = await agent.ainvoke({"messages" : [HumanMessage(content=message)]})
    # extract all SRR/ERR accessions in the message
    regex = re.compile(r"(?:SRR|ERR)\d{4,}")
    SRR_acc = regex.findall(response["messages"][-1].content)
    return {"SRR" : list(set(SRR_acc))}

def add2db(state: GraphState, config: RunnableConfig):
    """Add results to the records database"""
    # upload SRX metadata to the database
    data = [{
        "database": state["database"],
        "entrez_id": int(state["entrez_id"]),
        "srx_accession": state["SRX"],
        "is_illumina": state["is_illumina"],
        "is_single_cell": state["is_single_cell"],
        "is_paired_end": state["is_paired_end"],
        "lib_prep": state["lib_prep"],
        "tech_10x": state["tech_10x"],
        "cell_prep": state["cell_prep"],
        "organism": state["organism"],
        "tissue": fmt(state["tissue"]),
        "tissue_ontology_term_id": fmt(state["tissue_ontology_term_id"]),
        "disease": fmt(state["disease"]),
        "perturbation": fmt(state["perturbation"]),
        "cell_line": fmt(state["cell_line"]),
        "notes": "Metadata obtained by SRAgent"
    }]
    data = pd.DataFrame(data)
    if config.get("configurable", {}).get("use_database"):
        with db_connect() as conn:
            db_upsert(data, "srx_metadata", conn)

    # Upload SRR accessions to the database
    data = []
    for srr_acc in state["SRR"]:
        data.append({
            "srx_accession" : state["SRX"],
            "srr_accession" : srr_acc
        })
    if config.get("configurable", {}).get("use_database") and config.get("configurable", {}).get("no_srr") != True:
        with db_connect() as conn:
            db_upsert(pd.DataFrame(data), "srx_srr", conn)

def fmt(x: Union[str, List[str]]) -> str:
    """If a list, join them with a comma into one string"""
    if type(x) != list:
        return x
    return ",".join([str(y) for y in x])

def final_state(state: GraphState):
    """Provide the final state"""
    # get the metadata fields
    metadata = []
    for k,v in get_metadata_items("primary").items():
        metadata.append(f" - {v}: {state[k]}")
    for k,v in get_metadata_items("secondary").items():
        metadata.append(f" - {v}: {state[k]}")
    for k,v in get_metadata_items("tertiary").items():
        try:
            result = ", ".join(state[k])
        except ValueError:
            result = state[k]
        metadata.append(f" - {v}: {result}")
    # create the message
    message = "\n".join([
        "# SRX accession: " + state["SRX"],
        " - SRR accessions: " + fmt(state["SRR"]),
    ] + metadata)
    return {"messages": [HumanMessage(content=message)]}

def create_metadata_graph(db_add: bool=True) -> StateGraph:
    """
    Create a graph to extract metadata from an SRX accession
    Args:
        db_add: Add the results to the records database
    Return:
        A langgraph state graph object
    """
    #-- graph --#
    workflow = StateGraph(GraphState)

    # nodes
    workflow.add_node("sragent_agent_node", create_sragent_agent_node())
    workflow.add_node("get_metadata_node", create_get_metadata_node())
    workflow.add_node("router_node", create_router_node())
    workflow.add_node("bump_metadata_level_node", bump_metadata_level)
    workflow.add_node("tissue_ontology_node", create_tissue_ontology_node())
    workflow.add_node("SRX2SRR_node", invoke_SRX2SRR_sragent_agent_node)
    if db_add:
       workflow.add_node("add2db_node", add2db)
    workflow.add_node("final_state_node", final_state)

    # edges
    workflow.add_edge(START, "sragent_agent_node")
    workflow.add_edge("sragent_agent_node", "get_metadata_node")
    workflow.add_edge("get_metadata_node", "router_node")
    workflow.add_conditional_edges("router_node", route_retry_metadata)
    workflow.add_edge("bump_metadata_level_node", "sragent_agent_node")
    workflow.add_edge("tissue_ontology_node", "SRX2SRR_node")
    if db_add:
       workflow.add_edge("SRX2SRR_node", "add2db_node")
       workflow.add_edge("add2db_node", "final_state_node")
    else:
       workflow.add_edge("SRX2SRR_node", "final_state_node")

    # compile the graph
    return workflow.compile()

async def invoke_metadata_graph(
    state: GraphState, 
    graph: StateGraph,
    to_return: List[str] = list(PrimaryMetadataEnum.model_fields.keys()) + list(SecondaryMetadataEnum.model_fields.keys()),
    config: RunnableConfig=None,
) -> Annotated[dict, "Response from the metadata graph"]:
    """
    Invoke the graph to obtain metadata for an SRX accession
    Args:
        state: The graph state
        graph: The graph object
        to_return: The metadata items to return
    Return:
        A dictionary of the metadata items
    """
    response = await graph.ainvoke(state, config=config)
    # filter the response to just certain graph state fields
    filtered_response = {key: [response[key]] for key in to_return}
    return filtered_response

# main
if __name__ == "__main__":
    from functools import partial
    from Bio import Entrez

    #-- setup --#
    from dotenv import load_dotenv
    load_dotenv(override=True)
    Entrez.email = os.getenv("EMAIL")

    #-- graph --#
    async def main():
        entrez_id = 18060880
        SRX_accession = "SRX13201194"
        input = {
            "database": "sra",
            "entrez_id": entrez_id,
            "SRX": SRX_accession,
            #"metadata_level": "primary",
        }
        graph = create_metadata_graph(db_add=False)
        config = {"max_concurrency" : 3, "recursion_limit": 50, "configurable": {"organisms": ["mouse", "rat"]}}
        async for step in graph.astream(input, subgraphs=False, config=config):
            print(step)
    asyncio.run(main())
    exit();

    # Save the graph image
    # from SRAgent.utils import save_graph_image
    # save_graph_image(graph)
    # exit();

    ## invoke with graph object directly provided
    #invoke_metadata_graph = partial(invoke_metadata_graph, graph=graph)
    #print(invoke_metadata_graph(input))

    #-- nodes --#
    msg = """# The extracted metadata:
   - Is the dataset Illumina sequence data?: yes
   - Is the dataset single cell RNA-seq data?: yes
   - Is the dataset paired-end sequencing data?: yes
   - Which scRNA-seq library preparation technology?: 10x_Genomics
   - If 10X Genomics, which particular 10X technologies?: 5_prime_gex
   - Single nucleus or single cell RNA sequencing?: single_cell
 """
    state = {
        "messages" : [HumanMessage(content=msg)],
        "SRX": "SRX22716300",
        "SRR": [],
        "metadata_level": "primary",
        "is_illumina": "unsure",
        "is_single_cell": "unsure",
        "is_paired_end": "unsure",
        "lib_prep": "other",
        "tech_10x": "other",
        "cell_prep": "unsure",
        "organism": "other",
        "tissue": "other",
        "disease": "other",
        "perturbation": "other",
        "cell_line": "other", 
    }
    #node = create_router_node()
    #print(node(state)); exit();