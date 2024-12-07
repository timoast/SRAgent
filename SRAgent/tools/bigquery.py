# import
## batteries
import os
import json
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
import decimal
from google.cloud import bigquery
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage


def to_json(results, indent: int=None):
    """
    Convert a dictionary to a JSON string.
    Args:
        results: a bigquery query result object
    Returns:
        str: JSON string
    """
    def datetime_handler(obj):
        if hasattr(obj, 'isoformat'):
            return obj.isoformat()
        elif isinstance(obj, decimal.Decimal):
            return str(obj)
        raise TypeError(f'Object of type {type(obj)} is not JSON serializable')

    return json.dumps(
        [dict(row) for row in results],
        default=datetime_handler,
        indent=indent
    )

def join_accs(accessions: List[str]) -> str:
    """
    Join a list of accessions into a string.
    Args:
        accessions: list of accessions
    Returns:
        str: comma separated string of accessions
    """
    return ', '.join([f"'{acc}'" for acc in accessions])

def create_get_study_metadata(client):
    @tool
    def get_study_metadata(
        study_accessions: Annotated[List[str], "A list of SRA study accession numbers (SRP)"]
        ) -> Annotated[str, "JSON string of SRA experiment metadata"]:
        """
        Get study-level metadata for a list of SRA study accessions.
        The metadata fields returned:
        - sra_study: SRA study accession (the query accession)
        - bioproject: BioProject accession (parent of study)
        - experiments: Comma-separated list of associated experiment accessions (SRX)
        """
        query = f"""
        WITH distinct_values AS (
            SELECT DISTINCT
                m.sra_study,
                m.bioproject,
                m.experiment
            FROM `nih-sra-datastore.sra.metadata` as m
            WHERE m.sra_study IN ({join_accs(study_accessions)})
        )
        SELECT 
            sra_study,
            bioproject,
            STRING_AGG(experiment, ',') as experiments
        FROM distinct_values
        GROUP BY sra_study, bioproject
        """
        return to_json(client.query(query))
    return get_study_metadata

def create_get_experiment_metadata(client):
    @tool
    def get_experiment_metadata(
        experiment_accessions: Annotated[List[str], "A list of SRA experiment accession numbers (SRX)"]
        ) -> Annotated[str, "JSON string of SRA experiment metadata"]:
        """
        Get experiment-level metadata for a list of SRA experiment accessions.
        The metadata fields returned:
        - experiment: SRA experiment accession (the query accession)
        - sra_study: SRA study accession (parent of experiment)
        - library_name: Library name (e.g., 1, 2, 3)
        - librarylayout: Library layout (e.g., single, paired)
        - libraryselection: Library selection (e.g., random, PCR)
        - librarysource: Library source (e.g., transcriptomic, genomic)
        - platform: Sequencing platform (e.g., Illumina, PacBio)
        - instrument: Sequencing instrument (e.g., HiSeq, NovaSeq)
        - acc: Comma-separated list of associated run accessions (SRR)
        """
        query = f"""
        WITH distinct_values AS (
            SELECT DISTINCT
                m.experiment,
                m.sra_study,
                m.library_name, 
                m.librarylayout,
                m.libraryselection, 
                m.librarysource,
                m.platform,
                m.instrument,
                m.acc,
            FROM `nih-sra-datastore.sra.metadata` as m
            WHERE m.experiment IN ({join_accs(experiment_accessions)})
        )
        SELECT
            experiment,
            library_name,
            librarylayout,
            libraryselection,
            librarysource,
            platform,
            instrument,
            STRING_AGG(acc, ',') as acc
        FROM distinct_values
        GROUP BY experiment, library_name, librarylayout, libraryselection, librarysource, platform, instrument
        """
        return to_json(client.query(query))
    return get_experiment_metadata

def create_get_run_metadata(client):
    @tool
    def get_run_metadata(
        run_accessions: Annotated[List[str], "A list of SRA run accession numbers (SRR)"]
        ) -> Annotated[str, "JSON string of SRA run metadata"]:
        """
        Get run-level metadata for a list of SRA run accessions.
        The metadata fields returned:
        - acc: SRA run accession (the query accession)
        - experiment: SRA experiment accession (parent of run)
        - biosample: BioSample accession (parent of run)
        - organism: Organism name
        - assay_type: Assay type (e.g., RNA-Seq, ChIP-Seq)
        - mbases: Total bases sequenced (in megabases)
        - avgspotlen: Average spot length (in base pairs)
        - insertsize: Insert size (in base pairs)
        """
        query = f"""
        SELECT 
            m.acc,
            m.experiment,
            m.biosample,
            m.organism,
            m.assay_type,
            m.mbases,
            m.avgspotlen,
            m.insertsize,            
        FROM `nih-sra-datastore.sra.metadata` as m
        WHERE m.acc IN ({join_accs(run_accessions)})
        """
        return to_json(client.query(query))
    return get_run_metadata

def create_bigquery_agent(model_name="gpt-4o") -> Callable:
    # create model
    model = ChatOpenAI(model=model_name, temperature=0.1)

    # init client
    client = bigquery.Client()

    # set tools
    tools = [
        create_get_study_metadata(client),
        create_get_experiment_metadata(client),
        create_get_run_metadata(client)
    ]
  
    # state modifier
    state_mod = "\n".join([
        # Role and Purpose
        "You are an expert bioinformatician specialized in querying the Sequence Read Archive (SRA) database.",
        "Your purpose is to retrieve and analyze metadata across SRA's hierarchical structure: studies (SRP) → experiments (SRX) → runs (SRR).",
        # Tool Capabilities
        "You have access to three specialized tools:",
        " 1. get_study_metadata: Retrieves study and associated experiment accessions",
        " 2. get_experiment_metadata: Retrieves experiment details and associated run accessions",
        " 3. get_run_metadata: Retrieves detailed run-level information",
        # Metadata structure
        "If the task is to retrieve metadata for a specific accession type (SRP, SRX, or SRR), chain the tools as needed to gather complete information.",
        # Conversion Strategy
        "IMPORTANT - Follow this strategy for accession conversion tasks:",
        " - To go from study → run: Use get_study_metadata → get_experiment_metadata → extract run accessions",
        " - To go from run → study: Use get_run_metadata → extract experiment → get_experiment_metadata → extract study",
        "Always convert accessions when needed to provide complete information.",
        # Response Guidelines
        "When responding:",
        " - If the query mentions one accession type but asks about another, automatically perform the necessary conversions",
        " - Chain multiple tool calls when needed to gather complete information",
        " - If you receive an error, explain it clearly and suggest alternatives",
        # Output Format
        "Keep responses concise and structured:",
        " - Present metadata as key-value pairs",
        " - Group related information",
        " - Include accession IDs in outputs",
        " - No markdown formatting",
    ])

    # create agent
    agent = create_react_agent(
        model=model,
        tools=tools,
        state_modifier=state_mod
    )
    @tool
    def invoke_bigquery_agent(
        message: Annotated[str, "Message to send to the BigQuery agent"],
    ) -> Annotated[dict, "Response from the BigQuery agent"]:
        """
        Invoke the BigQuery agent with a message.
        The BigQuery agent will search the SRA database with BigQuery.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages" : [AIMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="bigquery_agent")]
        }
    return invoke_bigquery_agent

if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    client = bigquery.Client()

    # test agent
    bigquery_agent = create_bigquery_agent()
    #print(bigquery_agent.invoke({"message" : "Get study metadata for SRP548813"}))
    print(bigquery_agent.invoke({"message" : "Get experiment metadata for SRP548813"}))

    # test tools
    ## get_study_metadata
    # get_study_metadata = create_get_study_metadata(client)
    # print(get_study_metadata.invoke({"study_accessions" : ["SRP548813"]}))

    ## get_experiment_metadata
    # get_experiment_metadata = create_get_experiment_metadata(client)
    # print(get_experiment_metadata.invoke({"experiment_accessions" : ["SRX26939191"]}))

    ## get_run_metadata
    # get_run_metadata = create_get_run_metadata(client)
    # print(get_run_metadata.invoke({"run_accessions" : ["SRR31573627", "SRR31573628"]}))