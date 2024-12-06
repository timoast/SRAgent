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
        """
        query = f"""
        WITH distinct_values AS (
            SELECT DISTINCT
                m.sra_study,
                m.experiment,
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
        """
        query = f"""
        SELECT 
            m.experiment,
            m.acc,
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
        "You are a helpful bioinformatician assisting a researcher with a task involving the Sequence Read Archive (SRA) database.",
        "You can use your tools to search the SRA with BigQuery.",
        "\n",
        "The tools provide metadata for SRA studies, experiments, and runs.",
        "For instance, you can obtain the organism, assay type, and number of bases for an SRA run accession.",
        "\n",
        "Note that you can convert between studies, experiments, and runs with the output of each tool.",
        "For instance, the get_study_metadata tool provides a list of experiments for each study.",
        "You can then use the get_experiment_metadata tool to get metadata for each experiment.",
        "\n",
        "Continue calling tools until you successfully complete the task.",
        "If you encounter an error, provide the error message to the researcher.",
        "If you cannot complete the task, inform the researcher of the issue.",
        "\n",
        "Be very concise; provide simple lists when possible; do not include unnecessary wording.",
        "Write your output as plain text instead of markdown.",
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
    # bigquery_agent = create_bigquery_agent()
    # print(bigquery_agent.invoke({"message" : "Get study metadata for SRP548813"}))

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