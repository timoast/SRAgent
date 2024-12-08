# import
## batteries
import os
import sys
import time
import xml.etree.ElementTree as ET
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
## package
from SRAgent.tools.esearch import create_esearch_agent
from SRAgent.tools.esummary import create_esummary_agent
from SRAgent.tools.efetch import create_efetch_agent
from SRAgent.tools.elink import create_elink_agent
from SRAgent.tools.ncbi_fetch import create_ncbi_fetch_agent
from SRAgent.tools.entrez_db import which_entrez_databases
from SRAgent.tools.bigquery import create_bigquery_agent
from SRAgent.tools.utils import batch_ids

# functions
def elink(entrez_id, ntries=3):
    # Fetch detailed GEO record to get links to SRA
    error = None
    for _ in range(ntries):
        try:
            handle = Entrez.elink(
                id=entrez_id,
                dbfrom="gds",
                db="sra",
                retmode="xml"
            )
            batch_record = handle.read()
            handle.close()
        except Entrez.Parser.ValidationError:
            error = f"Failed to find links for IDs: {entrez_id}"
        except Exception as e:
            error = f"An error occurred: {e}"
        finally:
            try:
                handle.close()
            except:
                pass 

        if error is not None:
            time.sleep(1)
            continue

        # Parse XML
        try:
            root = ET.fromstring(batch_record)
        except ET.ParseError as e:
            error = f"XML Parsing Error: {e}"
            time.sleep(1)
            continue

        # Check for errors in the response
        error = root.find("ERROR")
        if error is not None:
            time.sleep(1)
            continue
        return root
    return error

@tool 
def geo2sra(
    entrez_ids: Annotated[List[str], "List of GEO (gds) database Entrez IDs"],
    )-> Annotated[List[str], "List linked SRA database Entrez IDs"]:
    """
    Convert GEO (gds) Entrez IDs to SRA Entrez IDs.
    """
    batch_size = 50
    id_map = []
    for entrez_id in entrez_ids:
        id_map.append(f"GEO Entrez ID: {entrez_id}; the associated SRA Entrez IDs:")

        # Fetch detailed GEO record to get links to SRA
        root = elink(entrez_id)
        if isinstance(root, str):
            id_map.append(root)
            continue
        
        # Extract links from the XML
        sra_ids = []
        for linksetdb in root.findall(".//LinkSetDb"):
            for link in linksetdb.findall("./Link"):
                sra_id = link.find("Id").text
                sra_ids.append(sra_id)

        # Add results to the output
        if not sra_ids:
            id_map.append(" - No SRA Entrez IDs found")
        else:
            id_map.append("\n".join([f" - {x}" for x in sra_ids]) + "\n")
        time.sleep(0.34) 

    # return SRA IDs
    return "\n".join(id_map)

def create_convert_agent(model_name="gpt-4o") -> Callable:
    # create model
    model_supervisor = ChatOpenAI(model=model_name, temperature=0.1)

    # set tools
    tools = [
        which_entrez_databases,
        create_ncbi_fetch_agent(),
        create_bigquery_agent(),
        create_esearch_agent(),
        create_esummary_agent(),
        create_efetch_agent(),
        create_elink_agent(),
    ]
  
    # state modifier
    state_mod = "\n".join([
        "You are a helpful senior bioinformatician assisting a researcher with a task involving Entrez databases.",
        "You have a team of agents who can perform specific tasks using Entrez tools.",
        "Provide guidance to the agents to help them complete the task successfully.",
        "\n",
        "Your main goal is to convert among Entrez IDs and accessions (e.g., SRX, SRR, GSE, GSM, PRJNA, SAMN).",
        "You can call the following agents to accomplish the task:",
        " - which_entrez_databases",
        " - bigquery",
        " - esearch", 
        " - esummary", 
        " - efetch",
        " - elink", 
        " - ncbi-fetch",
        "\n",
        "For SRA Entrez IDs, use the following strategy:",
        " 1. Use which_entrez_databases to determine which Entrez database contain the Entrez IDs (if not provided by the researcher).",
        " 2. Use esearch to find SRA accessions associated with the Entrez IDs.",
        " 3. Use bigquery to convert among SRA accessions in the SRA hierarchy: studies (SRP) → experiments (SRX) → runs (SRR)",
        " 4. If needed, use esummary and efetch to obtain detailed information on the accessions.",
        " 5. If needed, use elink to find related records in other databases.",
        "For GEO Entrez IDs, use the following strategy:",
         " 1. Use ncbi-fetch to find GEO accessions and SRA project IDs.",
         " 2. Use ncbi-fetch to find GEO accessions.",
        "\n",
        "Notes on strategy:\n",
        " - ENA records are equivalent to SRA records, so you can use the same strategy for both.",
        " - Be sure to provide context to the agents (e.g., \"Use esearch to find the Entrez ID for SRP123456\")."
        " - For Entrez agents (e.g., esearch or efetch), you should specify the database(s) to query (e.g., sra or gds).",
        " - If there are dozens of records, batch the IDs and call the agent multiple times to avoid rate limits and token count limits.",
        "\n",
        "Accession notes:",
        " - SRA accesssion prefixes: SRX, SRP, SRR",
        " - ENA accession prefixes: ERX, PRJNA, DRX, E-MTAB",
        " - GEO accession prefixes: GSE, GSM, GPL",
        " - BioProject accession prefixes: PRJNA, PRJEB, PRJDB",
        " - BioSample accession prefixes: SAMN, SAME",
        "\n",
        "Database notes:",
        " - Entrez databases: sra, gds, pubmed, biosample, bioproject",
        "\n",
        "Accession conversion workflows:",
        " - GSE → SRP → SRX → SRR",
        " - GSE → GSM → SRS → SRX → SRR",
        " - GSM → SRS → SRX → SRR",
        " - PRJNA → SRX → SRR",
        " - SAMN → SRX → SRR",
        " - ERP → SRP → SRX → SRR"
        "\n",
        "General notes:",
        " - Continue sending tasks to your agents until you successfully complete the task.",
        " - Be very concise; provide simple lists when possible; do not include unnecessary wording.",
        " - Simply list the converted accessions and Entrez IDs, unless told otherwise.",
        " - Write your output as plain text instead of markdown.",
        "\n",
    ])

    # create agent
    agent = create_react_agent(
        model=model_supervisor,
        tools=tools,
        state_modifier=state_mod
    )

    @tool
    def invoke_convert_agent(
        messages: List[BaseMessage]
    ) -> Annotated[str, "Response from the convert agent"]:
        """
        Invoke the convert agent to use Entrez tools to convert among Entrez IDs and accessions.
        Example 1: Convert SRX123456 to SRR accessions.
        Example 2: Convert Entrez ID 34747624 to SRX accessions.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": messages})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="convert_agent")]
        }
    return invoke_convert_agent


# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # test
    ## convert geo to sra
    #input = {'entrez_ids': ['200278601']}
    #input = {"entrez_ids" : ["35966237", "200254051"]}
    #input = {"entrez_ids" : ["200276533"]}
    #print(geo2sra.invoke(input)); 

    ## create agent
    agent = create_convert_agent()
    #message = "Convert SRX25716879 to SRR accessions"
    #message = "Convert Entrez ID 200276533 to SRX accessions. The Entrez ID is associated with the gds database."
    message = "Convert Entrez ID 200277303 to SRX accessions. The Entrez ID is associated with the gds database."
    input = {"messages" : [HumanMessage(content=message)]}
    print(agent.invoke(input))

