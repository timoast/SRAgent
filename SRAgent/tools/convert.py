# import
## batteries
import os
import sys
import time
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

# functions
@tool 
def geo2sra(
    entrez_ids: Annotated[List[str], "List of GEO (gds) database Entrez IDs"],
    )-> Annotated[List[str], "List linked SRA database Entrez IDs"]:
    """
    Convert GEO (gds) Entrez IDs to SRA Entrez IDs.
    """
    id_map = []
    for entrez_id in entrez_ids:
        # Fetch detailed GEO record to get links to SRA
        handle = Entrez.elink(dbfrom="gds", db="sra", id=entrez_id)
        links = Entrez.read(handle)
        handle.close()
        
        # header
        id_map.append(f"GEO Entrez ID: {entrez_id}; the associated SRA Entrez IDs:")

        # get SRA IDs
        IDs = []
        if links[0]['LinkSetDb']:
            try:
                IDs = [link['Id'] for link in links[0]['LinkSetDb'][0]['Link']]
            except (KeyError, IndexError) as e:
                print(f"Error: {e}", file=sys.stderr)
                pass
        if len(IDs) == 0:
            id_map.append(" - No SRA Entrez IDs found")
        else:
            id_map.append("\n".join([f" - {x}" for x in IDs]) + "\n")
        time.sleep(0.34) 

    # return SRA IDs
    return "\n".join(id_map)

def create_convert_agent(model_name="gpt-4o") -> Callable:
    # create model
    model_supervisor = ChatOpenAI(model=model_name, temperature=0.1)

    # set tools
    tools = [
        which_entrez_databases,
        geo2sra,
        create_ncbi_fetch_agent(),
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
        "You can call the ncbi-fetch, geo2sra, which_entrez_databases, esearch, esummary, efetch, elink agents to accomplish this task.",
        "Generally, start with the ncbi-fetch agent to check the NCBI website for infomation on the Entrez IDs and accessions.",
        "\n",
        "Be sure to provide context to the agents (e.g., \"Use esearch to find the Entrez ID for SRP123456\")."
        "Generally, you will want to specify the database(s) to query (e.g., sra, gds, or pubmed).",
        "If there are dozens of records, batch the IDs and call the agent multiple times to avoid rate limits and token count limits.",
        "\n",
        "Continue sending tasks to your agents until you successfully complete the task.",
        "Be very concise; provide simple lists when possible; do not include unnecessary wording.",
        "If you cannot complete the task, simply state that you were unable to do so.",
        "Write your output as plain text instead of markdown.",
        "\n",
        "#-- Accession notes --#",
        "SRA accesssion prefixes: SRX, SRP, SRR",
        "ENA accession prefixes: ERX, PRJNA, DRX, E-MTAB",
        "GEO accession prefixes: GSE, GSM, GPL",
        "BioProject accession prefixes: PRJNA, PRJEB, PRJDB",
        "BioSample accession prefixes: SAMN, SAME",
        "#-- Database notes --#",
        "Entrez databases: sra, gds, pubmed, biosample, bioproject",
        "#-- Accession conversion workflows --#",
        "GSE -> SRP -> SRX -> SRR",
        "GSE -> GSM -> SRS -> SRX -> SRR",
        "GSM -> SRS -> SRX -> SRR",
        "PRJNA -> SRX -> SRR",
        "SAMN -> SRX -> SRR",
        "ERP -> SRP -> SRX -> SRR",
        "#-- Example workflows --#",
        "# Task: Convert the Entrez ID 200249445 (gds database) to SRX accessions",
        "  1. geo2sra agent: Convert the GEO Entrez ID to SRA Entrez IDs",
        "  2. ncbi-fetch agent: check the NCBI website for basic information on the SRA Entrez IDs",
        "  3. esearch agent: eSearch of the SRA accessions to obtain the SRX accessions",
        "  4. efetch agent: eFetch of the SRA accessions to obtain the SRX accessions",
        "# Task: Convert the Entrez ID 200249001 (uknown database) to SRX accessions",
        "  1. which_entrez_databases agent: Check the available Entrez databases",
        "  2. geo2sra agent: If GEO Entrez ID, convert to SRA Entrez IDs",
        "  3. ncbi-fetch agent: check the NCBI website for basic information on the SRA Entrez IDs",
        "  4. esearch agent: eSearch of the SRA accessions to obtain the SRX accessions",
        "  5. efetch agent: eFetch of the SRA accessions to obtain the SRX accessions",
        "# Task: Convert GSE123456 to SRX accessions",
        "  1. esearch agent: eSearch of the GSE accession to obtain Entrez IDs",
        "  2. esummary agent: eSummary of the Entrez IDs to get the SRX accessions",
        "  3. ncbi-fetch agent: check the NCBI website for more information",
        "# Task: Obtain the SRR accessions for SRX4967527",
        "  1. esearch agent: eSearch of the SRX accession to obtain the Entrez ID",
        "  2. efetch agent: eFetch of the Entrez ID to obtain the SRR accessions",
        "  3. ncbi-fetch agent: check the NCBI website for more information",
        
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
    #print(geo2sra.invoke({"entrez_ids" : ["200254051", "200249445", "200276533"]})); 
    #exit();

    ## create agent
    agent = create_convert_agent()
    #message = "Convert SRX25716879 to SRR accessions"
    message = "Convert Entrez ID 200278601 to SRX accessions. The Entrez ID is associated with the gds database."
    input = {"messages" : [HumanMessage(content=message)]}
    print(agent.invoke(input))

