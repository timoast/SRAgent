# import
## batteries
import os
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
## package
from SRAgent.tools.esearch import create_esearch_agent
from SRAgent.tools.esummary import create_esummary_agent
from SRAgent.tools.efetch import create_efetch_agent
from SRAgent.tools.elink import create_elink_agent

# functions
def create_convert_agent(model_name="gpt-4o") -> Callable:
    # create model
    model_supervisor = ChatOpenAI(model=model_name, temperature=0.1)

    # set tools
    tools = [
        create_esearch_agent(),
        create_esummary_agent(),
        create_efetch_agent(),
        create_elink_agent()
    ]
  
    # state modifier
    state_mod = "\n".join([
        "You are a helpful senior bioinformatician assisting a researcher with a task involving Entrez databases.",
        "You have a team of agents who can perform specific tasks using Entrez tools.",
        "Provide guidance to the agents to help them complete the task successfully.",
        "\n",
        "Your main goal is to convert among Entrez IDs and accessions (e.g., SRX, SRR, GSE, GSM, PRJNA, SAMN).",
        "You can call the esearch, esummary, efetch, and elink agents to accomplish this task.",
        "\n",
        "Be sure to provide context to the agents (e.g., \"Use esearch to find the Entrez ID for SRP123456\")."
        "Generally, you will want to specify the database(s) to search (e.g., sra, gds, or pubmed).",
        "If there are dozens of records, batch the IDs and call the agent multiple times to avoid rate limits and token count limits.",
        "\n",
        "Continue sending tasks to your agents until you successfully complete the task.",
        "Be very concise; provide simple lists when possible; do not include unnecessary wording such as \"If you need further assistance\".",
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
        "# Task: Convert GSE123456 to SRX, SRP, or SRR accessions",
        "  1. esearch agent: eSearch of the GSE accession to obtain Entrez IDs",
        "  2. esummary agent: eSummary of the Entrez IDs to get the SRX accessions",
        "# Task: Obtain the SRR accessions for SRX4967527",
        "  1. esearch agent: eSearch of the SRX accession to obtain the Entrez ID",
        "  2. efetch agent: eFetch of the Entrez ID to obtain the SRR accessions",
    ])

    # create agent
    agent = create_react_agent(
        model=model_supervisor,
        tools=tools,
        state_modifier=state_mod
    )

    @tool
    def invoke_convert_agent(
        message: Annotated[str, "Message to the convert agent"]
    ) -> Annotated[str, "Response from the convert agent"]:
        """
        Invoke the convert agent to use Entrez tools to convert among Entrez IDs and accessions.
        Example 1: Convert SRX123456 to SRR accessions.
        Example 2: Convert Entrez ID 34747624 to SRX accessions.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
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
    # create agent
    agent = create_convert_agent()
    #message = "Convert SRX25716879 to SRR accessions"
    #print(agent(message))
    input = {"message" : "Convert SRX25716879 to SRR accessions"}
    print(agent.invoke(input))

