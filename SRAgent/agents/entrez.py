# import
## batteries
import os
from typing import Annotated, List, Dict, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
## package
from SRAgent.agents.esearch import create_esearch_agent
from SRAgent.agents.esummary import create_esummary_agent
from SRAgent.agents.efetch import create_efetch_agent
from SRAgent.agents.elink import create_elink_agent

# functions
def create_entrez_agent(model_name="gpt-4o") -> Callable:
    # create model
    model_supervisor = ChatOpenAI(model=model_name, temperature=0.1)

    # set tools
    tools = [
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
        "Generally, start with esearch to find Entrez records, then use efetch to get detailed information.",
        "Use esummary to obtain summary information on an Entrez record.",
        "Use elink to navigate between databases to find related records (e.g., GEO to SRA).",
        "\n",
        "Be sure to provide context to the agents (e.g., \"Use efetch to determine whether SRX4967527 is Illumina data.\")."
        "Generally, you will want to specify the database(s) to search (e.g., sra, gds, or pubmed).",
        "If there are dozens of records, batch the IDs and call the agent multiple times to avoid rate limits and token count limits.",
        "\n",
        "If the task involves accessions instead of Entrez IDs, you may need to convert them to Entrez IDs first.",
        "For example, convert SRX4967527 to the corresponding Entrez ID via eSearch of the SRA database.",
        "\n",
        "Continue sending tasks to your agents until you successfully complete the task.",
        "Be very concise; provide simple lists when possible; do not include unnecessary wording.",
        "Write your output as plain text instead of markdown.",
        "\n",
        "#Example workflows",
        "#@ Task: Convert GSE123456 to SRX, SRP, or SRR accessions",
        "  1. esearch agent: eSearch of the GSE accession to obtain Entrez IDs",
        "  2. esummary agent: eSummary of the Entrez IDs to get the SRX accessions",
        "#@ Task: Obtain the SRR accessions for SRX4967527",
        "  1. esearch agent: eSearch of the SRX accession to obtain the Entrez ID",
        "  2. efetch agent: eFetch of the Entrez ID to obtain the SRR accessions",
        "#@ Task: Is SRP309720 paired-end Illumina 10X Genomics data?",
        "  1. esearch agent: eSearch of SRP accession obtain the Entrez IDs",
        "  2. efetch agent: eFetch of the Entrez IDs to get the library preparation information",
        "#@ Task: Obtain the SRA study accessions for the Entrez ID 36098095",
        "  1. efetch agent: eFetch of the Entrez ID to obtain the SRA accessions",
    ])

    # create agent
    agent = create_react_agent(
        model=model_supervisor,
        tools=tools,
        state_modifier=state_mod
    )
    @tool
    def invoke_entrez_agent(
        message: Annotated[str, "Message to send to the Entrez agent"],
    ) -> Annotated[dict, "Response from the Entrez agent"]:
        """
        Invoke the Entrez agent with a message.
        The Entrez agent will perform a task using Entrez tools.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages" : [AIMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="entrez_agent")]
        }
    return invoke_entrez_agent

# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # create entrez agent
    agent = create_entrez_agent()
    
    # invoke agent
    #input = {"message": "Convert GSE121737 to SRX accessions"}
    #input = {"message": "Is SRX20554853 paired-end Illumina data?"}
    #input = {"message": "List the collaborators for the SRX20554853 dataset"}
    input = {"messages": [HumanMessage(content="How many bases per run in SRX20554853?")]}
    print(agent.invoke(input))