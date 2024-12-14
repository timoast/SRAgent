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
from SRAgent.tools.elink import elink
from SRAgent.tools.entrez_db import which_entrez_databases


def create_elink_agent(model_name: str="gpt-4o") -> Callable:
    """
    Create an agent that uses Entrez elink to help complete a task.
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.1)
    agent = create_react_agent(
        model=model,
        tools=[elink, which_entrez_databases],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "Based on the task provided by your supervisor, use Entrez elink to help complete the task.",
            "elink is useful for finding related entries between Entrez databases.",
            "elink requires Entrez IDs; if you are provided with SRA or GEO accessions, simply state that you need the Entrez IDs.",
            "If you are unsure of which database(s) to query (e.g., sra or gds), you can use which_entrez_databases to determine which databases contain the Entrez ID.",
            "Note that elink results are composed of Entrez IDs and not accessions (e.g., SRA accessions).",
            "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
        ])
    )

    @tool
    def invoke_elink_agent(
        message: Annotated[str, "Message to the elink agent"]
    ) -> Annotated[str, "Response from the elink agent"]:
        """
        Invoke the elink agent to run Entrez elink queries.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="elink_agent")]
        }
    return invoke_elink_agent


if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")
    Entrez.api_key = os.getenv('NCBI_API_KEY')

    # test elink agent
    agent = create_elink_agent()
    input = {"message": "Find SRA Entrez IDs for the GEO Entrez ID 200277303"}
    print(agent.invoke(input))
    
