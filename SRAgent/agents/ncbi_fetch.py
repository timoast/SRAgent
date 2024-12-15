# import
## batteries
import asyncio
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
## package
from SRAgent.tools.ncbi_fetch import fetch_geo_record, fetch_ncbi_record, fetch_pubmed_record


def create_ncbi_fetch_agent(model_name: str="gpt-4o-mini") -> Callable:
    """
    Create an agent that queries the NCBI website 
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.1)
    agent = create_react_agent(
        model=model,
        tools=[fetch_geo_record, fetch_ncbi_record, fetch_pubmed_record],
        state_modifier="\n".join([
            "# Instructions"
            " - You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            " - You will use tools that directly request data from the NCBI website.",
            " - You can query with both Entrez IDs and accessions.",
            "# Tools",
            " - fetch_ncbi_record: fetch NCBI information on SRA and GEO records (SRA and GEO accessions or Entrez IDs).",
            "   - If no records are found for the provided database, try the other database.",
            " - fetch_pubmed_record: fetch information on PubMed records (SRA acccessions or Entrez IDs).",
            " - fetch_geo_record: fetch information on GEO accessions (do not provide Entrez IDs).",
            "# Strategy",
            " - Continue calling tools until you have all the information you need.",
            " - Be sure to provide all important information each each tool, such as accessions, databases, or metadata fields.",
            "# Notes",
            " - Bulk RNA-seq is NOT the same as single-cell RNA-seq (scRNA-seq); be sure to distinguish between them.",
            "   - If you do not find evidence of single cell, do not assume it is scRNA-seq.",
            "# Response",
            " - Provide a concise summary of your findings.",
            " - Use lists when possible.",
            " - Do not include helpful wording",
        ])
    )

    @tool
    async def invoke_ncbi_fetch_agent(
        message: Annotated[str, "Message to the ncbi-fetch agent"]
    ) -> Annotated[str, "Response from the ncbi-fetch agent"]:
        """
        Invoke the ncbi-fetch agent to query the NCBI website for information 
        on Entrez IDs, SRA accessions, and GEO accessions.
        """
        # Invoke the agent with the message
        result = await agent.ainvoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="ncbi-fetch_agent")]
        }
    return invoke_ncbi_fetch_agent


# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()

    # test agent
    async def main():
        invoke_ncbi_fetch_agent = create_ncbi_fetch_agent()
        #message = "Fetch information for Entrez ID 35447314"
        message = "Fetch information for Entrez ID 200277303. The accession is associated with the gds database"
        result = await invoke_ncbi_fetch_agent.ainvoke(message)
        print(result)
    asyncio.run(main())