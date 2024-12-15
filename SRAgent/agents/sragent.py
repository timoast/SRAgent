# import
## batteries
import os
import asyncio
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
## package
from SRAgent.agents.entrez import create_entrez_agent
from SRAgent.agents.ncbi_fetch import create_ncbi_fetch_agent
from SRAgent.agents.bigquery import create_bigquery_agent
from SRAgent.agents.sequences import create_sequences_agent
from SRAgent.agents.utils import create_step_summary_chain

# functions
def create_sragent_agent(
    model_name="gpt-4o",
    return_tool: bool=True,
) -> Callable:
    # create model
    model = ChatOpenAI(model=model_name, temperature=0.1)

    # set tools
    tools = [
        create_entrez_agent(),
        create_ncbi_fetch_agent(),
        create_bigquery_agent(),
        create_sequences_agent()
    ]
  
    # state modifier
    state_mod = "\n".join([
        "# Introduction",
        " - You are a helpful senior bioinformatician assisting a researcher with a task involving NCBI Entrez databases.",
        " - You have a team of agents who can perform specific tasks using tools.",
        " - Your role is to coordinate these agents effectively to complete tasks, even if initial attempts fail.",
        "# Agents",
        "## Entrez Agent",
        "The Entrez agent can perform tasks using the Entrez tools:",
        "- Search for datasets in the NCBI SRA database",
        "- Retrieve metadata for datasets",
        "- Retrieve data for datasets",
        "- Link across databases",
        "Generally this agent works with Entrez IDs.",
        "## NCBI Fetch Agent",
        " - The NCBI Fetch agent can fetch data from the NCBI website.",
        " - The agent can use any set of Entrez IDs and/or SRA, ENA, or GEO accessions.",
        "## BigQuery Agent",
        " - The BigQuery agent can use the SRA BigQuery feature to find metadata on SRA and ENA datasets.",
        " - It can also be used to convert accessions among the SRA (and ENA) hierarchy: studies (SRP) → experiments (SRX) → runs (SRR).",
        "## Sequences Agent",
        " - The Sequences agent can fetch sequences and statistics on sequences (e.g., number of bases) from the NCBI database.",
        "# Strategy",
        " 1. Always try multiple approaches if the first attempt fails",
        " 2. Consider the following patterns:",
        "   - If searching for metadata, first try the BigQuery agent, and then try the Entrez and NCBI Fetch agents",
        "   - For accession conversions, try both BigQuery and Entrez agents",
        "   - For sequence data, combine Sequences agent with metadata from other agents",
        " 3. Keep track of what you've tried and what information you still need",
        " 4. If one agent returns partial information, use another agent to fill in the gaps",
        " 5. Cross-validate important information using multiple agents when possible",
        "# Execution rules",
        " 1. Before giving up, you MUST try at least two different agents or prompting approaches",
        " 2. Be sure to provide ALL important information each each agent, such as accessions, databases, or metadata fields",
        " 3. Do not simplify metadata fields; ask the agent for the EXACT metadata field",
        "   - For example, state \"Which 10X Genomics library preparation technology?\" to the agent instead of just \"10X Genomics data\"",
        " 4. If an agent returns no results, you MUST try a different prompt format or agent",
        " 5. Always verify if the obtained information fully answers the original question",
        "## Examples",
        "### Incorrect"
        "Fetch metadata for SRA accession SRX25994842 to find out if it is Illumina sequence data, 10X Genomics data, and the organism sequenced.",
        "### Correct",
        "Fetch the following metadata for SRA accession SRX25994842:",
        " - Illumina sequence data?",
        " - 10X Genomics data?",
        " - Which 10X Genomics library preparation technology?",
        " - The organism sequenced",
        "# Notes",
        " - Bulk RNA-seq is NOT the same as single-cell RNA-seq (scRNA-seq); be sure to distinguish between them.",
        "# Response",
        " 1. After each agent call, briefly analyze the response:",
        "   - What information was obtained?",
        "   - What information is still missing?",
        "   - Which agent should be tried next?",
        " 2. Be concise; use lists when possible",
        " 3. No markdown formatting",
        " 4. In final response, combine all gathered information"
    ])
    # create agent
    agent = create_react_agent(
        model=model,
        tools=tools,
        state_modifier=state_mod
    )

    # return agent instead of tool
    if not return_tool:
        return agent

    # create tool
    @tool
    async def invoke_sragent_agent(
        messages: Annotated[List[BaseMessage], "Messages to send to the SRAgent agent"],
    ) -> Annotated[dict, "Response from the SRAgent agent"]:
        """
        Invoke the SRAgent agent with a message.
        The SRAgent agent will perform a task using Entrez and other tools.
        """
        # Invoke the agent with the message
        result = await agent.ainvoke({"messages" : messages})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="sragent_agent")]
        }
    return invoke_sragent_agent


# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")


    async def main():
        # create entrez agent
        agent = create_sragent_agent()
    
        # invoke agent
        msg = "Convert GSE121737 to SRX accessions"
        # msg = "Is SRX20554853 paired-end Illumina data?"
        # msg = "Obtain all SRR accessions for SRX20554853"
        # msg = "List the collaborators for the SRX20554853 dataset"
        # msg = "\n".join([
        #     "For the SRA accession SRX25994842, find the following information:"
        #     " - Is the dataset Illumina sequence data?",
        #     " - Is the dataset single cell RNA-seq data?", 
        #     " - Is the dataset paired-end sequencing data?",
        #     " - Is the dataset 10X Genomics data?",
        #     " - Which 10X Genomics technology?",
        #     " - Which organism was sequenced?"
        # ])
        input = {"messages": [HumanMessage(content=msg)]}
        result = await agent.ainvoke(input)
        print(result)
        
    asyncio.run(main())