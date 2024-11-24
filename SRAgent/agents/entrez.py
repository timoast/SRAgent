# import
## batteries
import os
import sys
from functools import partial
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
from SRAgent.tools.seq import create_sequences_agent
from SRAgent.agents.utils import create_step_summary_chain

# functions
def create_entrez_agent(
    model_name="gpt-4o",
    return_tool: bool=True,
) -> Callable:
    # create model
    model_supervisor = ChatOpenAI(model=model_name, temperature=0.1)

    # set tools
    tools = [
        create_esearch_agent(),
        create_esummary_agent(),
        create_efetch_agent(),
        create_elink_agent(),
        create_ncbi_fetch_agent(),
        create_sequences_agent()
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
        "Use the ncbi-fetch agent to directly fetch data from the NCBI website (SRA, Pubmed, and GEO) and obtain more information.",
        "Use the sequences agent to obtain sequence data from the NCBI databases (e.g., check that reads are paired-end).",
        "\n",
        "Be sure to provide context to the agents (e.g., \"Use efetch to determine whether SRX4967527 is Illumina data.\")."
        "Generally, you will want to specify the database(s) to search (e.g., sra, gds, or pubmed).",
        "If there are dozens of records, batch the IDs and call the agent multiple times to avoid rate limits and token count limits.",
        "\n",
        "If the task involves accessions instead of Entrez IDs, you may need to convert them to Entrez IDs first.",
        "For example, convert SRX4967527 to the corresponding Entrez ID via eSearch of the SRA database.",
        "\n",
        "Continue sending tasks to your agents until you successfully complete the task.",
        "For instance, if you cannot determine paired-end state from the efetch agent, try using the ncbi-fetch worker.",
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
        "# Task: Is SRP309720 paired-end Illumina 10X Genomics data?",
        "  1. esearch agent: eSearch of SRP accession obtain the Entrez IDs",
        "  2. efetch agent: eFetch of the Entrez IDs to get the library preparation information",
        "  3. ncbi-fetch agent: ncbi-fetch of the SRP accession to get more details",
        "# Task: Obtain the SRA study accessions for the Entrez ID 36098095",
        "  1. efetch agent: eFetch of the Entrez ID to obtain the SRA accessions",
        "  2. ncbi-fetch agent: ncbi-fetch of the SRP accessions to get more details",
    ])

    # create agent
    agent = create_react_agent(
        model=model_supervisor,
        tools=tools,
        state_modifier=state_mod
    )

    # return agent
    if not return_tool:
        return agent

    @tool
    def invoke_entrez_agent(
        messages: List[BaseMessage],
        config: Optional[Dict[str, Any]] = None
    ) -> Annotated[dict, "Response from the Entrez agent"]:
        """
        Invoke the Entrez agent with a message.
        The Entrez agent will perform a task using Entrez tools.
        """
        input = {"messages": messages}        
        result = agent.invoke(input)
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="entrez_agent")]
        }
    return invoke_entrez_agent


async def create_entrez_agent_stream(input, config: dict={}, summarize_steps: bool=False) -> str:
    """
    Create an Entrez agent and stream the steps.
    Args:
        input: Input message to the agent.
        config: Configuration for the agent.
        summarize_steps: Whether to summarize the steps.
    Returns:
        The final step message.
    """
    # create entrez agent
    agent = create_entrez_agent(return_tool=False)

    # create step summary chain
    step_summary_chain = create_step_summary_chain() if summarize_steps else None
    
    # invoke agent
    step_cnt = 0
    final_step = ""
    async for step in agent.astream(input, stream_mode="values", config=config):
        step_cnt += 1
        final_step = step
        # summarize step
        if step_summary_chain:
            msg = step_summary_chain.invoke({"step": step})
            print(f"Step {step_cnt}: {msg.content}", file=sys.stderr)
        else:
            print(f"Step {step_cnt}: {step}", file=sys.stderr)
    try:
        final_step = final_step["agent"]["messages"][-1].content
    except KeyError:
        final_step = final_step["messages"][-1].content
    return final_step

# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # call streaming agent
    import asyncio
    input = {"messages": [HumanMessage(content="How many bases per run in SRX20554853?")]}
    config = {"max_concurrency" : 3, "recursion_limit": 40}
    results = asyncio.run(create_entrez_agent_stream(input, config, summarize_steps=True))
    print(results)
    
    