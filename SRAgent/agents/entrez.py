# import
## batteries
import os
import sys
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any
## 3rd party
from Bio import Entrez
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
## package
from SRAgent.agents.workers import create_worker_agent
from SRAgent.agents.utils import create_step_summary_chain

# functions
def create_entrez_agent():
    # create model
    model_supervisor = ChatOpenAI(model="gpt-4o", temperature=0.1)

    # set tools
    tools = []
    for x in ["ncbi-fetch", "esearch", "esummary", "efetch", "elink"]:
        tools.append(create_worker_agent(x)[1])

    # state modifier
    state_mod = "\n".join([
        "You are a helpful senior bioinformatician assisting a researcher with a task involving Entrez databases.",
        "You have a team of workers who can perform specific tasks using Entrez tools.",
        "Provide guidance to the workers to help them complete the task successfully.",
        "\n",
        "Generally, start with eSearch to find Entrez records, then use eFetch to get detailed information.",
        "Use eSummary to obtain summary information on an Entrez record.",
        "Use eLink to navigate between databases to find related records (e.g., GEO to SRA).",
        "Use the ncbi-fetch worker to directly fetch data from the NCBI website (SRA, Pubmed, and GEO) and obtain more information.",
        "\n",
        "Be sure to provide context to the workers (e.g., \"Use efetch to determine whether SRX4967527 is Illumina data.\")."
        "Generally, you will want to specify the database(s) to search (e.g., sra, gds, or pubmed).",
        "If there are dozens of records, batch the IDs and call the worker multiple times to avoid rate limits and token count limits.",
        "\n",
        "If the task involves accessions instead of Entrez IDs, you may need to convert them to Entrez IDs first.",
        "For example, convert SRX4967527 to the corresponding Entrez ID via eSearch of the SRA database.",
        "\n",
        "Continue sending tasks to your workers until you successfully complete the task.",
        "For instance, if you cannot determine paired-end state from the eFetch worker, try using the ncbi-fetch worker.",
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
        "  1. esearch worker: eSearch of the GSE accession to obtain Entrez IDs",
        "  2. esummary worker: eSummary of the Entrez IDs to get the SRX accessions",
        "# Task: Obtain the SRR accessions for SRX4967527",
        "  1. esearch worker: eSearch of the SRX accession to obtain the Entrez ID",
        "  2. efetch worker: eFetch of the Entrez ID to obtain the SRR accessions",
        "# Task: Is SRP309720 paired-end Illumina 10X Genomics data?",
        "  1. esearch worker: eSearch of SRP accession obtain the Entrez IDs",
        "  2. efetch worker: eFetch of the Entrez IDs to get the library preparation information",
        "  3. ncbi-fetch worker: ncbi-fetch of the SRP accession to get more details",
        "# Task: Obtain the SRA study accessions for the Entrez ID 36098095",
        "  1. efetch worker: eFetch of the Entrez ID to obtain the SRA accessions",
        "  2. ncbi-fetch worker: ncbi-fetch of the SRP accessions to get more details",
    ])

    # create agent
    agent = create_react_agent(
        model=model_supervisor,
        tools=tools,
        state_modifier=state_mod
    )
    return agent

def invoke_entrez_agent(
    input: dict,
    agent: Any,
    step_summary_chain: Optional[Any]=None,
    config: dict = {"max_concurrency" : 8, "recursion_limit": 50}
) -> dict:
    """
    Invoke the Entrez agent to perform a task.
    """
    final_step = ""
    for i,step in enumerate(agent.stream(input, config=config)):
        final_step = step
        if step_summary_chain:
            msg = step_summary_chain.invoke({"step": step})
            print(f"Step {i+1}: {msg.content}", file=sys.stderr)
        else:
            print(f"Step {i+1}: {step}", file=sys.stderr)
    try:
        print(final_step["agent"]["messages"][-1].content)
    except:
        pass
    return final_step

# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # test step summary chain
    #msg = {'tools': {'messages': [ToolMessage(content="{'messages': [HumanMessage(content='- **Entrez ID: 200121737**\\n  - **SRX Accessions**: Not directly available, but related SRA ID is **SRP167700**\\n  - **GSE Accession**: GSE121737\\n  - **Samples**:\\n    - GSM3444963\\n    - GSM3444962\\n    - GSM3444964\\n\\n- **Entrez ID: 100024679**\\n  - **SRX Accessions**: Not directly available\\n  - **GSE Accession**: GSE132325; GSE151535; GSE206234; GSE240796; GSE192477; GSE206238; GSE166916; GSE121737; GSE184948\\n\\n- **Entrez ID: 303444964**\\n  - **SRX Accessions**: **SRX4967529**\\n  - **GSM Accession**: GSM3444964\\n\\n- **Entrez ID: 303444963**\\n  - **SRX Accessions**: **SRX4967528**\\n  - **GSM Accession**: GSM3444963\\n\\n- **Entrez ID: 303444962**\\n  - **SRX Accessions**: **SRX4967527**\\n  - **GSM Accession**: GSM3444962', additional_kwargs={}, response_metadata={}, name='esummary worker')]}", name='invoke_esummary_worker', id='dac6ce94-900b-4a87-a6f6-e48292ab1a83', tool_call_id='call_vPNZmoRFcXpgwmzxzEFdTMHj')]}}
    #step_summary_chain = create_step_summary_chain()
    #print(step_summary_chain.invoke({"step": msg, "max_tokens": 25}).content)

    # create entrez agent
    agent = create_entrez_agent()
    step_summary_chain = create_step_summary_chain()

    # invoke agent
    #input = {"messages": [("user", "Convert GSE121737 to SRX accessions")]}
    #input = {"messages": [("user", "Is SRX20554853 paired-end Illumina data?")]}
    input = {"messages": [("user", "List the collaborators for the SRX20554853 dataset")]}
    invoke_entrez_agent(input, agent, None) #step_summary_chain)

