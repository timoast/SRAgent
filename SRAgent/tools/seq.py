# import
## batteries
import os
import time
import tempfile
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
## package
from SRAgent.tools.convert import create_convert_agent
from SRAgent.tools.utils import run_cmd, truncate_values, xml2json

# functions
@tool
def fastq_dump(
    SRR_accessions: Annotated[List[str], "List of SRA run accessions (e.g., SRR1234567)"],
    tries: Annotated[int, "Number of attempts to run fastq-dump"]=3
) -> str:
    """
    Use fastq-dump to download the first few lines from the fastq files of the given SRR accession.
    The tool is useful for quickly checking the fastq files of an SRR accession.
    """
    # check if accession is valid
    incorrect_accessions = [x for x in SRR_accessions if not x.startswith("SRR")]
    if len(incorrect_accessions) > 0:
        acc_str = ", ".join(incorrect_accessions)
        return f"Invalid SRA accession numbers {acc_str}. Please provide >=1 valid SRR accession number."

    # create temp directory
    temp_dir = tempfile.TemporaryDirectory()

    # create command
    cmd = ["fastq-dump", "--outdir", temp_dir.name, "--split-files", "--maxSpotId", 2] + SRR_accessions

    # run command
    for i in range(tries):
        return_code, output, error = run_cmd(cmd)
        if return_code == 0:
            break
        time.sleep(5 * (i + 1))
    if return_code != 0:
        return f"Error running fastq-dump: {error.decode('utf-8')}"

    # read in the files
    files = os.listdir(temp_dir.name)
    if len(files) == 0:
        return "No FASTQ files found."
    fastq_data = ""
    for file in files:
        file_name = os.path.basename(file)
        with open(os.path.join(temp_dir.name, file), "r") as f:
            fastq_data += f"#-- File: {file_name} --#\n"
            fastq_data += f.read() + "\n"
            #fastq_data[file_name] = f.read()

    # delete the temp directory
    temp_dir.cleanup()
    return str(fastq_data)
  
@tool
def sra_stat(
    accessions: Annotated[List[str], "List of GEO and/or SRA accessions (e.g., SRP359840, SRR1234567, or GSE12345)"],
    tries: Annotated[int, "Number of attempts to run sra-stat"]=3
    ) -> str: 
    """
    Run the sra-stat CLI command (SRA Tools) on a GEO or SRA accession.
    Use this tool to get information about all sequence data associated with the accession.
    """
    # check if accession is valid
    incorrect_accessions = [x for x in accessions if not x.startswith(("SRP", "SRX", "SRR", "GSE", "GSM"))]
    if len(incorrect_accessions) > 0:
        acc_str = ", ".join(incorrect_accessions)
        return f"Invalid GEO/SRA accession numbers {acc_str}. Please provide >=1 valid GEO and/or SRA accession."

    # run sra-stat
    cmd = ['sra-stat', '--xml', '--quick'] + accessions

    # run command
    for i in range(tries):
        return_code, output, error = run_cmd(cmd)
        if return_code == 0:
            break
        time.sleep(5 * (i + 1))
    if return_code != 0:
        return f"Error running fastq-dump: {error.decode('utf-8')}"
        
    # Decode the record if necessary
    if isinstance(output, bytes):
        try:
            output = output.decode("utf-8")
        except Exception as e:
            return f"Decoding error: {e}"
            
    # Truncate long values in the record
    output = truncate_values(output, max_length=1000)

    # convert to XML to JSON
    output = xml2json(output)
    return str(output)

def create_sequences_agent(model_name: str="gpt-4o") -> Callable:
    """
    Create an agent that uses sra-stat and fastq-dump tools to provide information about SRA and GEO accessions.
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.0)
    agent = create_react_agent(
        model=model,
        tools=[sra_stat, fastq_dump, create_convert_agent()],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "Based on the task provided by your supervisor, use fastq-dump and sra-stat to help complete the task.",
            "You can investige the sequence data (fastq files) associated with GEO and/or SRA accessions.",
            "Use sra-stat to obtain sequence data information (e.g., number of spots and bases) associated with GEO and/or SRA accessions.",
            "fastq-dump is useful for quickly checking the fastq files of SRR accessions (e.g., is the data actually paired-end?).",
            "Note: fastq-dump only works with SRR accessions; use the convert tool to first convert Entrez IDs and accessions to SRR accessions.",
            "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
        ])
    )

    @tool
    def invoke_sequences_agent(
        message: Annotated[str, "Message to the sequences agent"]
    ) -> Annotated[str, "Response from the sequences agent"]:
        """
        Invoke the sequences agent to run the sra-stat and fastq-dump tools
        and provide information the sequence data (fastq files).
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="sequence_agent")]
        }
    return invoke_sequences_agent

if __name__ == "__main__":
    # Create the sequences agent
    sequences_agent = create_sequences_agent()
    input = {"message": "I need information about the SRA accession SRP359840."}
    print(sequences_agent.invoke(input))

    # fastq-dump
    #input = {"SRR_accessions" : ["SRR13112659", "SRR13112660"]}
    input = {"SRR_accessions" : ["ERR12363157"]}
    #print(fastq_dump.invoke(input))

    # sra-stat
    input = {"accessions" : ["SRP359840", "GSE12345"]}
    #print(sra_stat.invoke(input))

