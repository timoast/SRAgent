# import 
import os
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any
## 3rd-party
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage, ToolMessage
from langchain_openai import ChatOpenAI
from langchain_core.tools import tool
from Bio import Entrez
## package
from SRAgent.tools.esearch import esearch
from SRAgent.tools.efetch import efetch
from SRAgent.tools.esummary import esummary
from SRAgent.tools.elink import elink
from SRAgent.tools.entrez_db import which_entrez_databases
from SRAgent.tools.seq import sra_stat, fastq_dump


# create agent
def create_worker_agent(agent_name: str="esearch"):
    # create model
    model_complex = ChatOpenAI(model="gpt-4o", temperature=0)
    model_simple = ChatOpenAI(model="gpt-4o-mini", temperature=0)

    # create agent and invoke tool
    agent = None
    invoke_tool = None
    if agent_name == "esearch":
        # create agent
        agent = create_react_agent(
            model=model_simple,
            tools=[esearch],
            state_modifier="\n".join([
                "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
                "Based on the task provided by your supervisor, use Entrez esearch to help complete the task.",
                "If the sra or gds database does not return findings, try the other database.",
                "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
            ])
        )
        # create tool
        @tool
        def invoke_esearch_worker(
            message: Annotated[str, "Message to the worker"],
        ) -> Annotated[str, "Response from the worker"]:
            """
            Invoke the esearch worker to perform a task.
            """
            result = agent.invoke({"messages": [("user", message)]})
            # just return the response
            return {
                "messages": [HumanMessage(content=result["messages"][-1].content, name="esearch worker")]
            }
        invoke_tool = invoke_esearch_worker
    elif agent_name == "esummary":
        # create agent
        agent = create_react_agent(
            model=model_simple,
            tools=[esummary, which_entrez_databases],
            state_modifier="\n".join([
                "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
                "Based on the task provided by your supervisor, use Entrez esummary to help complete the task.",
                "You can use which_entrez_databases to determine which databases to use for esummary queries.",
                "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
            ])
        )
        # create tool
        @tool
        def invoke_esummary_worker(
            message: Annotated[str, "Message to the worker. Be sure to provide Entrez IDs."],
        ) -> Annotated[str, "Response from the worker"]:
            """
            Invoke the esummary worker to run eSummary on the provided Entrez ID(s).
            """
            result = agent.invoke({"messages": [("user", message)]})
            # just return the final response
            return {
                "messages": [HumanMessage(content=result["messages"][-1].content, name="esummary worker")]
            }
        invoke_tool = invoke_esummary_worker
    elif agent_name == "efetch":
        # create agent
        agent = create_react_agent(
            model=model_simple,
            tools=[efetch, which_entrez_databases],
            state_modifier="\n".join([
                "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
                "Based on the task provided by your supervisor, use Entrez efetch to help complete the task.",
                "You can use which_entrez_databases to determine which databases to use for efetch queries.",
                "Note that \"PAIRED: null\" does not mean that the data is single-end; it just means a lack of information.",
                "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
            ])
        )
        # create tool
        @tool
        def invoke_efetch_worker(
            message: Annotated[str, "Message to the worker. Be sure to provide Entrez IDs."],
        ) -> Annotated[str, "Response from the worker"]:
            """
            Invoke the efetch worker to run eFetch on the provided Entrez ID(s).
            """
            result = agent.invoke({"messages": [("user", message)]})
            # just return the final response
            return {
                "messages": [HumanMessage(content=result["messages"][-1].content, name="efetch worker")]
            }
        invoke_tool = invoke_efetch_worker
    elif agent_name == "elink":
        # create agent
        agent = create_react_agent(
            model=model_complex,
            tools=[elink, which_entrez_databases],
            state_modifier="\n".join([
                "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
                "Based on the task provided by your supervisor, use Entrez elink to help complete the task.",
                "elink is useful for finding related entries between Entrez databases.",
                "Generally, you will want to use the which_entrez_databases tool to determine which databases to use for elink queries.",
                "Note that elink results are composed of Entrez IDs and not accessions (e.g., SRA accessions).",
                "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
            ])
        )
        # create tool
        @tool
        def invoke_elink_worker(
            message: Annotated[str, "Message to the worker. Be sure to provide Entrez IDs."],
        ) -> Annotated[str, "Response from the worker"]:
            """
            Invoke the elink worker to run eLink on the provided Entrez ID(s).
            """
            result = agent.invoke({"messages": [("user", message)]})
            # just return the final response
            return {
                "messages": [HumanMessage(content=result["messages"][-1].content, name="elink worker")]
            }
        invoke_tool = invoke_elink_worker
    elif agent_name == "sequences":
        # create agent
        agent = create_react_agent(
            model=model_simple,
            tools=[sra_stat, fastq_dump],
            state_modifier="\n".join([
                "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
                "Based on the task provided by your supervisor, use sra-stat and fastq-dump to help complete the task.",
                "You can investige the sequence data (fastq files) associated with GEO and/or SRA accessions.",
                "sra-stat provides information about the sequence data associated with GEO and/or SRA accessions.",
                "fastq-dump is useful for quickly checking the fastq files of SRR accessions.",
                "If you are provided with Entrez IDs instead of GEO/SRA accessions, just state that you require GEO and/or SRA accessions.",
                "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
            ])
        )
        @tool
        def invoke_sequences_worker(
            message: Annotated[str, "Message to the worker. Be sure to provide Entrez IDs"],
        ) -> Annotated[str, "Response from the worker"]:
            """
            Invoke the sequences worker to run fastq-dump and/or sra-stat to investigate the sequence data.
            """
            result = agent.invoke({"messages": [("user", message)]})
            # just return the final response
            return {
                "messages": [HumanMessage(content=result["messages"][-1].content, name="sequences worker")]
            }
        invoke_tool = invoke_sequences_worker
    else:
        raise ValueError(f"Invalid agent name: {agent_name}")
    
    return agent, invoke_tool


# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # test
    agent, invoke_tool = create_worker_agent("esearch")
    input = {"message" : "Investigate GSE121737"}
    print(invoke_tool(input))

