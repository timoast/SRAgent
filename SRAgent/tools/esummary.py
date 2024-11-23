# import
## batteries
import os
import time
import json
from pprint import pprint
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
#from langchain_core.runnables.config import RunnableConfig
from langchain_core.tools import InjectedToolArg
## package
from SRAgent.tools.utils import batch_ids, truncate_values, xml2json

# functions
@tool
def esummary(
    entrez_ids: Annotated[List[str], "List of Entrez IDs"],
    database: Annotated[str, "Database name (e.g., sra, gds, or pubmed)"],
) -> Annotated[str, "eSummary results in XML format"]:
    """
    Run an Entrez esummary query on Entrez IDs to obtain summary information for the records.
    Useful for obtaining summary information for specific records.
    """
    batch_size = 200  # Maximum number of IDs per request as per NCBI guidelines
    max_string_length = 500  # Maximum length of a string in the record
    records = []
    
    for id_batch in batch_ids(entrez_ids, batch_size):
        time.sleep(0.34)  # Respect NCBI's rate limits (no more than 3 requests per second)
        id_str = ",".join(id_batch)
        
        try:
            # Fetch summary record for the current batch
            handle = Entrez.esummary(db=database, id=id_str, retmode="xml")
            batch_record = handle.read()
            handle.close()
        except Entrez.Parser.ValidationError:
            batch_record = f"Failed to fetch summary for IDs: {id_str}. Check if the IDs exist."
        except Exception as e:
            batch_record = f"An error occurred: {e}"
        finally:
            try:
                handle.close()
            except:
                pass  # Handle cases where the handle might not be open
        
        # Decode the record if necessary
        if isinstance(batch_record, bytes):
            try:
                batch_record = batch_record.decode("utf-8")
            except Exception as e:
                print(f"Decoding error: {e}")
                continue
            
        # Truncate long values in the record
        batch_record = truncate_values(batch_record, max_length=500)

        # convert to XML to JSON
        batch_record = xml2json(batch_record)

        # Check for errors in the response
        if "ERROR" in batch_record.upper() or "INVALID_ID" in batch_record.upper():
            batch_record = f"Failed to fetch summary for IDs: {id_str}. Try a different database or verify the IDs."

        # Append the batch record to the list of records
        records.append(batch_record)
    
    # Combine all batch records into a single string
    combined_records = "\n".join(records)
    return combined_records

def create_esummary_agent(model_name: str="gpt-4o-mini") -> Callable:
    """
    Create an agent that uses Entrez esummary to help complete a task.
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.0)
    agent = create_react_agent(
        model=model,
        tools=[esummary],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "Based on the task provided by your supervisor, use Entrez esearch to help complete the task.",
            "If the sra or gds database does not return findings, try the other database.",
            "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
        ])
    )

    @tool
    def invoke_esummary_agent(
        message: Annotated[str, "Message to the esummary agent"]
    ) -> Annotated[str, "Response from the esummary agent"]:
        """
        Invoke the esearch agent to perform a task.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="esummary_agent")]
        }
    return invoke_esummary_agent


if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # test esummary
    #input = {"entrez_ids" : ["35966237"], "database" : "sra"}
    #input = {"entrez_ids" : ["200121737"], "database" : "sra"}
    input = {"entrez_ids" : ["27978912"], "database" : "sra"}
    print(esummary.invoke(input))
