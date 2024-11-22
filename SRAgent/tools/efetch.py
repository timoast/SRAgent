# import
## batteries
import os
import re
import time
from pprint import pprint
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
## 3rd party
from Bio import Entrez
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage
## package
from SRAgent.tools.utils import batch_ids, truncate_values, xml2json


@tool 
def efetch(
    entrez_ids: Annotated[List[str], "List of Entrez IDs"],
    database: Annotated[str, "Database name (e.g., sra, gds, or pubmed)"],
) -> Annotated[str, "eFetch results in XML format"]:
    """
    Run an Entrez efetch query on Entrez IDs to obtain metadata for the records.
    Useful for obtaining metadata for specific records.
    """
    batch_size = 200  # Maximum number of IDs per request as per NCBI guidelines
    records = []
    regex = re.compile(r'"(PAIRED|SINGLE)": +null')

    for id_batch in batch_ids(entrez_ids, batch_size):
        time.sleep(0.34)  # Respect the rate limit of 3 requests per second
        id_str = ",".join(id_batch)
        try:
            # Fetch the records for the current batch of IDs
            handle = Entrez.efetch(db=database, id=id_str, retmode="xml")
            batch_record = handle.read()
            handle.close()
        except Entrez.Parser.ValidationError:
            batch_record = f"Failed to fetch record for IDs: {id_str}"
        except Exception as e:
            batch_record = f"An error occurred: {e}"
        finally:
            try:
                handle.close()
            except:
                pass  # Handle cases where handle might not be open

        # Decode the record if necessary
        if isinstance(batch_record, bytes):
            try:
                batch_record = batch_record.decode("utf-8")
            except Exception as e:
                batch_record = f"Decoding error: {e}"

        # Truncate long values in the record
        batch_record = truncate_values(batch_record, max_length=1000)
            
        # convert to XML to JSON
        batch_record = xml2json(batch_record)

        # fix values for PAIRED and SINGLE
        batch_record = re.sub(regex, "\\1: \"yes\"", batch_record)

        # Check for errors in the response
        if "Error occurred: cannot get document summary" in batch_record:
            batch_record = f"Failed to fetch record for IDs: {id_str}. Try a different database."

        records.append(batch_record)

    # Combine all records into a single string
    return "\n".join(records)

def create_efetch_agent(model_name: str="gpt-4o-mini") -> Callable:
    """
    Create an agent that uses Entrez efetch to help complete a task.
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.0)
    agent = create_react_agent(
        model=model,
        tools=[efetch],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "Based on the task provided by your supervisor, use Entrez efetch to help complete the task.",
            "You can use which_entrez_databases to determine which databases to use for efetch queries.",
            "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
        ])
    )

    @tool
    def invoke_efetch_agent(
        message: Annotated[str, "Message to the efetch agent"]
    ) -> Annotated[str, "Response from the efetch agent"]:
        """
        Invoke the efetch agent to run Entrez efetch queries.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="efetch agent")]
        }
    return invoke_efetch_agent

if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email = os.getenv("EMAIL")

    # Test efetch
    input = { "entrez_ids" : ["35966237"], "database" : "sra"}
    input = {"entrez_ids" : ["200254051"], "database" : "gds"}
    #print(efetch.invoke(input))

    input = {"entrez_ids" : ["27978912"], "database" : "sra"}
    print(efetch.invoke(input))