import time
import requests
from bs4 import BeautifulSoup
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage


# functions
def _fetch_sra_record(
    term: Annotated[str, "The Entrez ID or SRA accession to fetch data for."],
) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/sra/?term={term}"
    response = requests.get(url)
    if response.status_code != 200:
        return f"Error: Unable to fetch data for accession {term}. Status code: {response.status_code}"
    
    # Locate the <p> tag with class 'details expand e-hidden'
    soup = BeautifulSoup(response.text, 'html.parser')
    section = soup.find('p', class_='details expand e-hidden')
    if section is None:
        return f"Error: Unable to locate details for accession {term}."
    parent_div = section.find_parent("div")
    if parent_div is None:
        return f"Error: Unable to locate details for accession {term}."
    return parent_div.text

@tool
def fetch_sra_record(
    terms: Annotated[List[str], "A list of Entrez IDs or SRA accessions to fetch data for."],
) -> str:
    """Fetches the NCBI SRA page for a given accession number."""
    data = []
    for term in terms:
        term = term.strip()
        data.append(f"#-- Query term: {term} --#")
        data.append(_fetch_sra_record(term) + "\n")
        time.sleep(0.34)
    return "\n\n".join(data)
    
def _fetch_pubmed_record(
    term: Annotated[List[str], "The Entrez ID or PubMed ID to fetch data for."],
) -> str:
    """Fetches the NCBI PubMed page for a given PubMed ID."""
    url = f"https://pubmed.ncbi.nlm.nih.gov/{term}"
    response = requests.get(url)
    if response.status_code != 200:
        return f"Error: Unable to fetch data for PubMed ID {term}. Status code: {response.status_code}"
    
    # Locate the <div> tag with class 'abstract-content selected'
    soup = BeautifulSoup(response.text, 'html.parser')
    section = soup.find('div', class_='abstract-content selected')
    if section is None:
        return f"Error: Unable to locate details for PubMed ID {term}."
    return section.text.strip()

@tool
def fetch_pubmed_record(
    terms: Annotated[List[str], "A list of Entrez IDs or PubMed IDs to fetch data for."],
) -> str:
    """Fetches the NCBI PubMed page for a given PubMed ID."""
    data = []
    for term in terms:
        term = term.strip()
        data.append(f"#-- Query term: {term} --#")
        data.append(_fetch_pubmed_record(term) + "\n")
        time.sleep(0.34)
    return "\n\n".join(data)

def _fetch_geo_record(
    GEO_accession: Annotated[str, "The GEO accession to fetch data for."]
) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO_accession}"
    response = requests.get(url)
    if response.status_code != 200:
        return f"Error: Unable to fetch data for GEO accession {GEO_accession}. Status code: {response.status_code}"

    # extract the sections of interest
    section_names = [
        "Status", "Title", "Organism", "Experiment type", "Summary", "Overall design", 
        "Contributor(s)", "Citation(s)", "Platforms", "Samples", "BioProject", "SRA"
    ]
    sections = []
    soup = BeautifulSoup(response.text, 'html.parser')
    for section_name in section_names:
        for section in soup.find_all('tr'):
            text = section.text.strip()
            if text.startswith(section_name):
                sections.append("# " + text)
    if not sections:
        sections.append("No data found for GEO accession " + GEO_accession)
    return "\n\n".join(sections)

@tool
def fetch_geo_record(
    GEO_accessions: Annotated[List[str], "A list of GEO accessions to fetch data for."],
) -> str:
    """Fetches the NCBI PubMed page for a given PubMed ID."""
    data = []
    for acc in GEO_accessions:
        acc = acc.strip()
        data.append(f"#-- Query term: {acc} --#")
        data.append(_fetch_geo_record(acc) + "\n")
        time.sleep(0.34)
    return "\n\n".join(data)

def create_ncbi_fetch_agent(model_name: str="gpt-4o-mini") -> Callable:
    """
    Create an agent that queries the NCBI website 
    """
    model = ChatOpenAI(model_name=model_name, temperature=0.0)
    agent = create_react_agent(
        model=model,
        tools=[fetch_geo_record, fetch_sra_record, fetch_pubmed_record],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "You will use tools that directly request data from the NCBI website.",
            "You can query with both Entrez IDs and accessions.",
            "fetch_sra_record is useful for fetching information on SRA records (SRA accessions or Entrez IDs).",
            "fetch_pubmed_record is useful for fetching information on PubMed records (SRA acccessions or Entrez IDs).",
            "fetch_geo_record is useful for fetching information on GEO accessions (Entrez IDs are not).",
            "Provide a concise summary of your findings; use lists when possible; do not include helpful wording.",
        ])
    )

    @tool
    def invoke_ncbi_fetch_agent(
        message: Annotated[str, "Message to the ncbi-fetch agent"]
    ) -> Annotated[str, "Response from the ncbi-fetch agent"]:
        """
        Invoke the ncbi-fetch agent to query the NCBI website for information 
        on Entrez IDs, SRA accessions, and GEO accessions.
        """
        # Invoke the agent with the message
        result = agent.invoke({"messages": [HumanMessage(content=message)]})
        return {
            "messages": [AIMessage(content=result["messages"][-1].content, name="ncbi-fetch_agent")]
        }
    return invoke_ncbi_fetch_agent

# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()

    # test GEO fetch
    input = {"GEO_accessions" : ["GSE110878"]}
    print(fetch_geo_record.invoke(input))


    # test fetch_sra_record
    input = {"terms" : ["27978912"]}
    #print(fetch_sra_record.invoke(input))

    # test fetch_pubmed_record
    input = {"terms" : ["34747624"]}
    #print(fetch_pubmed_record.invoke(input))

    