import re
import time
import requests
from bs4 import BeautifulSoup
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Callable
from langchain_core.tools import tool
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import HumanMessage, AIMessage


# functions
def _fetch_ncbi_record(
    term: Annotated[str, "The Entrez ID or SRA accession to fetch data for."],
    database: Annotated[str, "The NCBI database to fetch data from (sra or gds)."] = "sra",
) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/{database}/?term={term}"
    response = requests.get(url)
    if response.status_code != 200:
        return f"Error: Unable to fetch data for accession {term}. Status code: {response.status_code}"
    
    # final data
    gds_data = []

    # Locate the <p> tag with class 'details expand e-hidden'
    soup = BeautifulSoup(response.text, 'html.parser')
    section = soup.find('p', class_='details expand e-hidden')
    if section is None:
        section = soup.find('div', id='maincontent')
        # find all hrefs in the section
        if section is not None:
            hrefs = section.find_all('a')
            # find href starting with /gds/
            for href in hrefs:
                if href.get('href').startswith('/geo/query/'):
                    # extract the url and GEO accession
                    url = href.get("href")
                    GEO_accession = str(re.sub(r".+=", "", url))
                    url = f"https://www.ncbi.nlm.nih.gov{url}"
                    # fetch the page
                    response = requests.get(url)
                    if response.status_code != 200:
                        return f"Error: Unable to fetch data for accession {term}. Status code: {response.status_code}"
                    # parse the page
                    gds_section = _extract_geo_sections(response, GEO_accession)
                    gds_data.append(
                        f"\n\n#-- GEO accession for {term}: {GEO_accession} --#\n" + "\n\n".join(gds_section) 
                    )
                    break
    else:
        section = section.find_parent("div")
    if section is None:
        return f"Error: Unable to locate details for accession {term}."
    return re.sub(r"\n\n+", "\n\n", section.text.strip()) + "\n\n".join(gds_data)

@tool
def fetch_ncbi_record(
    terms: Annotated[List[str], "A list of Entrez IDs or SRA accessions to fetch data for."],
    database: Annotated[str, "The NCBI database to fetch data from (sra or gds)."] = "sra",
) -> str:
    """Fetches the NCBI SRA or GEO (gds) pages for a given accession number or Entrez ID."""
    data = []
    for term in terms:
        term = term.strip()
        data.append(f"#-- Query term: {term} --#")
        data.append(_fetch_ncbi_record(term, database) + "\n")
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
    return section.text

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

def _extract_geo_sections(response, GEO_accession):
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
    return sections

def _fetch_geo_record(
    GEO_accession: Annotated[str, "The GEO accession to fetch data for."]
) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO_accession}"
    response = requests.get(url)
    if response.status_code != 200:
        return f"Error: Unable to fetch data for GEO accession {GEO_accession}. Status code: {response.status_code}"
    return "\n\n".join(_extract_geo_sections(response, GEO_accession))

@tool
def fetch_geo_record(
    GEO_accessions: Annotated[List[str], "A list of GEO accessions to fetch data for."],
) -> str:
    """Fetches the NCBI GEO page for a given GEO accession number."""
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
        tools=[fetch_geo_record, fetch_ncbi_record, fetch_pubmed_record],
        state_modifier="\n".join([
            "You are an expert in bioinformatics and you are working on a project to find information about a specific dataset.",
            "You will use tools that directly request data from the NCBI website.",
            "You can query with both Entrez IDs and accessions.",
            "Use the following tools:",
            " - fetch_ncbi_record: fetch NCBI information on SRA and GEO records (SRA and GEO accessions or Entrez IDs).",
            "   - If no records are found for the provided database, try the other database.",
            " - fetch_pubmed_record: fetch information on PubMed records (SRA acccessions or Entrez IDs).",
            " - fetch_geo_record: fetch information on GEO accessions (do not provide Entrez IDs).",
            "General instructions:",
            " - Continue calling tools until you have all the information you need.",
            " - Provide a concise summary of your findings.",
            " - Use lists when possible.",
            " - Do not include helpful wording",
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

    # test agent
    invoke_ncbi_fetch_agent = create_ncbi_fetch_agent()
    #message = "Fetch information for Entrez ID 35447314"
    #message = "Fetch information for Entrez ID 200277303. The accession is associated with the gds database"
    #print(invoke_ncbi_fetch_agent.invoke(message))

    # test GEO fetch
    input = {"GEO_accessions" : ["GSE110878"]}
    #print(fetch_geo_record.invoke(input))

    # test fetch_sra_record
    input = {"terms" : ["200277303"], "database" : "gds"}
    # print(fetch_ncbi_record.invoke(input))

    # test fetch_pubmed_record
    input = {"terms" : ["34747624"]}
    #print(fetch_pubmed_record.invoke(input))

    