import os
from typing import List
import pandas as pd
import gspread


def filter_by_db(entrez_ids: List[str]) -> List[str]:
    """
    Filter out the Entrez IDs that are already in the gsheet database.
    Args:
        entrez_ids: list of Entrez IDs
    Return:
        list of Entrez IDs that are not in the database
    """
    # Authenticate and open the Google Sheet
    db_name = "SRAgent_database"
    gc = gspread.service_account(filename=os.getenv("GOOGLE_APPLICATION_CREDENTIALS"))
    sheet = gc.open(db_name)
    worksheet = sheet.worksheet("database")
    # Read existing data with the first row as header
    db = worksheet.get_all_values()
    db = pd.DataFrame(db[1:], columns=db[0])
    # Extract Entrez IDs
    db_entrez_ids = db["entrez_id"].tolist()
    # Filter 
    return [ID for ID in entrez_ids if ID not in db_entrez_ids]


def save_graph_image(graph, outfile: str="graph_image.png") -> None:
    """
    Save the langgraph graph as a PNG image.
    Args:
        graph: langgraph graph object
        outfile: path to save the image
    """
    with open(outfile, "wb") as file:
        png_data = graph.get_graph().draw_mermaid_png()
        file.write(png_data)