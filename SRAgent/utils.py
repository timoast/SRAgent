import os
from typing import List
import pandas as pd
import gspread
from gspread_dataframe import set_with_dataframe


def filter_by_db(entrez_ids: List[str]) -> List[str]:
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