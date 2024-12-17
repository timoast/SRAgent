# import
## batteries
import os
import warnings
from importlib import resources
from typing import List, Dict, Any, Tuple, Optional
## 3rd party
import pandas as pd
import psycopg2
from pypika import Query, Table, Field, Column, Criterion
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
from dynaconf import Dynaconf

# Suppress the specific warning
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")

# functions
def db_connect() -> connection:
    """Connect to the sql database"""
    s_path = None
    with resources.path("SRAgent", "settings.yml") as settings_path:
        s_path = str(settings_path)
    settings = Dynaconf(
        settings_files=["settings.yml", s_path], 
        environments=True, 
        env_switcher="DYNACONF"
    )
    # connect to db
    db_params = {
        'host': settings.db_host,
        'database': settings.db_name,
        'user': settings.db_user,
        'password': os.environ["GCP_SQL_DB_PASSWORD"],
        'port': settings.db_port,
        'connect_timeout': settings.db_timeout
    }
    return psycopg2.connect(**db_params)


# main
if __name__ == '__main__':
    from dotenv import load_dotenv
    load_dotenv()

    with db_connect() as conn:
        print(conn)
    