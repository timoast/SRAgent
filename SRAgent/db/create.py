# import
## batteries
import os
import warnings
from typing import List, Dict, Any, Tuple, Optional
## 3rd party
import psycopg2
import pandas as pd
from pypika import Query, Table, Field, Column, Criterion
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
## package
from SRAgent.db.utils import execute_query

# functions
def create_updated_at_trigger(tbl_name: str, conn: connection) -> None:

    # Define the raw SQL for the trigger function and trigger
    trigger_function_sql = """
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;
    """

    trigger_sql = f"""
CREATE TRIGGER set_updated_at
BEFORE UPDATE ON {tbl_name}
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();
    """

    # Execute the SQL statements
    with conn.cursor() as cur:
        cur.execute(trigger_function_sql)
        cur.execute(trigger_sql)
        conn.commit()

def create_srx_metadata(conn: connection) -> None:
    tbl_name = "srx_metadata"
    stmt = Query \
        .create_table(tbl_name) \
        .columns(
            Column("database", "VARCHAR(20)", nullable=False),
            Column("entrez_id", "INT", nullable=False),
            Column("srx_accession", "VARCHAR(20)"),
            Column("is_illumina", "VARCHAR(10)"),
            Column("is_single_cell", "VARCHAR(10)"),
            Column("is_paired_end", "VARCHAR(10)"),
            Column("lib_prep", "VARCHAR(30)"),
            Column("tech_10x", "VARCHAR(30)"),
            Column("cell_prep", "VARCHAR(30)"),
            Column("organism", "VARCHAR(80)"),
            Column("tissue", "VARCHAR(80)"),
            Column("disease", "VARCHAR(100)"),
            Column("purturbation", "VARCHAR(100)"),
            Column("cell_line", "VARCHAR(100)"),
            Column("notes", "TEXT"),
            Column("created_at", "TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"),
            Column("updated_at", "TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"),
        ) \
        .unique("database", "entrez_id")
    
    execute_query(stmt, conn)
    create_updated_at_trigger(tbl_name, conn)

def create_srx_srr(conn: connection) -> None:
    tbl_name = "srx_srr"
    stmt = Query \
        .create_table(tbl_name) \
        .columns(
            Column("srx_accession", "VARCHAR(20)", nullable=False),
            Column("srr_accession", "VARCHAR(20)", nullable=False),
            Column("created_at", "TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"),
            Column("updated_at", "TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"),
        ) \
        .unique("srx_accession", "srr_accession")
    execute_query(stmt, conn)
    create_updated_at_trigger(tbl_name, conn)

def create_eval(conn: connection) -> None:
    tbl_name = "eval"
    stmt = Query \
        .create_table(tbl_name) \
        .columns(
            Column("dataset_id", "VARCHAR(30)", nullable=False),
            Column("database", "VARCHAR(20)", nullable=False),
            Column("entrez_id", "INT", nullable=False),
            Column("srx_accession", "VARCHAR(20)"),
            Column("is_illumina", "VARCHAR(10)"),
            Column("is_single_cell", "VARCHAR(10)"),
            Column("is_paired_end", "VARCHAR(10)"),
            Column("lib_prep", "VARCHAR(30)"),
            Column("tech_10x", "VARCHAR(30)"),
            Column("organism", "VARCHAR(80)"),
            Column("cell_prep", "VARCHAR(30)"),
            Column("created_at", "TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"),
            Column("updated_at", "TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP"),
        ) \
        .unique("dataset_id", "database", "entrez_id")
    execute_query(stmt, conn)
    create_updated_at_trigger(tbl_name, conn)

def create_table(table_name: str, conn: connection) -> None:
    if table_name == "srx_metadata":
        create_srx_metadata(conn)
    elif table_name == "srx_srr":
        create_srx_srr(conn)
    elif table_name == "eval":
        create_eval(conn)
    else:
        raise ValueError(f"Table {table_name} not recognized")

# main
if __name__ == "__main__":
    from dotenv import load_dotenv
    from SRAgent.db.connect import db_connect
    load_dotenv()

    # connect to db
    with db_connect() as conn:
        # create tables
        #create_srx_metadata()
        create_entrez_status()