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
def create_srx_metadata():
    # SRX_metadata
    stmt = Query \
        .create_table("srx_metadata") \
        .columns(
            #Column("id", "SERIAL", nullable=False),
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
        ) \
        .unique("database", "entrez_id")
    execute_query(stmt, conn)

def create_srr_srx():
    # SRX_SRR
    stmt = Query \
        .create_table("srx_srr") \
        .columns(
            Column("id", "SERIAL", nullable=False),
            Column("srx_accession", "VARCHAR(20)", nullable=False),
            Column("srr_accession", "VARCHAR(20)", nullable=False)
        ) \
        .unique("srx_accession", "srr_accession") \
        .primary_key("id")
    execute_query(stmt, conn)

def create_eval():
    # ground truth
    stmt = Query \
        .create_table("eval") \
        .columns(
            Column("id", "SERIAL", nullable=False),
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
        ) \
        .unique("dataset_id", "database", "entrez_id") \
        .primary_key("id")
    execute_query(stmt, conn)

# main
if __name__ == "__main__":
    from dotenv import load_dotenv
    from SRAgent.db.connect import db_connect
    load_dotenv()

    # connect to db
    with db_connect() as conn:
        # create tables
        create_srx_metadata()