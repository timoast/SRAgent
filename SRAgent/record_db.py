# import
## batteries
import os
import warnings
## 3rd party
import psycopg2
import pandas as pd
from pypika import Query, Table, Field, Column

# Suppress the specific warning
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")

# functions
def db_connect():
    host = os.path.join(os.path.expanduser("~"), "cloudsql", os.environ["GCP_SQL_DB_HOST"])
    db_params = {
        'host': host,
        'database': os.environ["GCP_SQL_DB_NAME"],
        'user': os.environ["GCP_SQL_DB_USERNAME"],
        'password': os.environ["GCP_SQL_DB_PASSWORD"],
        'port': '5432',
        'connect_timeout': 10 
    }
    return psycopg2.connect(**db_params)

def db_list_tables(conn):
    tables = Table('tables', schema='information_schema')
    query = Query.from_(tables).select('table_name').where(tables.table_schema == 'public')
    with conn.cursor() as cur:
        cur.execute(str(query))
        tables = cur.fetchall()
        return tables

def db_glimpse_tables(conn):
    for table in db_list_tables(conn):
        table_name = table[0]
        print(f"#-- Table: {table[0]} --#")
        df = pd.read_sql(f"SELECT * FROM {table_name} LIMIT 5", conn)
        print(df)

def execute_query(stmt, conn):
    try:
        with conn.cursor() as cur:
            cur.execute(str(stmt))
            conn.commit() 
            return cur.fetchall()
    except psycopg2.errors.DuplicateTable as e:
        print(f"Table already exists: {e}")

def db_get_processed_entrez_ids(conn, database='sra'):
    # SRX_SRR
    srx_metadata = Table("srx_metadata")
    stmt = Query \
        .from_(srx_metadata) \
        .select(srx_metadata.entrez_id) \
        .distinct() \
        .where((srx_metadata.processed != "complete") | (srx_metadata.processed.isnull())) \
        .where(srx_metadata.database == database)
        
    # Fetch the results and return a list of entrez_id values
    return [row[0] for row in execute_query(stmt, conn)]

# main
if __name__ == '__main__':
    # connect to database
    from dotenv import load_dotenv
    load_dotenv()

    # connect
    conn = db_connect()
    
    # glimpse tables
    #print(db_glimpse_tables(conn))

    # get processed entrez ids
    print(db_get_processed_entrez_ids(conn))

    # close connection
    conn.close()