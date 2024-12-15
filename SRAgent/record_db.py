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
from SRAgent.secret import get_secret

# Suppress the specific warning
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")

# functions
def db_connect() -> connection:
    """Connect to the sql database"""
    host = get_secret("GCP_SQL_DB_HOST", False)
    host = os.path.join(os.path.expanduser("~"), "cloudsql", host)
    db_params = {
        'host': host,
        'database': get_secret("GCP_SQL_DB_NAME"),
        'user':  get_secret("GCP_SQL_DB_USERNAME"),
        'password': get_secret("GCP_SQL_DB_PASSWORD"),
        'port': '5432',
        'connect_timeout': 10 
    }
    return psycopg2.connect(**db_params)

def db_list_tables(conn: connection) -> List[Tuple[str]]:
    """
    List all tables in the public schema of the database.
    Args:
        conn: Connection to the database.
    Returns:
        List of table names in the public schema.
    """
    tables = Table('tables', schema='information_schema')
    query = Query.from_(tables).select('table_name').where(tables.table_schema == 'public')
    with conn.cursor() as cur:
        cur.execute(str(query))
        tables = cur.fetchall()
        return tables

def db_glimpse_tables(conn: connection) -> None:
    """
    Print the first 5 rows of each table in the database.
    Args:
        conn: Connection to the database.
    """
    for table in db_list_tables(conn):
        table_name = table[0]
        print(f"#-- Table: {table[0]} --#")
        df = pd.read_sql(f"SELECT * FROM {table_name} LIMIT 5", conn)
        print(df)

def execute_query(stmt, conn: connection) -> Optional[List[Tuple]]:
    """
    Execute a query and return the results, if any.
    Args:
        stmt: Query to execute.
        conn: Connection to the database.
    Returns:
        Results of the query, if any.
    """
    try:
        with conn.cursor() as cur:
            cur.execute(str(stmt))
            conn.commit() 
            # Return the results of the query, if any
            try:
                return cur.fetchall()
            except psycopg2.ProgrammingError:
                return None
    except psycopg2.errors.DuplicateTable as e:
        print(f"Table already exists: {e}")

def db_get_srx_records(conn: connection, column: str="entrez_id", database: str="sra") -> List[int]:
    """
    Get the entrez_id values of all SRX records that have not been processed.
    Args:
        conn: Connection to the database.
        database: Name of the database to query.
    Returns:
        List of entrez_id values of SRX records that have not been processed.
    """
    srx_metadata = Table("srx_metadata")
    target_column = getattr(srx_metadata, column)
    stmt = Query \
        .from_(srx_metadata) \
        .select(target_column) \
        .distinct() \
        .where(srx_metadata.database == database)
        
    # Fetch the results and return a list of {target_column} values
    return [row[0] for row in execute_query(stmt, conn)]

def db_get_unprocessed_records(conn: connection, database: str="sra", max_records: int=3) -> pd.DataFrame:
    """
    Get all suitable unprocessed SRX records
    Args:
        conn: Connection to the database.
        database: Name of the database to query.
    Returns:
        Table of unprocessed SRX records.
    """
    srx_metadata = Table("srx_metadata")
    srx_srr = Table("srx_srr")

    stmt = Query \
        .from_(srx_metadata) \
        .inner_join(srx_srr) \
        .on(srx_metadata.srx_accession == srx_srr.srx_accession) \
        .where(Criterion.all([
            srx_metadata.screcounter_status.isnull(),
            srx_metadata.database == database,
            srx_metadata.is_illumina == "yes",
            srx_metadata.is_single_cell == "yes",
            srx_metadata.is_paired_end == "yes",
            srx_metadata.is_10x.isin(["yes", "unsure"])
        ])) \
        .select(
            srx_metadata.srx_accession.as_("sample"),
            srx_srr.srr_accession.as_("accession"),
            srx_metadata.entrez_id.as_("entrez_id"),
            srx_metadata.tech_10x.as_("lib_prep_method"),
            srx_metadata.organism.as_("organism")
        ) \
        .limit(max_records)
        
    # fetch as pandas dataframe
    return pd.read_sql(str(stmt), conn)
    
def get_unique_columns(table: str, conn: connection) -> List[str]:
    """
    Get all unique constraint columns for a table from the database schema.
    Prioritizes composite unique constraints over primary keys.
    
    Args:
        table: Name of the table
        conn: Database connection
        
    Returns:
        List of column names that form the most appropriate unique constraint
    """
    query = """
    SELECT c.contype, ARRAY_AGG(a.attname ORDER BY array_position(c.conkey, a.attnum)) as columns
    FROM pg_constraint c
    JOIN pg_class t ON c.conrelid = t.oid
    JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(c.conkey)
    WHERE t.relname = %s 
    AND c.contype IN ('p', 'u')  -- primary key or unique constraint
    GROUP BY c.conname, c.contype
    ORDER BY c.contype DESC;  -- 'u'nique before 'p'rimary key
    """
    
    with conn.cursor() as cur:
        cur.execute(query, (table,))
        constraints = cur.fetchall()
        
    if not constraints:
        raise ValueError(f"No unique constraints found in table {table}")
    
    # Prefer composite unique constraints over single-column primary keys
    for constraint_type, columns in constraints:
        if len(columns) > 1 or constraint_type == 'u':
            return columns
    
    # Fall back to primary key if no other suitable constraint found
    return constraints[0][1]

def db_add_update(data_list: List[Dict[str, Any]], table: str, conn: connection) -> None:
    """
    Add or update records in a database table. If a record with matching unique columns exists,
    it will be updated; otherwise, a new record will be inserted.

    Args:
        data_list: List of dictionaries to add to the table. If a key is missing from a dictionary,
                  the corresponding column in the table will be NULL.
        table: Name of the table to add the data to.
        conn: Connection to the database.
    """
    if not data_list:
        return

    # Get unique columns from schema
    unique_columns = get_unique_columns(table, conn)

    # Table and columns
    srx_metadata = Table(table)
    columns = list(data_list[0].keys())

    # Values for batch insert
    values = [tuple(d.get(col) for col in columns) for d in data_list]

    # Build the query with ON CONFLICT clause
    set_columns = [col for col in columns if col not in unique_columns]
    update_set = ", ".join(f"{col} = EXCLUDED.{col}" for col in set_columns)
    
    query = f"""
        INSERT INTO {srx_metadata} ({', '.join(columns)})
        VALUES %s
        ON CONFLICT ({', '.join(unique_columns)})
    """

    # Add DO UPDATE SET clause only if there are non-key columns to update
    if set_columns:
        update_set = ", ".join(f"{col} = EXCLUDED.{col}" for col in set_columns)
        query += f" DO UPDATE SET {update_set}"
    else:
        # If all columns are part of the unique constraint, do nothing on conflict
        query += " DO NOTHING"
    
    # Batch insert/update
    try:
        with conn.cursor() as cur:
            execute_values(cur, query, values)  
            conn.commit()
    except Exception as e:
        conn.rollback()
        raise e
    
def update_srx_record_status(conn: connection, srx_accession: str, status: str, log_msg: str) -> None:
    """
    Update the screcounter_status of a given SRX record.
    Args:
        conn: Connection to the database.
        srx_accession: SRX accession number of the record to update.
        status: New status of the record.
    """
    srx_metadata = Table("srx_metadata")
    stmt = Query \
        .update(srx_metadata) \
        .set(srx_metadata.screcounter_status, status) \
        .set(srx_metadata.screcounter_log, log_msg) \
        .where(srx_metadata.srx_accession == srx_accession)

    with conn.cursor() as cur:
        cur.execute(str(stmt))
        conn.commit()


def upsert_df(df: pd.DataFrame, table_name: str, conn: connection) -> None:
    """
    Upload a pandas DataFrame to PostgreSQL, performing an upsert operation.
    If records exist (based on unique constraints), update them; otherwise insert new records.
    
    Args:
        df: pandas DataFrame to upload
        table_name: name of the target table
        conn: psycopg2 connection object
    """    
    # Get DataFrame columns
    columns = list(df.columns)
    
    # Create ON CONFLICT clause based on unique constraints
    # For your ground_truth table, the unique constraints are (dataset_id, database, entrez_id)
    unique_columns = ["dataset_id", "database", "entrez_id"]

    # Drop duplicate records based on unique columns
    df.drop_duplicates(subset=unique_columns, keep='first', inplace=True)

        # Convert DataFrame to list of tuples
    values = [tuple(x) for x in df.to_numpy()]
    
    # Create the INSERT statement with ON CONFLICT clause
    insert_stmt = f"""
        INSERT INTO {table_name} ({', '.join(columns)})
        VALUES %s
        ON CONFLICT ({', '.join(unique_columns)})
        DO UPDATE SET
        {', '.join(f"{col} = EXCLUDED.{col}" 
                   for col in columns 
                   if col not in unique_columns + ['id'])}
    """
    
    try:
        with conn.cursor() as cur:
            # Use execute_values for better performance with multiple records
            execute_values(cur, insert_stmt, values)
            conn.commit()
    except Exception as e:
        conn.rollback()
        raise Exception(f"Error uploading data to {table_name}: {str(e)}")


# main
if __name__ == '__main__':
    # connect to database
    from dotenv import load_dotenv
    load_dotenv()

    # glimpse tables
    with db_connect() as conn:
        db_glimpse_tables(conn)
    #exit();

    # get processed entrez ids
    #with db_connect() as conn:
        #print(db_get_srx_records(conn, "srx_accession"))
        #print(db_get_unprocessed_records(conn))

    # update srx record status
    # with db_connect() as conn:
    #     update_srx_record_status(conn, "SRX26727599", "complete", "test log")

    # check unique columns
    # with db_connect() as conn:
    #     print(get_unique_columns("srx_metadata", conn))
    # exit();

    # update database
    # with db_connect() as conn:
        # db_add_update([{"database": "sra", "entrez_id": 36106630, "is_illumina": "yes"}], "srx_metadata", conn)
        # db_add_update([{"srx_accession": "SRX22716300", "srr_accession": "SRR27024456"}], "srx_srr", conn)
    # exit();

    data = [{
        "database": "sra",
        "entrez_id": 123456,
        "srx_accession": "test",
        "is_illumina": "unsure",
        "is_single_cell": "unsure",
        "is_paired_end": "unsure",
        "is_10x": "unsure",
        "tech_10x": "other",
        "organism": "other"
    }]
    #with db_connect() as conn:
    #    db_add_update(data, "srx_metadata", conn)

    df = pd.read_csv("data/ground_truth1.csv")
    df["dataset_id"] = "ground_truth1"
    # with db_connect() as conn:
    #     upsert_df(df, "ground_truth", conn)