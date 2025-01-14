# import
## batteries
import os
from typing import List, Dict, Any, Tuple, Optional
## 3rd party
import pandas as pd
from pypika import Query, Table, Field, Column, Criterion
import psycopg2
from psycopg2.extras import execute_values
from psycopg2.extensions import connection
## package
from SRAgent.db.utils import get_unique_columns


def db_update(df: pd.DataFrame, table_name: str, conn: connection) -> None:
    """
    Update existing records in a PostgreSQL table based on unique constraints.
    Args:
        df: pandas DataFrame with updated records
        table_name: name of the target table
        conn: psycopg2 connection object
    """
    if df.empty:
        return
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    # Filter columns to only those that exist in the table
    df = df[list(set(get_table_columns(table_name, conn)).intersection(df.columns))]

    # Sanitize integers, drop duplicates, etc.
    df = sanitize_int_columns(df.copy())
    unique_columns = get_unique_columns(table_name, conn)
    
    # Remove "id" 
    if "id" in df.columns:
        df = df.drop(columns=["id"])
    
    # Get non-unique columns
    columns = list(df.columns)
    non_unique_cols = [c for c in columns if c not in unique_columns]

    # Nothing to update if there are no non-unique columns or no rows
    if not non_unique_cols or df.empty:
        return

    # Convert DataFrame rows to tuples
    values = [tuple(x) for x in df.to_numpy()]
    
    # Build the WITH data(...) clause
    # E.g. WITH data(col_a, col_b, col_c) AS (VALUES %s)
    with_data_cols = ", ".join(columns)
    with_clause = f"WITH data({with_data_cols}) AS (VALUES %s)"

    # Build the UPDATE ... SET ... FROM data ... WHERE ...
    set_clause = ", ".join(f"{col} = data.{col}" for col in non_unique_cols)
    join_condition = " AND ".join(f"t.{uc} = data.{uc}" for uc in unique_columns)

    update_stmt = f"""
    {with_clause}
    UPDATE {table_name} AS t
    SET {set_clause}
    FROM data
    WHERE {join_condition}
    """

    try:
        with conn.cursor() as cur:
            execute_values(cur, update_stmt, values)
        conn.commit()
    except Exception as e:
        conn.rollback()
        raise Exception(f"Error updating data in {table_name}: {str(e)}")