#!/usr/bin/env python
# import
## batteries
import os
import sys
import argparse
from datetime import datetime
from typing import List, Dict, Literal, Any
## 3rd party
from dotenv import load_dotenv
import pandas as pd
from psycopg2.extensions import connection
from pypika import Query, Table, Criterion, functions as fn
## package
from SRAgent.record_db import db_connect, upsert_df, db_list_tables

# get current date in YYYY-MM-DD format
today = datetime.today().strftime('%Y-%m-%d')

# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Database tools'
epi = """DESCRIPTION:

"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('--dump', action="store_true", default=False, 
                    help='Dump all tables to CSV files')
parser.add_argument('--dump-dir', type=str, default=os.path.join("db", today),
                    help='Directory to dump CSV files')


def dump_all_tables(dump_dir: str, conn: connection) -> List[str]:
    """Dump all tables to CSV files
    Args:
        dump_dir: Directory to dump CSV files
        conn: Database connection
    Returns:
        List of paths to dumped CSV files
    """
    db_tables = db_list_tables(conn)
    os.makedirs(dump_dir, exist_ok=True)
    outfiles = []
    for table in db_tables:
        df = pd.read_sql(f"SELECT * FROM {table}", conn)
        outfile = os.path.join(dump_dir, f"{table}.csv")
        df.to_csv(outfile, index=False)
        print(f"Dumped {table} to {dump_dir}/{table}.csv")
        outfiles.append(outfile)
    return outfiles

# functions
def main(args):
    # set pandas options
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.width", 100)

    if args.dump:
        with db_connect() as conn:
            dump_all_tables(args.dump_dir, conn)


# Example usage
if __name__ == "__main__":
    args = parser.parse_args()
    load_dotenv()
    main(args)