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
from SRAgent.db.connect import db_connect
from SRAgent.db.upsert import db_upsert
from SRAgent.db.utils import db_list_tables, db_glimpse_tables, execute_query

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
parser.add_argument('--dump-dir', type=str, default=os.path.join("db_bkup", today),
                    help='Directory to dump CSV files')
parser.add_argument('--list', action="store_true", default=False,
                    help='List all tables in database')
parser.add_argument('--glimpse', action="store_true", default=False,
                    help='Glimpse all tables in database')
parser.add_argument('--upsert-csv', type=str, default=None,
                    help='CSV file to upsert into database')
parser.add_argument('--upsert-target', type=str, default=None,
                    help='Table to upsert into database')
parser.add_argument('--drop', type=str, default=None, nargs='+',
                    help='>=1 table to delete from database')



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

    # list tables
    if args.list:
        with db_connect() as conn:
            print("Tables in the database:")
            print("\n".join(db_list_tables(conn)))

    # glimpse tables
    if args.glimpse:
        with db_connect() as conn:
            db_glimpse_tables(conn)

    # dump tables
    if args.dump:
        with db_connect() as conn:
            dump_all_tables(args.dump_dir, conn)

    # delete tables
    if args.drop:
        with db_connect() as conn:
            tbl_names = db_list_tables(conn)
            for table in args.drop:
                if table not in tbl_names:
                    print(f"Table {table} not found in database")
                    continue
                #execute_query(f"DROP TABLE {table}", conn)
                with conn.cursor() as cur:
                    cur.execute(f"DROP TABLE {table}")
                    conn.commit()
                print(f"Dropped: {table}")

    # upsert tables
    if args.upsert_csv and args.upsert_target:
        with db_connect() as conn:
            tbl_names = db_list_tables(conn)
            if args.upsert_target not in tbl_names:
                print(f"Table {args.upsert_target} not found in database")
                sys.exit(1)
            df = pd.read_csv(args.upsert_csv)
            db_upsert(df, args.upsert_target, conn)
            print(f"Upserted {args.upsert_csv} into {args.upsert_target}")


# Example usage
if __name__ == "__main__":
    args = parser.parse_args()
    load_dotenv()
    main(args)