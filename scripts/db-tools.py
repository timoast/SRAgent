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
from SRAgent.db.get import db_find_srx
from SRAgent.db.create import create_table, create_table_router

def parse_cli_args():

    # get current date in YYYY-MM-DD format
    today = datetime.today().strftime('%Y-%m-%d')

    # argparse
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):
        pass

    desc = 'Database tools'
    epi = """DESCRIPTION:
# list all tables
db-tools.py --list

# glimpse all tables
db-tools.py --glimpse

# dump all tables
db-tools.py --dump

# drop >=1 table 
db-tools.py --drop srx_metadata

# delete >=1 SRX accession
db-tools.py --delete-srx SRR1234567

# delete >=1 SRX accession from the scRecounter tables
db-tools.py --delete-srx-screcounter SRR1234567

# upsert from a csv
db-tools.py --upsert-target srx_metadata --upsert-csv db_bkup/2024-12-17/srx_metadata.csv
    """
    parser = argparse.ArgumentParser(
        description=desc, epilog=epi, formatter_class=CustomFormatter
    )
    parser.add_argument(
        '--dump', action="store_true", default=False, 
        help='Dump all tables to CSV files'
    )
    parser.add_argument(
        '--dump-dir', type=str, default=os.path.join("db_bkup", today), 
        help='Directory to dump CSV files'
    )
    parser.add_argument(
        '--list', action="store_true", default=False,
        help='List all tables in database'
    )
    parser.add_argument(
        '--glimpse', action="store_true", default=False,
        help='Glimpse all tables in database'
    )
    parser.add_argument(
        '--view', type=str, default=None,
        help='View table in database'
    )
    parser.add_argument(
        '--create', type=str, default=None, nargs='+',
        choices=list(create_table_router().keys()) + ["ALL"],
        help='Create >=1 table in the database'
    )
    parser.add_argument(
        '--find-srx', type=str, default=None, nargs='+', 
        help='Find SRX accessions in the database'
    )
    parser.add_argument(
        '--upsert-csv', type=str, default=None, 
        help='CSV file to upsert into database'
    )
    parser.add_argument(
        '--upsert-target', type=str, default=None, 
        help='Table to upsert into database'
    )
    parser.add_argument(
        '--drop', type=str, default=None, nargs='+', 
        help='>=1 table to delete from database.'
    )
    parser.add_argument(
        '--delete-srx', type=str, default=None, nargs='+', 
        help='>=1 SRX accession to delete from the entire database, except the Eval table.'
    )
    parser.add_argument(
        '--delete-srx-screcounter', type=str, default=None, nargs='+', 
        help='>=1 SRX accession to delete from the scRecounter tables in the database.'
    )
    parser.add_argument(
        '--tenant', type=str, default=os.getenv("DYNACONF"), 
        choices = ["test", "prod"],
        help='Database tenant to connect to. Defaults to DYNACONF env variable'
    )
    return parser.parse_args()

# functions
def dump_all_tables(dump_dir: str, tenant: str, conn: connection) -> List[str]:
    """Dump all tables to CSV files
    Args:
        dump_dir: Directory to dump CSV files
        tenant: Database tenant
        conn: Database connection
    Returns:
        List of paths to dumped CSV files
    """
    db_tables = db_list_tables(conn)
    dump_dir = os.path.join(dump_dir, tenant)
    os.makedirs(dump_dir, exist_ok=True)
    outfiles = []
    for table in db_tables:
        df = pd.read_sql(f"SELECT * FROM {table}", conn)
        outfile = os.path.join(dump_dir, f"{table}.csv")
        df.to_csv(outfile, index=False)
        print(f"Dumped {table} to {outfile}")
        outfiles.append(outfile)
    return outfiles

# functions
def main(args):
    # set pandas options
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.width", 100)

    # set tenant
    if args.tenant:
        os.environ["DYNACONF"] = args.tenant

    # list tables
    if args.list:
        with db_connect() as conn:
            print("Tables in the database:")
            print("\n".join(db_list_tables(conn)))

    # glimpse tables
    if args.glimpse:
        with db_connect() as conn:
            db_glimpse_tables(conn)

    # view table
    if args.view:
        with db_connect() as conn:
            tbl_names = db_list_tables(conn)
            if args.view not in tbl_names:
                print(f"Table {args.view} not found in database")
                sys.exit(1)
            df = pd.read_sql(f"SELECT * FROM {args.view}", conn)
            df.to_csv(sys.stdout, index=False)

    # find SRX accessions
    if args.find_srx:
        with db_connect() as conn:
            df = db_find_srx(args.find_srx, conn)
            df.to_csv(sys.stdout, index=False)

    # create table
    if args.create:
        with db_connect() as conn:
            for table_name in args.create:
                create_table(table_name, conn)
                print(f"Created table: {table_name}")

    # dump tables
    if args.dump:
        with db_connect() as conn:
            dump_all_tables(args.dump_dir, args.tenant, conn)

    # upsert tables
    if args.upsert_csv:
        if not args.upsert_target:
            print("Please provide --upsert-target")
            sys.exit(1)
        with db_connect() as conn:
            tbl_names = db_list_tables(conn)
            if args.upsert_target not in tbl_names:
                print(f"Table {args.upsert_target} not found in database")
                sys.exit(1)
            # determine separator from file extension
            df = pd.read_csv(args.upsert_csv, sep=get_sep(args.upsert_csv))
            db_upsert(df, args.upsert_target, conn)
            print(f"Upserted {args.upsert_csv} into {args.upsert_target}")

    # drop tables
    if args.drop:
        with db_connect() as conn:
            tbl_names = db_list_tables(conn)
            for table in args.drop:
                print(f"Attempting to drop: \"{table}\"...")
                if table not in tbl_names:
                    print(f"  Table {table} not found in database")
                    continue
                with conn.cursor() as cur:
                    cur.execute(f"DROP TABLE {table}")
                    conn.commit()
                print(f"  Dropped: {table}")
    
    # delete SRX accessions from all tables, except eval
    if args.delete_srx:
        with db_connect() as conn:
            for srx in args.delete_srx:
                # srx metadata & srx_srr
                for tbl_name in ["srx_metadata", "srx_srr"]:
                    with conn.cursor() as cur:
                        cur.execute(f"DELETE FROM {tbl_name} WHERE srx_accession = '{srx}'")
                        conn.commit()
                # screcounter_log & screcounter_star
                for tbl_name in ["screcounter_log", "screcounter_star_params", "screcounter_star_results"]:
                    with conn.cursor() as cur:
                        cur.execute(f"DELETE FROM {tbl_name} WHERE sample = '{srx}'")
                        conn.commit()
                print(f"Deleted: {srx}")
    
    # delete SRX accessions from scRecounter tables
    if args.delete_srx_screcounter:
        with db_connect() as conn:
            for srx in args.delete_srx_screcounter:
                # screcounter_log & screcounter_star
                for tbl_name in ["screcounter_log", "screcounter_star_params", "screcounter_star_results"]:
                    with conn.cursor() as cur:
                        cur.execute(f"DELETE FROM {tbl_name} WHERE sample = '{srx}'")
                        conn.commit()
                print(f"Deleted: {srx}")

def get_sep(infile: str) -> str:
    """Determine separator from file extension
    Args:
        infile: Input file
    Returns:
        Separator
    """
    sep = ","
    if infile.endswith(".csv"):
        sep = ","
    elif infile.endswith(".tsv"):
        sep = "\t"
    else:
        print("Input file must be CSV or TSV")
        sys.exit(1)
    return sep

# Example usage
if __name__ == "__main__":
    load_dotenv()
    args = parse_cli_args()
    main(args)