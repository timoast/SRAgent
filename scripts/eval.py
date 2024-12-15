#!/usr/bin/env python
# import
## batteries
import os
import sys
import argparse
from typing import List, Dict, Literal, Any
## 3rd party
from dotenv import load_dotenv
import pandas as pd
from tabulate import tabulate
from pypika import Query, Table, Criterion, functions as fn
## package
from SRAgent.record_db import db_connect, upsert_df


# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Evaluate a the accuracy of the data'
epi = """DESCRIPTION:

"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('--eval-dataset', type=str, help='Evaluation dataset ID')
parser.add_argument('--list-datasets', action='store_true', default=False,
                    help='List available evaluation datasets')
parser.add_argument('--add-dataset', type=str, default=None,
                    help='Provide a dataset csv to add to the database')
parser.add_argument('--outfile', type=str, default="incorrect.tsv",
                    help='Output file for incorrect predictions')


# functions
def add_suffix(columns: list, suffix: str="_x") -> list:
    """
    Add suffixes to duplicate column names.
    Args:
        columns: List of column names.
        suffix: Suffix to add to duplicate column names.
    Return:
        List of column names with suffixes added.
    """
    seen = set()
    result = []
    
    for col in columns:
        if col in seen:
            result.append(col + suffix)
        else:
            seen.add(col) 
            result.append(col)  
    return result

def eval(
    df: pd.dataframe, 
    exclude_cols: List[str]=["database", "entrez_id", "srx_accession"], 
    outfile: str="incorrect.tsv"
    ) -> None:
    """
    Evaluate the accuracy of the predictions.
    Args:
        df: DataFrame of the evaluation dataset.
        exclude_cols: Columns to exclude from evaluation.
        outfile: Output file for incorrect predictions.
    """
    # Get base columns (those without _pred suffix)
    #base_cols = [col for col in df.columns if not col.endswith('_pred') and col not in ['screcounter_status', 'screcounter_log']]
    base_cols = [col.replace("_pred", "") for col in df.columns if col.endswith('_pred')]

    # Create comparison for each column pair
    accuracy = {} 
    idx = set()
    for col in base_cols:
        if col in exclude_cols:
            continue
        pred_col = f"{col}_pred"
        if pred_col in df.columns:  # Check if prediction column exists
            # Compare values and show where they differ
            mismatches = df[df[col] != df[pred_col]]
            idx.update(mismatches.index)
            
            # Calculate mismatch percentage
            mismatch_pct = (len(mismatches) / len(df)) * 100
            accuracy[col] = 100.0 - mismatch_pct
            
            print(f"\n#-- {col} --#")
            print(f"# Total mismatches: {len(mismatches)} ({mismatch_pct:.2f}%)")
            
            if len(mismatches) > 0:
                # Display count of each 
                print("\n# Mismatches")
                df_mm = mismatches.groupby([col, pred_col]).size().reset_index(name="count")
                print(tabulate(df_mm.values, headers=df_mm.columns, tablefmt="github"))

    # convert to dataframe
    accuracy = pd.DataFrame(accuracy.items(), columns=["column", "accuracy_percent"])
    accuracy["accuracy_percent"] = accuracy["accuracy_percent"].round(2)
    accuracy["count"] = df.shape[0]

    # write to stdout
    print("\n#-- Accuracy Table --#")
    #accuracy.to_csv(sys.stdout, index=False, sep="\t")
    print(tabulate(accuracy.values, headers=df.columns, tablefmt="github"))

    # overall accuracy
    overall_accuracy = accuracy["accuracy_percent"].mean()
    print(f"\nOverall accuracy: {overall_accuracy:.2f}%")

    # print out the mismatch records
    print("\n#-- Mismatch Records --#")
    df_wrong = df.iloc[list(idx)]
    outdir = os.path.dirname(outfile)
    if outdir and outdir != ".":
        os.makedirs(outdir, exist_ok=True)
    df_wrong.to_csv(outfile, sep="\t", index=False)
    print(f"Saved mismatch records to: {outfile}")

def list_datasets() -> pd.DataFrame:
    """
    List available datasets in the database.
    Return:
        DataFrame  of dataset IDs and record counts.
    """
    with db_connect() as conn:
        tbl = Table("ground_truth")
        stmt = Query \
            .from_(tbl) \
            .select(tbl.dataset_id, fn.Count(tbl.dataset_id).as_("record_count")) \
            .groupby(tbl.dataset_id)
        datasets = pd.read_sql(str(stmt), conn)
        return datasets

def add_update_dataset(csv_file: str) -> None:
    """
    Add or update a dataset in the database.
    Args:
        csv_file: Path to the dataset CSV file.
    """
    # check if file exists
    if not os.path.exists(csv_file):
        print(f"File not found: {csv_file}")
        return None
    # load csv
    df = pd.read_csv(csv_file)
    dataset_id = os.path.splitext(os.path.split(csv_file)[1])[0]
    df["dataset_id"] = dataset_id
    # does dataset exist?
    existing_datasets = list_datasets()["dataset_id"].tolist()
    action = "Updated existing" if dataset_id in existing_datasets else "Added new"
    # add to database
    with db_connect() as conn:
        upsert_df(df, "ground_truth", conn)
    print(f"{action} dataset: {dataset_id}")

def load_eval_dataset(eval_dataset: str) -> pd.DataFrame:
    """
    Load the evaluation dataset from the database.
    Args:
        eval_dataset: Evaluation dataset ID.
    Return:
        DataFrame of the evaluation dataset
    """
    df = None
    with db_connect() as conn:
        tbl_eval = Table("ground_truth")
        tbl_pred = Table("srx_metadata")
        stmt = Query \
            .from_(tbl_eval) \
            .where(tbl_eval.dataset_id == eval_dataset) \
            .join(tbl_pred) \
            .on(
                (tbl_eval.database == tbl_pred.database) & 
                (tbl_eval.entrez_id == tbl_pred.entrez_id) &
                (tbl_eval.srx_accession == tbl_pred.srx_accession)
            ) \
            .select("*") 
        df = pd.read_sql(str(stmt), conn).drop("id", axis=1)
        df.columns = add_suffix(df.columns, "_pred")
    return df

def main(args):
    # add evaluation dataset
    if args.add_dataset:
        add_update_dataset(args.add_dataset)
        return None

    # list available datasets
    if args.list_datasets:
        print(list_datasets())
        return None

    if not args.eval_dataset:
        print("Please provide an evaluation dataset ID")

    # load evaluation dataset
    df = load_eval_dataset(args.eval_dataset)
    eval(df, outfile=args.outfile)


# Example usage
if __name__ == "__main__":
    args = parser.parse_args()
    load_dotenv()
    main(args)