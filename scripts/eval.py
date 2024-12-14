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
from pypika import Query, Table, Criterion
## package
from SRAgent.record_db import db_connect


# argparse
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'Evaluate a the accuracy of the data'
epi = """DESCRIPTION:

"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
parser.add_argument('eval_dataset', type=str, help='Evaluation dataset ID')


# functions
def add_suffix(columns, suffix="_x"):
    """Add suffixes to duplicate column names."""
    seen = set()
    result = []
    
    for col in columns:
        if col in seen:
            result.append(col + suffix)
        else:
            seen.add(col) 
            result.append(col)  
    return result

def eval(df, exclude_cols=["database", "entrez_id", "srx_accession"]):
    # Get base columns (those without _pred suffix)
    #base_cols = [col for col in df.columns if not col.endswith('_pred') and col not in ['screcounter_status', 'screcounter_log']]
    base_cols = [col.replace("_pred", "") for col in df.columns if col.endswith('_pred')]

    # Create comparison for each column pair
    accuracy = {} 
    for col in base_cols:
        if col in exclude_cols:
            continue
        pred_col = f"{col}_pred"
        if pred_col in df.columns:  # Check if prediction column exists
            # Compare values and show where they differ
            mismatches = df[df[col] != df[pred_col]]
            
            # Calculate mismatch percentage
            mismatch_pct = (len(mismatches) / len(df)) * 100
            accuracy[col] = 100.0 - mismatch_pct
            
            print(f"\nComparison for {col}:")
            print(f"Total mismatches: {len(mismatches)} ({mismatch_pct:.2f}%)")
            
            if len(mismatches) > 0:
                # Display sample of mismatches
                print("\n--Mismatches--")
                mismatches[[col, pred_col]].head(10).transpose().to_csv(sys.stdout, sep="\t", header=False)

    # convert to dataframe
    accuracy = pd.DataFrame(accuracy.items(), columns=["column", "accuracy_percent"])
    accuracy["accuracy_percent"] = accuracy["accuracy_percent"].round(2)
    accuracy["count"] = df.shape[0]

    # write to stdout
    print("\n#-- Accuracy Table --#")
    accuracy.to_csv(sys.stdout, index=False, sep="\t")

    # overall accuracy
    overall_accuracy = accuracy["accuracy_percent"].mean()
    print(f"\nOverall accuracy: {overall_accuracy:.2f}%")

    # TODO: include which records are problematic

def main(args):
    # load evaluation dataset
    df = None
    with db_connect() as conn:
        tbl_eval = Table("ground_truth")
        #eval_columns = pd.read_sql("SELECT * FROM ground_truth LIMIT 0", conn).columns

        tbl_pred = Table("srx_metadata")
        stmt = Query \
            .from_(tbl_eval) \
            .where(tbl_eval.dataset_id == args.eval_dataset) \
            .join(tbl_pred) \
            .on(
                (tbl_eval.database == tbl_pred.database) & 
                (tbl_eval.entrez_id == tbl_pred.entrez_id) &
                (tbl_eval.srx_accession == tbl_pred.srx_accession)
            ) \
            .select("*") 
        df = pd.read_sql(str(stmt), conn).drop("id", axis=1)
        df.columns = add_suffix(df.columns, "_pred")
    
    eval(df)

    #compare_dfs(df, df, ["database", "entrez_id", "srx_accession"])



    #print(df.columns)

    # load predictions with overlapping database, entrez_id, and srx_accession columns
    # pred_data = None
    # with db_connect() as conn:
    #     tbl = Table("srx_metadata")
    #     stmt = Query \
    #         .from_(tbl) \
    #         .where(Criterion.all([
    #             tbl.database.isin(eval_data.database.unique()),
    #             tbl.entrez_id.isin(eval_data.entrez_id.unique()),
    #             tbl.srx_accession.isin(eval_data.srx_accession.unique())
    #         ])) \
    #         .select("*") 
    #     pred_data = pd.read_sql(str(stmt), conn)
    # print(pred_data)


# Example usage
if __name__ == "__main__":
    args = parser.parse_args()
    load_dotenv()
    main(args)