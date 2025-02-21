#!/usr/bin/env python
# import
## batteries
import os
import io
import sys
import argparse
from typing import List, Dict, Literal, Any
## 3rd party
from dotenv import load_dotenv
import gcsfs
import scanpy as sc
import pandas as pd

# functions
def parse_cli_args():
    # argparse
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.RawDescriptionHelpFormatter):
        pass

    desc = 'Database tools'
    epi = """DESCRIPTION:
    Example:
    ./h5ad-to-parquet.py gs://arc-ctc-scbasecamp/2025-02-25/
    """
    parser = argparse.ArgumentParser(
        description=desc, epilog=epi, formatter_class=CustomFormatter
    )
    parser.add_argument(
        'gcp_path', type=str, help='GCP path to h5ad file'
    )
    return parser.parse_args()

# functions
def get_h5ad_files(gcp_path, fs):
    print(f"Listing files in {gcp_path}...")

    # List all files in the bucket
    all_files = fs.glob(os.path.join(gcp_path, "**"))

    # Filter files with the specified extension
    h5ad_files = [file for file in all_files if file.endswith(".h5ad.gz")]

    # Convert to dataframe: organism & file_path
    df = pd.DataFrame(
        [[os.path.basename(os.path.dirname(file)),file] for file in h5ad_files], 
        columns=["organism", "file_path"]
    )
    # status
    print(f"Found {len(df)} h5ad files")
    return df

def process_org(organism: str, group: pd.DataFrame, fs, gcp_path: str):
    obs_metadata = []
    for idx, row in group.iterrows():
        # load h5ad file
        with fs.open(row["file_path"], "rb") as f:
            df = sc.read_h5ad(f).obs
            df["cell_barcode"] = df.index
            obs_metadata.append(df)

    # write to parquet to GCS bucket
    outfile = os.path.join(gcp_path, "metadata", organism, "obs_metadata.parquet.gz")
    with fs.open(outfile, 'wb') as f:
        f.write(pd.concat(obs_metadata).to_parquet(None, compression="gzip"))
    print(f"Written: {outfile}")

def main(args):
    fs = gcsfs.GCSFileSystem()

    # get h5ad files
    h5ad_files = get_h5ad_files(args.gcp_path, fs)

    # group by organism and process each group
    for organism, group in h5ad_files.groupby("organism"):
        print(f"Processing {organism} (file count: {group.shape[0]})...")
        process_org(organism, group, fs, args.gcp_path)
            


# Example usage
if __name__ == "__main__":
    load_dotenv(override=True)
    args = parse_cli_args()
    main(args)