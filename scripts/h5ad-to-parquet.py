#!/usr/bin/env python
import os
import io
import argparse
import logging
from typing import List, Tuple
from dotenv import load_dotenv
import gcsfs
import scanpy as sc
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# logging with timestamp
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logging.getLogger("scanpy").setLevel(logging.WARNING)   
logging.getLogger("anndata").setLevel(logging.WARNING)
logging.getLogger("gcsfs").setLevel(logging.WARNING)

# Functions
def parse_cli_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments containing GCP path and number of workers.
    """
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description='Database tools',
        epilog="""DESCRIPTION:
        Example:
        ./h5ad-to-parquet.py gs://arc-ctc-scbasecamp/2025-02-25/ --workers 4""",
        formatter_class=CustomFormatter
    )
    parser.add_argument(
        'gcp_path', type=str, help='GCP path to h5ad files'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS', 
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'], 
        help='Feature type to process'
    )
    parser.add_argument(
        '--workers', type=int, default=1, help='Number of workers for parallel processing'
    )
    return parser.parse_args()

def get_h5ad_files(gcp_path: str) -> pd.DataFrame:
    """
    List all h5ad.gz files in the specified GCP bucket path.

    Args:
        gcp_path (str): GCP bucket path.
    Returns:
        pd.DataFrame: DataFrame containing organism and file paths.
    """
    print(f"Listing files in {gcp_path}...")
    fs = gcsfs.GCSFileSystem()

    # List all files in the bucket
    all_files = fs.glob(os.path.join(gcp_path, "**"))

    # Filter files with the .h5ad.gz extension
    h5ad_files = [f"gs://{f}" for f in all_files if f.endswith(".h5ad.gz")]

    # Convert to dataframe: organism & file_path
    df = pd.DataFrame(
        [[os.path.basename(os.path.dirname(f)), f] for f in h5ad_files],
        columns=["organism", "file_path"]
    )

    # Display the number of found files
    print(f"Found {len(df)} h5ad files")
    return df

def process_org(organism: str, group: pd.DataFrame, gcp_path: str, feature_type: str) -> None:
    """
    Process a group of h5ad files for a specific organism.

    Args:
        organism (str): Organism name.
        group (pd.DataFrame): DataFrame containing file paths for the organism.
        gcp_path (str): GCP bucket path.
    """
    # initialize GCS connection
    fs = gcsfs.GCSFileSystem()
    
    # Process each file in the group
    obs_metadata = []
    for idx, row in group.iterrows():
        logging.info(f"  Reading: {row['file_path']}")
        with fs.open(row["file_path"], "rb") as f:
            df = sc.read_h5ad(io.BytesIO(f.read())).obs
            df["cell_barcode"] = df.index
            obs_metadata.append(df.reset_index(drop=True))
        
        # Debugging: Limit to 3 files per organism (Remove in production)
        #if idx >= 2:
        #    break

    # Combine and write to parquet format in GCS bucket
    outfile = os.path.join(gcp_path, "metadata", organism, feature_type, "obs_metadata.parquet.gz")
    with fs.open(outfile, "wb") as f:
        f.write(pd.concat(obs_metadata).to_parquet(None, compression="gzip"))
    logging.info(f"  Metadata written to {outfile}")

def process_group(args_tuple: Tuple[str, pd.DataFrame, str]) -> None:
    """
    Wrapper for processing groups in parallel.

    Args:
        args_tuple (Tuple[str, pd.DataFrame, str]): Tuple containing organism, group DataFrame, and GCP path.
    """
    organism, group, gcp_path, feature_type = args_tuple
    logging.info(f"Processing {organism} (file count: {group.shape[0]})...")
    process_org(organism, group, gcp_path, feature_type)

def main(args: argparse.Namespace) -> None:
    """
    Main function to handle GCS connection, file listing, and parallel processing.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    # Get h5ad files
    h5ad_files = get_h5ad_files(args.gcp_path)

    # Group by organism and prepare for parallel processing
    groups = [(org, group, args.gcp_path, args.feature_type) for org, group in h5ad_files.groupby("organism")]

    # Process each group 
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            executor.map(process_group, groups)
    else:
        for group in groups:
            process_group(group)

if __name__ == "__main__":
    load_dotenv(override=True)
    args = parse_cli_args()
    main(args)