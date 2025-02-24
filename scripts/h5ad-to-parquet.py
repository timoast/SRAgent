#!/usr/bin/env python
import os
import io
import argparse
import logging
from typing import List, Tuple
from dotenv import load_dotenv
from google.cloud import storage
import scanpy as sc
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

# logging with timestamp
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logging.getLogger("scanpy").setLevel(logging.WARNING)
logging.getLogger("anndata").setLevel(logging.WARNING)


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


def parse_gcp_path(gcp_path: str) -> Tuple[str, str]:
    """
    Parse the GCP path to extract the bucket name and prefix.

    Args:
        gcp_path (str): GCP bucket path.
    Returns:
        Tuple[str, str]: Bucket name and prefix.
    """
    parts = gcp_path.replace("gs://", "").split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix


def get_h5ad_files(gcp_path: str, feature_type: str) -> pd.DataFrame:
    """
    List all h5ad.gz files in the specified GCP bucket path.

    Args:
        gcp_path (str): GCP bucket path.
        feature_type (str): Feature type to filter by.
    Returns:
        pd.DataFrame: DataFrame containing organism and file paths.
    """
    print(f"Listing files in {gcp_path}...")
    bucket_name, prefix = parse_gcp_path(gcp_path)
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # Adjust the prefix to match your folder structure
    full_prefix = os.path.join(prefix, "h5ad", feature_type)
    blobs = bucket.list_blobs(prefix=full_prefix)

    h5ad_files = [
        f"gs://{bucket_name}/{blob.name}"
        for blob in blobs
        if blob.name.endswith(".h5ad.gz")
    ]

    def split_path(f: str) -> Tuple[str, str, str]:
        p = f.split("/")
        return p[-3], p[-2], f

    df = pd.DataFrame(
        [split_path(f) for f in h5ad_files],
        columns=["feature_type", "organism", "file_path"]
    )
    
    # filter to target feature type
    df = df[df["feature_type"] == feature_type]
    if not df["feature_type"].eq(feature_type).all():
        raise ValueError(f"Feature type mismatch: {df['feature_type'].unique()}")

    ## check empty
    if df.empty:
        raise ValueError(f"No files found for feature type: {feature_type}")

    print(f"Found {len(df)} h5ad files")
    return df


def read_h5ad_from_gcs(file_path: str) -> pd.DataFrame:
    """
    Read h5ad file from GCS and return its observation metadata.

    Args:
        file_path (str): GCS path to the h5ad file.
    Returns:
        pd.DataFrame: Observation metadata as a DataFrame.
    """
    bucket_name, blob_name = parse_gcp_path(file_path)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)

    # Read h5ad file into memory
    file_bytes = blob.download_as_bytes()
    df = sc.read_h5ad(io.BytesIO(file_bytes)).obs
    df["cell_barcode"] = df.index
    return df.reset_index(drop=True)


def process_org(organism: str, group: pd.DataFrame, gcp_path: str, feature_type: str) -> None:
    """
    Process a group of h5ad files for a specific organism.

    Args:
        organism (str): Organism name.
        group (pd.DataFrame): DataFrame containing file paths for the organism.
        gcp_path (str): GCP bucket path.
    """
    obs_metadata = []
    for idx, row in group.iterrows():
        logging.info(f"  Reading: {row['file_path']}")
        obs_metadata.append(read_h5ad_from_gcs(row["file_path"]))

        # Debugging: Limit to 3 files per organism (Remove in production)
        #if idx >= 2:
        #    break

    # Combine and write to parquet format in GCS bucket
    outfile = os.path.join(gcp_path, "metadata", feature_type, organism, "obs_metadata.parquet.gz")
    bucket_name, blob_name = parse_gcp_path(outfile)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    blob.upload_from_string(
        pd.concat(obs_metadata).to_parquet(None, compression="gzip"),
        content_type="application/octet-stream"
    )
    logging.info(f"  Metadata written to {outfile}")


def process_group(args_tuple: Tuple[str, pd.DataFrame, str, str]) -> None:
    """
    Wrapper for processing groups in parallel.

    Args:
        args_tuple (Tuple[str, pd.DataFrame, str, str]): Tuple containing organism, group DataFrame, GCP path, and feature type.
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
    h5ad_files = get_h5ad_files(args.gcp_path, args.feature_type)

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
    multiprocessing.set_start_method("spawn")  # Use "spawn" to avoid forking issues
    args = parse_cli_args()
    main(args)