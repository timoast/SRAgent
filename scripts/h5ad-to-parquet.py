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

# Configure logging with timestamps
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logging.getLogger("scanpy").setLevel(logging.WARNING)
logging.getLogger("anndata").setLevel(logging.WARNING)

def parse_cli_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments containing GCP path, feature type, and number of workers.
    """
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description='Convert h5ad files to combined parquet metadata per organism on GCP',
        epilog="""Example:
    ./h5ad-to-parquet.py gs://bucket/path/ --workers 4""",
        formatter_class=CustomFormatter
    )
    parser.add_argument(
        'gcp_path', type=str, help='GCP base path for input files and output metadata'
    )
    parser.add_argument(
        '--feature-type', default='GeneFull_Ex50pAS',
        choices=['Gene', 'GeneFull', 'GeneFull_Ex50pAS', 'GeneFull_ExonOverIntron', 'Velocyto'],
        help='Feature type to process'
    )
    parser.add_argument(
        '--max-files', type=int, default=0,
        help='Maximum number of files to process (0 for all files)'
    )
    parser.add_argument(
        '--workers', type=int, default=1,
        help='Number of workers for parallel processing per organism'
    )
    return parser.parse_args()

def parse_gcp_path(gcp_path: str) -> Tuple[str, str]:
    """
    Parse a GCP path to extract the bucket name and the blob name (prefix).
    
    Args:
        gcp_path (str): GCP path (e.g., gs://bucket/path)
    Returns:
        Tuple[str, str]: Bucket name and blob name.
    """
    parts = gcp_path.replace("gs://", "").split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix

def get_h5ad_files(gcp_path: str, feature_type: str) -> pd.DataFrame:
    """
    List all .h5ad.gz files in the specified GCP bucket path for a given feature type.
    
    Args:
        gcp_path (str): Base GCP path.
        feature_type (str): Feature type to filter by.
    Returns:
        pd.DataFrame: DataFrame with columns for feature_type, organism, and file_path.
    """
    print(f"Listing files in {gcp_path}...")
    bucket_name, prefix = parse_gcp_path(gcp_path)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    # Construct the full prefix based on the folder structure
    full_prefix = os.path.join(prefix, "h5ad", feature_type)
    blobs = bucket.list_blobs(prefix=full_prefix)
    h5ad_files = [
        f"gs://{bucket_name}/{blob.name}"
        for blob in blobs if blob.name.endswith(".h5ad.gz")
    ]

    def split_path(f: str) -> Tuple[str, str, str]:
        parts = f.split("/")
        # Assume folder structure: .../{feature_type}/{organism}/{filename}
        return parts[-3], parts[-2], f

    df = pd.DataFrame(
        [split_path(f) for f in h5ad_files],
        columns=["feature_type", "organism", "file_path"]
    )
    # Filter to ensure the feature type matches the target
    df = df[df["feature_type"] == feature_type]
    if not df["feature_type"].eq(feature_type).all():
        raise ValueError(f"Feature type mismatch: {df['feature_type'].unique()}")
    if df.empty:
        raise ValueError(f"No files found for feature type: {feature_type}")
    print(f"Found {len(df)} h5ad files")
    return df

def read_h5ad_from_gcs(file_path: str) -> pd.DataFrame:
    """
    Read an h5ad file from GCS and return its observation metadata as a DataFrame.
    
    Args:
        file_path (str): GCS path to the h5ad file.
    Returns:
        pd.DataFrame: Observation metadata.
    """
    bucket_name, blob_name = parse_gcp_path(file_path)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    file_bytes = blob.download_as_bytes()
    adata = sc.read_h5ad(io.BytesIO(file_bytes))
    df = adata.obs
    df["cell_barcode"] = df.index
    return df.reset_index(drop=True)

def process_org(organism: str, group: pd.DataFrame, gcp_path: str, feature_type: str, workers: int) -> None:
    """
    Process all h5ad files for a given organism in parallel, combine their observation metadata,
    and write the combined DataFrame as a parquet file to GCS.
    
    Args:
        organism (str): Organism name.
        group (pd.DataFrame): DataFrame containing file paths for the organism.
        gcp_path (str): Base GCP path for outputs.
        feature_type (str): Feature type being processed.
        workers (int): Number of workers for parallel file processing.
    """
    logging.info(f"Processing organism: {organism} (file count: {group.shape[0]})")
    file_paths: List[str] = group["file_path"].tolist()
    # Read each file in parallel
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            result_list = list(executor.map(read_h5ad_from_gcs, file_paths))
    else:
        result_list = [read_h5ad_from_gcs(fp) for fp in file_paths]
    # Combine all observation metadata DataFrames
    combined_df = pd.concat(result_list, ignore_index=True)
    # Define output file path: gs://bucket/path/metadata/{feature_type}/{organism}/obs_metadata.parquet.gz
    outfile = os.path.join(gcp_path, "metadata", feature_type, organism, "obs_metadata.parquet.gz")
    bucket_name, blob_name = parse_gcp_path(outfile)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    parquet_bytes = combined_df.to_parquet(None, compression="gzip")
    blob.upload_from_string(parquet_bytes, content_type="application/octet-stream")
    logging.info(f"Metadata written for organism {organism} to {outfile}")

def main(args: argparse.Namespace) -> None:
    """
    Main function to list h5ad files, group them by organism, process files in parallel per organism,
    combine the results, and write the combined metadata to GCS.
    
    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    # read h5ad files
    h5ad_df = get_h5ad_files(args.gcp_path, args.feature_type)

    # filter
    if args.max_files > 0:
        h5ad_df = h5ad_df.head(args.max_files)

    # group by organism and process in parallel
    for organism, group in h5ad_df.groupby("organism"):
        process_org(organism, group, args.gcp_path, args.feature_type, args.workers)

if __name__ == "__main__":
    load_dotenv(override=True)
    multiprocessing.set_start_method("spawn")  # Use "spawn" to avoid forking issues
    args = parse_cli_args()
    main(args)