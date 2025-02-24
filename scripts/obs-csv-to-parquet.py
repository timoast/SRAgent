#!/usr/bin/env python
import os
import io
import argparse
import logging
from typing import List, Tuple
import pandas as pd
from dotenv import load_dotenv
from google.cloud import storage
from concurrent.futures import ProcessPoolExecutor

# Set up logging with timestamp
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

def parse_cli_args() -> argparse.Namespace:
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description="Database tools",
        epilog="""DESCRIPTION:
Example:
./obs-csv-to-parquet.py gs://arc-ctc-scbasecamp/2025-02-25/metadata_TMP/ --workers 4
""",
        formatter_class=CustomFormatter
    )
    parser.add_argument("gcp_path", type=str, help="GCP path to h5ad files")
    parser.add_argument(
        "--feature-type",
        default="GeneFull_Ex50pAS",
        choices=["Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "Velocyto"],
        help="Feature type to process"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="gs://arc-ctc-scbasecamp/2025-02-25/",
        help="Output directory for parquet files on GCS"
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Limit the number of files to process (0 for all)"
    )
    parser.add_argument(
        "--workers", type=int, default=1, help="Number of workers for parallel processing"
    )
    return parser.parse_args()

def parse_gcp_path(gcp_path: str) -> Tuple[str, str]:
    """Parse a GCS path into bucket name and prefix."""
    parts = gcp_path.replace("gs://", "").split("/", 1)
    bucket_name = parts[0]
    prefix = parts[1] if len(parts) > 1 else ""
    return bucket_name, prefix

def get_csv_files(gcp_path: str, feature_type: str) -> List[str]:
    """List all csv.gz files from a given GCS path and feature type."""
    logging.info(f"Listing files in {gcp_path}...")
    bucket_name, prefix = parse_gcp_path(gcp_path)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    full_prefix = os.path.join(prefix, feature_type)
    blobs = bucket.list_blobs(prefix=full_prefix)
    csv_files = [f"gs://{bucket_name}/{blob.name}" for blob in blobs if blob.name.endswith(".csv.gz")]
    return csv_files

def read_csv_file(gcs_path: str) -> pd.DataFrame:
    """Read a gzipped csv file from GCS and return a DataFrame."""
    logging.info(f"Reading {gcs_path}...")
    bucket_name, blob_name = parse_gcp_path(gcs_path)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    content = blob.download_as_string()
    return pd.read_csv(io.BytesIO(content), compression="gzip")

def save_parquet_to_gcs(df: pd.DataFrame, outdir: str, feature_type: str, organism: str) -> None:
    """
    Write the DataFrame to a parquet file and upload it to GCS.
    The file is saved under the path:
      <outdir>/<feature_type>/<organism>/obs_metadata.parquet.gz
    """
    organism_str = organism.replace(" ", "_")
    bucket_name, prefix = parse_gcp_path(outdir)
    blob_name = os.path.join(prefix, "metadata", feature_type, organism_str, "obs_metadata.parquet.gz")
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)
    # Write DataFrame to an in-memory buffer in parquet format with gzip compression
    buffer = io.BytesIO()
    df.drop(columns=["organism"]).to_parquet(buffer, index=False, compression="gzip")
    buffer.seek(0)
    blob.upload_from_file(buffer, content_type="application/octet-stream")
    logging.info(f"Saved metadata for {organism} to gs://{bucket_name}/{blob_name}")

def main(args: argparse.Namespace) -> None:
    # List all csv files to process
    csv_files = get_csv_files(args.gcp_path, args.feature_type)
    logging.info(f"Found {len(csv_files)} csv files to process.")
    
    # Limit the number of files to process if specified
    if args.max_files > 0:
        csv_files = csv_files[:args.max]

    # Read CSV files, using parallel processing if more than one worker is specified
    if args.workers == 1:
        dataframes = [read_csv_file(f) for f in csv_files]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(read_csv_file, f) for f in csv_files]
            dataframes = [future.result() for future in futures]
    metadata = pd.concat(dataframes, ignore_index=True)
    # Group metadata by organism and save each group as a parquet file to GCS
    for organism, df in metadata.groupby("organism"):
        logging.info(f"Processing metadata for {organism}...")
        save_parquet_to_gcs(df, args.outdir, args.feature_type, organism)

if __name__ == "__main__":
    load_dotenv(override=True)
    args = parse_cli_args()
    main(args)