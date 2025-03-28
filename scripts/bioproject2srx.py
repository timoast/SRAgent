#!/usr/bin/env python3
import sys
import argparse
from Bio import Entrez

def fetch_srx_accessions(bioproject_id, email):
    # set email
    Entrez.email = email

    # run esearch
    print(f"esearch of bioproject for: {bioproject_id}", file=sys.stderr)
    handle = Entrez.esearch(db="bioproject", term=bioproject_id)
    record = Entrez.read(handle)
    if not record["IdList"]:
        raise ValueError(f"No UID found for BioProject {bioproject_id}")
    bioproject_uid = record["IdList"][0]

    # run elink
    print(f"elink of bioproject for: {bioproject_id}", file=sys.stderr)
    handle = Entrez.elink(dbfrom="bioproject", id=bioproject_uid, db="sra")
    linkset = Entrez.read(handle)
    if not linkset[0]["LinkSetDb"]:
        raise ValueError(f"No linked SRA records found for {bioproject_id}")
    sra_ids = [link["Id"] for link in linkset[0]["LinkSetDb"][0]["Link"]]
    print(f"  Total SRA records: {len(sra_ids)}", file=sys.stderr)

    # fetch SRX accessions
    srx_list = set()
    for sra_id in sra_ids:
        print(f"efetch of sra for: {sra_id}", file=sys.stderr)
        handle = Entrez.efetch(db="sra", id=sra_id, rettype="runinfo", retmode="text")
        lines = handle.read().decode('utf-8').splitlines()
        # get header
        header = {x:i for i,x in enumerate(lines[0].split(","))}

        for line in lines[1:]:
            line = line.split(",")
            try:
                srx_list.add(line[header["Experiment"]])
            except KeyError:
                continue
    print(f"  Total SRX records: {len(srx_list)}", file=sys.stderr)
    return sorted(srx_list)

def main():
    parser = argparse.ArgumentParser(description="Fetch all SRX accessions from a BioProject")
    parser.add_argument("bioproject_id", help="NCBI BioProject accession (e.g., PRJNA123456)")
    parser.add_argument("--email", required=True, help="Your email for NCBI Entrez access")
    args = parser.parse_args()

    try:
        srx_accessions = fetch_srx_accessions(args.bioproject_id, args.email)
        for srx in srx_accessions:
            print(srx)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()