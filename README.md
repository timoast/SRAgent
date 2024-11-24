SRAgent
=======

LLM agents for working with the SRA and associated bioinformatics databases.


# Install 
    
```bash
pip install .
```

## Environmental variables

* `OPENAI_API_KEY` = API key for using the OpenAI API
  * **required**
* `EMAIL` = email for using the Entrez API
  * optional, but HIGHLY recommended
* `NCBI_API_KEY` = API key for using the Entrez API
  * optional, increases rate limits

# Development

## Install

```bash
pip install -e .
```

# Usage

## Entrez Agent

Example accession conversion:

```bash
SRAgent entrez "Convert GSE121737 to SRX accessions"
```

Example of obtaining pubmed articles associated with a dataset accession:

```bash
SRAgent entrez "Obtain any available publications for GSE196830"
```

## Metadata agent

Example of querying metadata for an SRA dataset (Entrez ID 36178506):

```bash
SRAgent metadata 36178506
```

Example of querying metadata for a GEO dataset (Entrez ID 200254051):

```bash
SRAgent metadata --database gds 200254051
```


# Workflow

## Option 1

* `esearch` to find Entrez IDs for new datasets
* `elink` for each GEO ID, find associated SRA IDs
* Merge Entrez IDs
* Filter IDs that have been processed
  * Database lookup
* Filter IDs that appear to not be scRNA-seq
  * Entrez-agent => `"Is Entrez ID {ID} a scRNA-seq dataset?"`
* For each remaining Entrez ID, obtain SRX accessions
* For each SRX accession:
  * Sequencing platform metadata:
    * `Illumina sequencing?`
    * `Paired-end?`
    * `10X Genomics library prep?`
    * `Species`
  * Obtain SRR accessions
    * For each SRR accesion:
      * Check sequence data files

## Option 2

Agents all the way down...

* Top agent
  * Input: dataset entrez ID
  * Tools
    * `convert`
    * `metadata`
    * `add2db`

  

# TODO

* [X] Handle conversion of GEO to SRA
  * For use of `geo2sra` if gds database
* [ ] Handle conversion of all final accessions to SRR
  * Adapt from https://github.com/ArcInstitute/scRecounter/blob/main/scripts/acc2srr.py
