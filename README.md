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

> The metadata is stored in the 
[SRAgent_database](https://docs.google.com/spreadsheets/d/1dkFvBYTX7DQLxLQKjwQxMvo5dx6fijSh4TRCX2xChlA/edit?usp=sharing)
Google Sheet by default.

Example of querying metadata for an SRA dataset (Entrez ID 36178506):

```bash
SRAgent metadata 36178506
```

Example of querying metadata for a GEO dataset (Entrez ID 200254051):

```bash
SRAgent metadata --database gds 200254051
```

# About

## Tools

The following tools are available for interacting with NCBI databases:

### esearch

Search NCBI databases using query terms:
* Search for specific accessions or terms across databases (sra, gds, pubmed)
* Specialized search for recent single-cell RNA-seq studies
* Returns Entrez IDs for matching records

### efetch

Fetch detailed metadata records:
* Retrieve full metadata for specific Entrez IDs
* Supports multiple databases (sra, gds, pubmed)
* Returns detailed XML/JSON format records

### esummary

Get summary information:
* Retrieve concise summaries for specific Entrez IDs
* Supports multiple databases (sra, gds, pubmed)
* Returns summarized record information

### elink

Find related records across databases:
* Link records between different NCBI databases
* Find associated BioProject, BioSample, or publication records
* Map relationships between different types of records

### ncbi_fetch

Direct web scraping of NCBI pages:
* Fetch detailed information from SRA, GEO, and PubMed web pages
* Extract structured data from HTML responses
* Useful for getting human-readable descriptions

### seq

Tools for working with sequence data:
* Use fastq-dump to preview FASTQ file contents
* Get sequence statistics using sra-stat
* Validate sequence data format and quality
* Check paired-end vs single-end status

## Agents

### Entrez Agent

A ReAct agent that coordinates NCBI database queries using the available tools:
* Converts between different accession types (GEO, SRA, BioProject)
* Retrieves metadata from various NCBI databases
* Follows multi-step workflows to gather comprehensive information
* Handles rate limits and batches large queries

### Convert Agent

Specialized agent for converting between different accession types:
* Focuses on obtaining SRX accessions from other identifiers
* Works with Entrez IDs, GEO accessions, and BioProject IDs
* Validates accession formats and handles edge cases
* Uses retry logic when conversions require multiple steps

### Metadata Agent

LangGraph workflow for extracting standardized metadata:
* Determines sequencing platform (Illumina vs other)
* Identifies single-cell vs bulk RNA-seq protocols
* Validates paired-end vs single-end sequencing
* Detects 10X Genomics library preparation
* Maps organism taxonomy
* Supports both SRA and GEO databases
* Uses multiple approaches to resolve uncertain metadata

## Workflows

### Metadata Workflow

Multi-stage workflow for processing sequencing datasets:
* Converts database records to SRA accessions (SRX/ERX)
* Processes each accession in parallel using the Metadata Agent
* Extracts standardized metadata fields for each sample
* Validates and consolidates results across all samples
* Optionally stores results in a tracking database
* Handles both SRA and GEO database records
* Supports batched processing of large datasets

