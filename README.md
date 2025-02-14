SRAgent
=======

LLM agents for working with the SRA and associated bioinformatics databases.



# Install

Create a conda environment [optional]: 

```bash
mamba create -n sragent-env -y python=3.12 sra-tools=3.1 \
  && conda activate sragent-env
```

Clone the repository:

```bash
git clone git@github.com:ArcInstitute/SRAgent.git \
  && cd SRAgent
```

Install the package:
    
```bash
pip install .
```

## Environmental variables

* `OPENAI_API_KEY` = API key for using the OpenAI API
  * **required**
* `EMAIL` = email for using the Entrez API
  * optional, but **HIGHLY** recommended
* `NCBI_API_KEY` = API key for using the Entrez API
  * optional, increases rate limits
* `DYNACONF` = switch between "test" and "prod" environments
  * optional, default is "prod"

# Usage

## find-datasets agent

Find datasets in the SRA database and then process them with the SRX-info agent. 

#### Examples

```bash
SRAgent find-datasets "Obtain recent single cell RNA-seq datasets in the SRA database"
```

#### Target specific organisms

```bash
SRAgent find-datasets --no-summaries --max-datasets 1 --organisms rat -- \
  "Obtain recent single cell RNA-seq datasets in the SRA database"
```

<details>
  <summary><strong>Available organisms</strong></summary>

- **Mammals**
  - Human (*Homo sapiens*)
  - Mouse (*Mus musculus*)
  - Rat (*Rattus norvegicus*)
  - Macaque (*Macaca mulatta*)
  - Marmoset (*Callithrix jacchus*)
  - Horse (*Equus caballus*)
  - Dog (*Canis lupus*)
  - Bovine (*Bos taurus*)
  - Sheep (*Ovis aries*)
  - Pig (*Sus scrofa*)
  - Rabbit (*Oryctolagus cuniculus*)
  - Naked mole-rat (*Heterocephalus glaber*)
  - Chimpanzee (*Pan troglodytes*)
  - Gorilla (*Gorilla gorilla*)

- **Birds**
  - Chicken (*Gallus gallus*)

- **Amphibians**
  - Frog (*Xenopus tropicalis*)

- **Fish**
  - Zebrafish (*Danio rerio*)

- **Invertebrates**
  - Fruit fly (*Drosophila melanogaster*)
  - Roundworm (*Caenorhabditis elegans*)
  - Mosquito (*Anopheles gambiae*)
  - Blood fluke (*Schistosoma mansoni*)

- **Plants**
  - Thale cress (*Arabidopsis thaliana*)
  - Rice (*Oryza sativa*)
  - Tomato (*Solanum lycopersicum*)
  - Corn (*Zea mays*)

</details>


#### Using an SQL database to store results

Using the `test` database:

```bash
SRAgent find-datasets --use-database --no-summaries --max-datasets 1 --organisms rat -- \
  "Obtain recent single cell RNA-seq datasets in the SRA database"
```


## SRX-info agent

Obtain SRX metadata for >=1 SRA or GEO dataset.

> The metadata is stored in the 
[SRAgent_database](https://docs.google.com/spreadsheets/d/1dkFvBYTX7DQLxLQKjwQxMvo5dx6fijSh4TRCX2xChlA/edit?usp=sharing)
Google Sheet by default.

#### Examples

A single SRA dataset:

```bash
SRAgent srx-info 25576380
```

Multiple SRA datasets:

```bash
SRAgent srx-info 36404865 36106630 32664033
```

Use the SQL database to filter existing:

```bash
SRAgent srx-info --use-database 18060880 27454880 27454942 27694586
```

## Metadata agent

Provide a CSV of Entrez IDs and their associated SRX accessions to obtain metadata. 
Useful for when you already have the SRX accessions, instead of providing the Entrez IDs to `SRAgent srx-info`.

The CSV should have the header: `entrez_id,srx_accession`.

#### Examples

```bash
SRAgent metadata "entrez-id_srx-accession.csv"
```

## SRAgent agent

General tool for working with the SRA database via Entrez tools, SRA BigQuery, fetching NCBI webpages, and other methods.

#### Example of converting a GEO accession to SRX accessions:

```bash
SRAgent sragent "Convert GSE121737 to SRX accessions"
```

#### Example of obtaining metadata for a specific SRX accession:

```bash
SRAgent sragent "Obtain any available publications for GSE196830"
```

#### Example of obtaining specific metadata fields for a dataset:

```bash
SRAgent sragent "Which 10X Genomics technology was used for ERX11887200?"
```

## Entrez Agent

General agent for working specifically with Entrez tools (esearch, efetch, esummary, elink).
Usually, the SRAgent agent will be more useful.

#### Example accession conversion:

```bash
SRAgent entrez "Convert GSE121737 to SRX accessions"
```

#### Example of obtaining pubmed articles associated with a dataset accession:

```bash
SRAgent entrez "Obtain any available publications for GSE196830"
```


# Evaluations

See the `eval.py` script for running evaluations.

# Workflows

* Obtain studes
  * esearch
* Convert to SRX
  * entrez
  * ncbi-fetch
  * sra-bigquery
* Get SRX metadata
  * entrez
  * ncbi-fetch
  * sra-bigquery
  * seq
* Get SRR accessions per SRX
  * entrez
  * ncbi-fetch
  * sra-bigquery


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


# resources

* https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/


***

# OLD

## Network proxy

Install via (assuming `${HOME}/bin` is in your path):

```bash
mkdir -p ${HOME}/bin/ \
  && curl -o ${HOME}/bin/cloud-sql-proxy \
    https://storage.googleapis.com/cloud-sql-connectors/cloud-sql-proxy/v2.14.1/cloud-sql-proxy.linux.amd64 \
  && chmod u+x ${HOME}/bin/cloud-sql-proxy \
  && mkdir -p ${HOME}/cloudsql
```

Run via:

```bash
SERVICE_ACCOUNT_JSON="c-tc-429521-6f6f5b8ccd93.json"
PROXY_NAME="c-tc-429521:us-east1:sragent"
rm -rf ${HOME}/cloudsql/${PROXY_NAME}
cloud-sql-proxy ${PROXY_NAME} \
  --unix-socket ${HOME}/cloudsql \
  --credentials-file ${HOME}/.gcp/${SERVICE_ACCOUNT_JSON}
```