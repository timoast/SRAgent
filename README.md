SRAgent
=======

Agentic workflows for obtaining data from the Sequence Read Archive.

# Manuscript

[scBaseCamp: an AI agent-curated, uniformly processed, and continually expanding single cell data repository](https://arc-website-git-ben-virtual-cell-atlas-tool-arc-institute.vercel.app/manuscripts/scBaseCamp)

# Install

Create a conda environment [optional]: 

```bash
mamba create -n sragent-env -y python=3.12 sra-tools=3.1 \
  && conda activate sragent-env
```

Clone the repository:

```bash
git clone https://github.com/ArcInstitute/SRAgent.git \
  && cd SRAgent
```

Install the package:
    
```bash
pip install .
```

## Environmental variables

* `OPENAI_API_KEY` = API key for using the OpenAI API
  * **required**
  * currently, no other models are supported besides OpenAI
* `EMAIL` = email for using the Entrez API
  * optional, but **HIGHLY** recommended
* `NCBI_API_KEY` = API key for using the Entrez API
  * optional, increases rate limits
* `DYNACONF` = switch between "test" and "prod" environments
  * optional, default is "prod"
  * this only affects the SQL database used, and no database is used by default

# Testing

```bash
pip install pytest
```

```bash
pytest tests/
```

# Usage

## SQL database

Components of SRAgent can use an SQL database to store the results.

This was crucial for the scBaseCamp project, in order to:
* track which datasets had been processed 
* quickly assess the progress of the project

However, for most users, the SQL database is not necessary.
SRAgent does not use the SQL database by default.

> Note: currently only a GCP Postgresql database is supported.

To set up the database, see [Setting up the SQL Database](#setting-up-the-sql-database).

## Entrez Agent

The lowest-level agent in the SRAgent hiearchy.
The agent can call various Entrez tools (`esearch`, `efetch`, `esummary`, and `elink`).
Usually, the SRAgent agent will be more useful, since it includes more tools, including calling the Entrez agent.

#### Example accession conversion:

```bash
SRAgent entrez "Convert GSE121737 to SRX accessions"
```

#### Example of obtaining pubmed articles associated with a dataset accession:

```bash
SRAgent entrez "Obtain any available publications for GSE196830"
```

## SRAgent agent

A general tool for extracting data from the SRA database.
The tools available:

* Entrez agent (see above)
* SRA BigQuery
* scraping NCBI webpage HTML
* sra-stat and fastq-dump (directly assessing sequence data)

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

## SRX-info agent

Obtain specific metadata for >=1 SRA dataset.

* Input: >=1 Entrez ID
* Output metadata fields:
  * SRX accession for the Entrez ID
  * SRR accessions for the SRX accession
  * Is the dataset Illumina sequence data?
  * Is the dataset single cell RNA-seq data?
  * Is the dataset paired-end sequencing data?
  * Which scRNA-seq library preparation technology?
  * If 10X Genomics, which particular 10X technologies?
  * Single nucleus or single cell RNA sequencing?
  * Which organism was sequenced?
  * Which tissue was sequenced?
  * Any disease information?
  * Any treatment/purturbation information?
  * Any cell line information?
* Workflow
  * The agent converts the Entrez IDs to SRX accessions
  * For each SRX accession, the agent obtains metadata
  * The agent consolidates the metadata into a single report

> As of now, the metadata fields are hard-coded into the agent.
> If you need alternative metadata fields, you will have to modify [metadata.py](./workflows/metadata.py).

#### Examples

A single SRA dataset:

```bash
SRAgent srx-info 25576380
```

Multiple SRA datasets:

```bash
SRAgent srx-info 36404865 36106630 32664033
```

Use the SQL database to filter out already-processed datasets:

```bash
SRAgent srx-info --use-database 18060880 27454880 27454942 27694586
```


## Metadata agent

Similar to the `SRX-info` agent, but you can provide SRX accessions directly, instead of Entrez IDs.
This saves compute time, since the agent does not need to convert the Entrez IDs to SRX accessions.

Provide a CSV of Entrez IDs and their associated SRX accessions to obtain metadata. 
Useful for when you already have the SRX accessions, instead of providing the Entrez IDs to `SRAgent srx-info`.

The CSV should have the header: `entrez_id,srx_accession`.

The metadata fields are the same as the `SRX-info` agent.

#### Examples

```bash
SRAgent metadata "entrez-id_srx-accession.csv"
```

## find-datasets agent

A high-level agent for finding datasets in the SRA via `esearch` and then
processing them with the `SRX-info` agent. 

* Input: a search query
* Output: metadata fields for the datasets found (same as `SRX-info` agent)
* Workflow
  * The agent uses `esearch` to find datasets
  * The agent processes the datasets with the `SRX-info` agent
  * The agent consolidates the metadata into a single report

#### Examples

```bash
SRAgent find-datasets "Obtain recent single cell RNA-seq datasets in the SRA database"
```

#### Target specific organisms

```bash
SRAgent find-datasets --no-summaries --max-datasets 1 --organisms pig -- \
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

# Setting up the SQL database

* Create a GCP Postgresql database. See the [docs](https://cloud.google.com/sql?hl=en).
* Required secrets:
  * `GCP_SQL_DB_PASSWORD`
  * Store in the `.env` file or GCP Secret Manager
  * If using GCP Secret Manager, you must also provide:
    * `GOOGLE_APPLICATION_CREDENTIALS`
    * `GCP_PROJECT_ID`
* Update the [settings.py](./SRAgent/settings.yml) file with the database information.

# Evaluations

See the [eval.py](./scripts/eval.py) script for running evaluations.


# Contributing

Feel free to fork the repository and submit a pull request.
However, the top priority is to keep SRAgent functioning 
for the ongoing scBaseCamp project.