SRAgent
=======

LLM agents for working with the SRA and associated bioinformatics databases.


# Install 
    
```bash
pip install .
```

# Development

## Install

```bash
pip install -e .
```

# Usage

## Entrez Agent

Example

```bash
SRAgent entrez-agent "Convert GSE121737 to SRX accessions"
```


# TODO

* [X] Handle conversion of GEO to SRA
  * For use of `geo2sra` if gds database
* [ ] Handle conversion of all final accessions to SRR
  * Adapt from https://github.com/ArcInstitute/scRecounter/blob/main/scripts/acc2srr.py
