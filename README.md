Genomics Guide v2
=================

Chatbot for answering questions about genomics at the Arc Institute.

# Install 
    
```bash
pip install .
```




# Development

## Install

```bash
pip install -e .
```


# Agent structure

* Supervisor
* Workers
  * `Database query`
    * tools
      * `Entrez_esearch`
        * search terms
        * organism
        * Notes: the results must be parsed
      * `Entrez_efetch`
        * database
        * id
        * Notes: the results must be parsed
      * `geo2sra`
        * geo_id
  * `Critic`
    * tools
      * `web search for accession`?
      * `read abstract`?
