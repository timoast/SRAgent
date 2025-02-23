import pytest
from Bio import Entrez
from SRAgent.tools.esearch import esearch

def test_esearch_successful():
    """Test successful esearch with actual API call"""
    input = {"esearch_query" : "35447314", "database" : "sra"}
    ret = esearch.invoke(input)
    
    assert "'Count': '1'" in ret
    assert "'IdList': ['35447314']" in ret
    
def test_esearch_invalid():
    """Test successful esearch with actual API call"""
    input = {"esearch_query" : "35447314", "database" : "geo"}
    ret = esearch.invoke(input)
    assert "Error searching" in ret
    