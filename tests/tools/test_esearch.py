import pytest
from Bio import Entrez
from SRAgent.tools.esearch import esearch
from SRAgent.tools.esearch import to_sci_name
from SRAgent.organisms import OrganismEnum

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

def test_to_sci_name_with_all_organisms():
    """Test to_sci_name function with all values from OrganismEnum"""
    # Test each enum value by converting its name to lowercase and checking the result
    for organism in OrganismEnum:
        # Get the enum name (e.g., 'HUMAN') and convert to lowercase
        organism_name = organism.name.lower()
        # Call to_sci_name with the lowercase name
        result = to_sci_name(organism_name)
        # Verify the result is the quoted scientific name from the enum
        expected = f'"{organism.value}"'
        assert result == expected, f"Failed for {organism_name}: expected {expected}, got {result}"
    
    # Test with spaces in names (replacing underscores with spaces)
    for organism in OrganismEnum:
        if "_" in organism.name:
            # Convert enum name to lowercase and replace underscores with spaces
            organism_name = organism.name.lower().replace("_", " ")
            # Call to_sci_name with the spaced name
            result = to_sci_name(organism_name)
            # Verify the result
            expected = f'"{organism.value}"'
            assert result == expected, f"Failed for '{organism_name}': expected {expected}, got {result}"

def test_to_sci_name_unknown_organism():
    """Test to_sci_name function with an unknown organism name"""
    with pytest.raises(ValueError, match="Organism 'nonexistent' not found in OrganismEnum"):
        to_sci_name("nonexistent")
