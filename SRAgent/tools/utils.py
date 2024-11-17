# import
## batteries
import json
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any
import xml.etree.ElementTree as ET
from xml.parsers.expat import ExpatError
## 3rd party
import xmltodict

# functions
def batch_ids(ids: List[str], batch_size: int) -> List[List[str]]:
    """
    Batch a list of IDs into smaller lists of a given size.
    """
    for i in range(0, len(ids), batch_size):
        yield ids[i:i + batch_size]

def truncate_values(record, max_length):
    # truncate long values in the record
    try:
        root = ET.fromstring(record)
    except ET.ParseError:
        return record
    for item in root.findall(".//Item"):
        if item.text and len(item.text) > max_length:
            item.text = item.text[:max_length] + "...[truncated]"
    # convert back to string
    return ET.tostring(root, encoding="unicode")

def xml2json(record: str) -> Dict[str, Any]:
    """
    Convert an XML record to a JSON object.
    """
    try:
        return json.dumps(xmltodict.parse(record), indent=2)
    except ExpatError:
        return record