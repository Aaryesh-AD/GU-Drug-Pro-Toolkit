from rdkit import Chem
import pubchempy as pcp
import requests
import re

# Precompile the regex pattern for efficiency
cas_reg = re.compile(r'^(\d{1,7})-(\d{2})-(\d)$')
def is_valid_cas_number(cas_number):
    # Use the compiled regex object directly for matching
    match = cas_reg.match(cas_number)
    if not match:
        return False
    part1, part2, checksum = match.groups()
    calculated_checksum = sum((i + 1) * int(digit) for i, digit in enumerate(reversed(part1 + part2))) % 10
    return calculated_checksum == int(checksum)

def check_internet_connection():
    """Checks for an active internet connection."""
    try:
        response = requests.get('http://www.google.com', timeout=5)
        return response.status_code == 200
    except requests.ConnectionError:
        return False

def chem_name_info(smiles):
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:  # If RDKit cannot parse the SMILES string, return an error
        return {"Error": "Invalid SMILES string"}
    # Prepare initial parameters with InChi and InChiKey
    params = {
        "SMILES": smiles,
        "InChi": Chem.MolToInchi(mol),
        "InChiKey": Chem.MolToInchiKey(mol),
        
    }
    
    # Check internet connection before proceeding
    connection_status = check_internet_connection()
    if not connection_status:
        return params
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            compound = compounds[0]
            params["IUPAC Name"] = compound.iupac_name
            params["Synonyms"] = ', '.join(compound.synonyms[:10])
            cas_number = next((syn for syn in compound.synonyms if is_valid_cas_number(syn) or 'CAS' in syn), "Not found")
            params["CAS-Number"] = cas_number
        
    except Exception as e:
        # Log the exception details here, if possible

        params["Error"] = "Failed to retrieve additional data from PubChem"

    return params

