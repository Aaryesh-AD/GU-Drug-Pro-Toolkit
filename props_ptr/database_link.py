import urllib.parse, re
import webbrowser, urllib
import pubchempy as pcp
from rdkit import Chem

cas_reg = re.compile(r'^(\d{1,7})-(\d{2})-(\d)$')
def is_valid_cas_number(cas_number):
    # Use the compiled regex object directly for matching
    match = cas_reg.match(cas_number)
    if not match:
        return False
    part1, part2, checksum = match.groups()
    calculated_checksum = sum((i + 1) * int(digit) for i, digit in enumerate(reversed(part1 + part2))) % 10
    return calculated_checksum == int(checksum)

def get_database_links(smiles):
    # Encode the SMILES code for URL use
    encoded_smiles = urllib.parse.quote(smiles)
    mol = Chem.MolFromSmiles(smiles)
    inchikey = Chem.MolToInchiKey(mol)
    compounds = pcp.get_compounds(smiles, namespace='smiles')
    if compounds:
            compound = compounds[0]
            cas_number = next((syn for syn in compound.synonyms if is_valid_cas_number(syn) or 'CAS' in syn), None)
    # URLs for database searches, excluding ChemDB due to its specific requirements
    links = {
        "PubChem": f"https://pubchem.ncbi.nlm.nih.gov/#query={encoded_smiles}",
        "ZINC": f"http://zinc15.docking.org/substances/search/?q={encoded_smiles}",
        "ChEMBL": f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={encoded_smiles}",
        "ChemSpider": f"http://www.chemspider.com/Search.aspx?q={encoded_smiles}",
        "Drug Central": f"https://drugcentral.org/substructure?m={encoded_smiles}",
        "UniChem": f"https://www.ebi.ac.uk/unichem/compoundsources?type=inchikey&compound={inchikey}",
        "Drug Bank": f"https://go.drugbank.com/unearth/q?searcher=drugs&query={cas_number}&button="
    }

    return links
