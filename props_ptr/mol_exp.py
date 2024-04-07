from rdkit import Chem
import pubchempy as pcp
import requests

def check_internet_connection():
    """Checks for an active internet connection."""
    try:
        response = requests.get('http://www.google.com', timeout=5)
        return response.status_code == 200
    except requests.ConnectionError:
        return False
    

list_prop = [
    'elements',
    'atoms',
    'bonds',
    'molecular_formula',
    'molecular_weight',
    'exact_mass',
    'monoisotopic_mass',
    'h_bond_acceptor_count',
    'h_bond_donor_count',
    'heavy_atom_count',
    'isotope_atom_count',
    'rotatable_bond_count',
    'atom_stereo_count',
    'bond_stereo_count',
    'covalent_unit_count',
    'tpsa',
    'xlogp',
    'volume_3d',
    'charge',
    'complexity',
    'coordinate_type',
    'defined_atom_stereo_count',
    'defined_bond_stereo_count',
    'undefined_atom_stereo_count',
    'undefined_bond_stereo_count',
    'effective_rotor_count_3d',
    'conformer_id_3d',
    'conformer_rmsd_3d',
    'feature_selfoverlap_3d',
    'mmff94_energy_3d',
    'mmff94_partial_charges_3d',
    'multipoles_3d',
    'pharmacophore_features_3d',
    'shape_fingerprint_3d',
    'shape_selfoverlap_3d'
]

def exp_info(smiles):
    params = {}  # Initialize params dictionary
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:  # If RDKit cannot parse the SMILES string, return an error
        return {"Error": "Invalid SMILES string"}

    # Check internet connection before proceeding
    connection_status = check_internet_connection()
    if not connection_status:
        params["Error"] = "Failed to retrieve additional data from PubChem"
        return params
    
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            compound = compounds[0]
            for prop in list_prop:
                params[prop] = getattr(compound, prop, 'N/A')  # Fetch each property
    except Exception as e:
        params["Error"] = "Failed to retrieve additional data from PubChem"

    return params


