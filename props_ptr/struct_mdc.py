from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

# Define the function to get PAINS matches
def get_pains_matches(mol):
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog.FilterCatalog(params)
    
    matches = []
    for entry in pains_catalog.GetFilterMatches(mol):
        pattern_name = entry.GetDescription()
        matches.append(pattern_name)
    
    return matches

def check_leadlikeness(mol):
    mw = round(Descriptors.MolWt(mol),3)
    logp = round(Crippen.MolLogP(mol),3)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    violations = []
    
    if not (250 <= mw <= 350):
        violations.append(f"Molecular weight {mw} is out of range (250-350) ")
    if not (logp <= 3.5):
        violations.append(f"LogP {logp} is greater than 3.5 ")
    if not (rotatable_bonds <= 7):
        violations.append(f"Rotatable bonds {rotatable_bonds} is greater than 7")
    
    return violations

def struct_exp(smiles):
    results = {}
    mol = Chem.MolFromSmiles(smiles)
    pains = get_pains_matches(mol)
    
    if pains:
        results['PAINS Alerts'] = ', '.join(pains)
    else:
        results['PAINS Alerts'] = "None"
    
    lead_likeness = check_leadlikeness(mol)
    
    if not lead_likeness:  # If the list is empty, it passes the criteria
        results['Lead Likeness'] = "Pass - Satisfies all Teague criteria"
        
    else:
        results['Lead Likeness'] = "Fail - (Violations): "+','.join(lead_likeness)
        

    return results

