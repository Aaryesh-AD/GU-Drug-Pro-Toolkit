from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from collections import Counter

def basic_mol_props(smiles):
    # Convert SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    mol_hyd = Chem.AddHs(mol)

    # Calculating molecular formula and weights
    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
    mol_wt = round(Descriptors.MolWt(mol), 4)
    exact_mol_wt = Descriptors.ExactMolWt(mol)
    
    # Calculate elemental composition
    element_counts = Counter(atom.GetSymbol() for atom in mol_hyd.GetAtoms())
    element_percentages = {
        element: f"{round((Chem.PeriodicTable.GetAtomicWeight(Chem.GetPeriodicTable(), element) * count / mol_wt) * 100, 3)}%"
        for element, count in element_counts.items()
    }

    #Parameters output Info
    params = {
        "Molecular Formula": mol_formula, 
        "Molecular Weight": f"{round(mol_wt,3)} g/mol",
        "Exact Molecular Mass": f"{round(exact_mol_wt,4)} Da",
        "Composition": str(element_percentages)
    }
    return params



