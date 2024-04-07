from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, Crippen
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem

# Load the trained model
model_path = 'models/optimized_svm_logD_model.pkl'
# model_path = 'random_forest_logD_model.pkl'
model = joblib.load(model_path)

# Function to compute Morgan fingerprints for a given SMILES string
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))

def count_carbon_atoms(mol):
    return sum(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())

def compute_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    descriptors = compute_descriptors(smiles)
    predicted_logD = model.predict([descriptors])[0]
    pred_logD = round(predicted_logD,3) 
    return {
        "HDonors": Lipinski.NumHDonors(mol),
        "HAcceptors": Lipinski.NumHAcceptors(mol),
        "MW": round(Descriptors.MolWt(mol), 3),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "AtomCount": mol.GetNumAtoms(),
        "MolarRefractivity": round(Crippen.MolMR(mol), 3),
        "RotatableBonds": Lipinski.NumRotatableBonds(mol),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "RingCount": Lipinski.RingCount(mol),
        "NumCarbon": count_carbon_atoms(mol),
        "NumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "logD": pred_logD
    }

def drug_likeness_properties(smiles):
    properties = compute_properties(smiles)
    
    rules = {
        "Lipinski's Rule": {"HDonors": 5, "HAcceptors": 10, "MW": 500, "LogP": 5},
        "Ghose's Rule": {"LogP": (-0.4, 5.6), "MW": (160, 480), "AtomCount": (20, 70), "MolarRefractivity": (40, 130)},
        "Veber's Rule": {"RotatableBonds": 10, "TPSA": 140},
        "Egan's Rule": {"LogP": 5.88, "TPSA": 131.6},
        "Muegge's Rule": {"MW": (200, 600), "LogP": (-2, 5), "RingCount": (0, 7), "NumCarbon": (4, None), "NumHeteroatoms": (1, None), "RotatableBonds": (0, 15), "HAcceptors": (0, 10), "HDonors": (0, 5)},
        "Pfizer Rule": {"TPSA": 75, "LogP": (3, None)},
        "GSK Rule": {"MW": 400, "LogP": 4},
        "Golden Triangle": {"MW":(200,500),"logD":(-2,5)} 
    }

    violations = {}
    for rule_name, criteria in rules.items():
        rule_violations = []
        for prop, threshold in criteria.items():
            prop_value = properties[prop]
            if isinstance(threshold, tuple):
                lower, upper = threshold
                if (lower is not None and prop_value < lower) or (upper is not None and prop_value > upper):
                    rule_violations.append((prop, prop_value))
            elif prop_value > threshold:
                rule_violations.append((prop, prop_value))
        if rule_violations:
            violations[rule_name] = {"Status": "Fail", "Details": rule_violations}
        else:
            violations[rule_name] = {"Status": "Pass", "Details": "Satisfies all criteria"}

    return violations

def rules_drg(smiles):
    results = {}
    violations = drug_likeness_properties(smiles)

    # Update in the rules_thresholds dictionary
    rules_thresholds = {
        # Update all rules to use "<=" or ">=" instead of "≤" and "≥"
        "Lipinski's Rule": {"HDonors": "<= 5", "HAcceptors": "<= 10", "MW": "<= 500", "LogP": "<= 5"},
        "Ghose's Rule": {"LogP": "(between -0.4 and 5.6)", "MW": "(between 160 and 480)", "AtomCount": "(between 20 and 70)", "MolarRefractivity": "(between 40 and 130)"},
        "Veber's Rule": {"RotatableBonds": "<= 10", "TPSA": "<= 140"},
        "Egan's Rule": {"LogP": "<= 5.88", "TPSA": "<= 131.6"},
        "Muegge's Rule": {
            "MW": "(between 200 and 600)", "LogP": "(between -2 and 5)", "RingCount": "(between 0 and 7)",
            "NumCarbon": ">= 4", "NumHeteroatoms": ">= 1", "RotatableBonds": "(0, 15)",
            "HAcceptors": "(between 0 and 10)", "HDonors": "(between 0 and 5)"
        },
        "Pfizer Rule": {"TPSA": "<= 75", "LogP": ">= 3"},
        "GSK Rule": {"MW": "<= 400", "LogP": "<= 4"},
        "Golden Triangle": {"MW":"between 200 and 500","logD":"between -2 and 5"} 
    }


    for rule, info in violations.items():
        if info["Status"] == "Fail":  # If the rule is violated
            # Prepare the violation details with expected ranges
            violation_details = ", ".join([
                f"{prop}: {value} (Expected: {rules_thresholds[rule].get(prop, 'N/A')})"
                for prop, value in info["Details"]
            ])
            # Include the number of violations in the message
            num_violations = len(info["Details"])
            results[rule] = f"Fail ({num_violations} violations) - {violation_details}"
        else:
            # If the rule is passed, just indicate it satisfies all criteria
            results[rule] = "Pass - Satisfies all criteria"

    return results

