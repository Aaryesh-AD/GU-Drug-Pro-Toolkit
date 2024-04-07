import numpy as np
import pandas as pd
import pickle, os
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from pathlib import Path
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem
from props_ptr.pka import get_pKa


# Load the trained model
model_path = 'models/optimized_svm_logD_model.pkl'
# model_path = 'random_forest_logD_model.pkl'
model = joblib.load(model_path)

# Function to compute Morgan fingerprints for a given SMILES string
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))

# Function to calculate the proportion of aromatic atoms in a molecule
def AromaticProportion(molecule):
    aromatic_atoms = [molecule.GetAtomWithIdx(i).GetIsAromatic() for i in range(molecule.GetNumAtoms())]
    aromatic_atom_count = sum(aromatic_atoms)
    heavy_atom_count = Descriptors.HeavyAtomCount(molecule)
    aromatic_proportion = aromatic_atom_count / heavy_atom_count
    return aromatic_proportion

# Function to generate molecular descriptors from a list of SMILES strings
def generate_descriptors(smiles_list, verbose=False):
    if isinstance(smiles_list, str):
        smiles_list = [smiles_list]  # Ensure it's a list even for a single SMILES string

    molecule_data = []
    for smiles in smiles_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule:
            molecule_data.append(molecule)

    descriptors_list = []
    for molecule in molecule_data:
        descriptors = [
            Descriptors.MolLogP(molecule),
            Descriptors.MolWt(molecule),
            Descriptors.NumRotatableBonds(molecule),
            AromaticProportion(molecule)
        ]
        descriptors_list.append(descriptors)

    column_names = ["MolLogP", "MolWt", "NumRotatableBonds", "AromaticProportion"]
    descriptors_df = pd.DataFrame(descriptors_list, columns=column_names)
    return descriptors_df

# Function to assign solubility class based on log S value
def assign_solubility_class(logS):
    if logS <= -10:
        return "Insoluble"
    elif logS <= -6:
        return "Poorly Soluble"
    elif logS <= -4:
        return "Moderately Soluble"
    elif logS < -2:
        return "Soluble"
    elif logS < 0:
        return "Very Soluble"
    else:
        return "Highly Soluble"

# Load the model outside of your function to improve efficiency
model_path = "models/ESOL_solubility_model.pkl" # Relative path from the script to the model
loaded_model = pickle.load(open(model_path, 'rb'))

def molar_sol(logS, molecular_weight):
    molar_solubility = 10 ** logS  # mol/L
    return molar_solubility
    
def convert_logS_to_mg_per_ml(logS, molecular_weight):
    molar_solubility = 10 ** logS  # mol/L
    solubility_mg_per_ml = molar_solubility * molecular_weight  # mg/mL
    return solubility_mg_per_ml

# Example usage

def predict_sol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    descriptors = generate_descriptors([smiles])
    prediction = loaded_model.predict(descriptors)
    logS = round(float(prediction), 3)
    solubility_class = assign_solubility_class(logS)
    molecular_weight = descriptors['MolWt'].iloc[0]  # assuming the molecular weight is part of the descriptors
    solubility_mg_per_ml = convert_logS_to_mg_per_ml(logS, molecular_weight)
    molar_solubility = molar_sol(logS, molecular_weight)
    descriptors = compute_descriptors(smiles)
    predicted_logD = model.predict([descriptors])[0]
    params = {"LogS":logS,
          "Molar Solubility (mol/l)": round(molar_solubility,5),
          "Solubility (mg/ml)": round(solubility_mg_per_ml,5),
          "Solubility Class":solubility_class,
          "WLogP" : round(Crippen.MolLogP(mol),3),
          "LogD": round(predicted_logD,3),
          "Nature": get_pKa(smiles)
          }
    return params
