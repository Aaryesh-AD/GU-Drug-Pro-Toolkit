import joblib
from rdkit import Chem
from rdkit.Chem import AllChem

# Load the trained model
model_path = 'model/optimized_svm_logD_model.pkl'
# model_path = 'random_forest_logD_model.pkl'
model = joblib.load(model_path)

# Function to compute Morgan fingerprints for a given SMILES string
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))

def predictor_logD(smiles):
    descriptors = compute_descriptors(smiles)
    predicted_logD = model.predict([descriptors])[0]
    pred_logD = round(predicted_logD,3) 
    return pred_logD

