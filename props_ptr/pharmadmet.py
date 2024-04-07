import joblib
from rdkit import Chem
from rdkit.Chem import AllChem

# List of model names
model_names = [
    "AMES", "BBB_Martins", "Bioavailability_Ma", "CYP1A2_Veith", "CYP2C19_Veith",
    "CYP2C9_Substrate_CarbonMangels", "CYP2C9_Veith", "CYP2D6_Substrate_CarbonMangels",
    "CYP2D6_Veith", "CYP3A4_Substrate_CarbonMangels", "CYP3A4_Veith", "Carcinogens_Lagunin",
    "ClinTox", "DILI", "HIA_Hou", "NR-AR-LBD", "NR-AR", "NR-AhR", "NR-Aromatase", "NR-ER-LBD",
    "NR-ER", "NR-PPAR-gamma", "PAMPA_NCATS", "Pgp_Broccatelli", "SR-ARE", "SR-ATAD5", "SR-HSE",
    "SR-MMP", "SR-p53", "Skin_Reaction", "hERG", "Caco2_Wang", "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ", "Half_Life_Obach", "HydrationFreeEnergy_FreeSolv", "LD50_Zhu",
    "Lipophilicity_AstraZeneca", "PPBR_AZ", "Solubility_AqSolDB", "VDss_Lombardo"
]

# Function to compute Morgan fingerprints for a given SMILES string
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))

# Function to predict with a given model
def predict_with_model(model_path, descriptors):
    model = joblib.load(model_path)
    predicted_value = model.predict([descriptors])[0]
    return round(predicted_value, 3)

# Main prediction function
def predict_for_all_models(smiles):
    descriptors = compute_descriptors(smiles)
    predictions = {}
    model_names = [
    "AMES", "BBB_Martins", "Bioavailability_Ma", "CYP1A2_Veith", "CYP2C19_Veith",
    "CYP2C9_Substrate_CarbonMangels", "CYP2C9_Veith", "CYP2D6_Substrate_CarbonMangels",
    "CYP2D6_Veith", "CYP3A4_Substrate_CarbonMangels", "CYP3A4_Veith", "Carcinogens_Lagunin",
    "ClinTox", "DILI", "HIA_Hou", "NR-AR-LBD", "NR-AR", "NR-AhR", "NR-Aromatase", "NR-ER-LBD",
    "NR-ER", "NR-PPAR-gamma", "PAMPA_NCATS", "Pgp_Broccatelli", "SR-ARE", "SR-ATAD5", "SR-HSE",
    "SR-MMP", "SR-p53", "Skin_Reaction", "hERG", "Caco2_Wang", "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ", "Half_Life_Obach", "HydrationFreeEnergy_FreeSolv", "LD50_Zhu",
    "Lipophilicity_AstraZeneca", "PPBR_AZ", "Solubility_AqSolDB", "VDss_Lombardo"
]

    property_display_names = [
        "Ames: Ames Mutagenicity",
        "BBB_Martins: Blood-Brain Barrier Penetration (Martins)",
        "Bioavailability_Ma: Bioavailability (Ma)",
        "CYP1A2_Veith: CYP1A2 Inhibition (Veith)",
        "CYP2C19_Veith: CYP2C19 Inhibition (Veith)",
        "CYP2C9_Substrate_CarbonMangels: CYP2C9 Substrate (CarbonMangels)",
        "CYP2C9_Veith: CYP2C9 Inhibition (Veith)",
        "CYP2D6_Substrate_CarbonMangels: CYP2D6 Substrate (CarbonMangels)",
        "CYP2D6_Veith: CYP2D6 Inhibition (Veith)",
        "CYP3A4_Substrate_CarbonMangels: CYP3A4 Substrate (CarbonMangels)",
        "CYP3A4_Veith: CYP3A4 Inhibition (Veith)",
        "Carcinogens_Lagunin: Carcinogenicity (Lagunin)",
        "ClinTox: Clinical Toxicity",
        "DILI: Drug-Induced Liver Injury",
        "HIA_Hou: Human Intestinal Absorption (Hou)",
        "NR-AR-LBD: Nuclear Receptor AR Ligand Binding Domain",
        "NR-AR: Nuclear Receptor AR",
        "NR-AhR: Nuclear Receptor AhR",
        "NR-Aromatase: Nuclear Receptor Aromatase",
        "NR-ER-LBD: Nuclear Receptor ER Ligand Binding Domain",
        "NR-ER: Nuclear Receptor ER",
        "NR-PPAR-gamma: Nuclear Receptor PPAR-gamma",
        "PAMPA_NCATS: Parallel Artificial Membrane Permeability Assay (NCATS)",
        "Pgp_Broccatelli: P-glycoprotein Substrate (Broccatelli)",
        "SR-ARE: Stress Response - Antioxidant Response Element",
        "SR-ATAD5: Stress Response - ATAD5",
        "SR-HSE: Stress Response - Heat Shock Response Element",
        "SR-MMP: Stress Response - Mitochondrial Membrane Permeability",
        "SR-p53: Stress Response - p53 Response Element",
        "Skin_Reaction: Skin Reaction",
        "hERG: hERG Inhibition",
        "Caco2_Wang: Caco-2 Permeability (Wang)",
        "Clearance_Hepatocyte_AZ: Hepatocyte Clearance (AstraZeneca)",
        "Clearance_Microsome_AZ: Microsomal Clearance (AstraZeneca)",
        "Half_Life_Obach: Half-Life (Obach)",
        "HydrationFreeEnergy_FreeSolv: Hydration Free Energy (FreeSolv)",
        "LD50_Zhu: LD50 (Zhu)",
        "Lipophilicity_AstraZeneca: Lipophilicity (AstraZeneca)",
        "PPBR_AZ: Plasma Protein Binding (AstraZeneca)",
        "Solubility_AqSolDB: Solubility (AqSolDB)",
        "VDss_Lombardo: Volume of Distribution (Lombardo)"
    ]

    for idx, model_name in enumerate(model_names):
        model_path = f'models/{model_name}_model.pkl'  # Adjust path as needed
        predicted_value = predict_with_model(model_path, descriptors)
        predictions[property_display_names[idx]] = predicted_value
    return predictions

