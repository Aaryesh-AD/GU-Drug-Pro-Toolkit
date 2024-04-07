from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, MACCSkeys
import pandas as pd
import os

def compute_descriptors(smiles):
    # Convert SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")

    # Dictionary to hold our descriptors
    descriptors = {}

    # 2D Descriptors
    descriptors_2D = {desc_name: desc_func(mol) for desc_name, desc_func in Descriptors.descList}
    descriptors['2D'] = descriptors_2D

    # Morgan Fingerprints (also known as ECFP)
    for radius in [1, 2, 3]:  # Corresponds to ECFP2, ECFP4, ECFP6
        morgan_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=2048)
        descriptors[f'ECFP{2*radius}'] = {f'Bit_{i}': morgan_fp.GetBit(i) for i in range(2048)}

    # MACCS Keys
    maccs_keys = MACCSkeys.GenMACCSKeys(mol)
    descriptors['MACCS'] = {f'Bit_{i}': maccs_keys.GetBit(i) for i in range(167)}

    # MQN Descriptors
    mqn = rdMolDescriptors.MQNs_(mol)
    descriptors['MQN'] = {f'MQN_{i+1}': val for i, val in enumerate(mqn)}

    return descriptors


def export_xls(smiles):
    descriptors = compute_descriptors(smiles)
    # Export 2D Descriptors
    df_2D = pd.DataFrame([descriptors['2D']])
    df_2D.to_excel('Export_2D.xlsx', index=False)
    # Export ECFP Descriptors
    for radius in [1, 2, 3]:  # Corresponds to ECFP2, ECFP4, ECFP6
        df_ECFP = pd.DataFrame([descriptors[f'ECFP{2*radius}']])
        df_ECFP.to_excel(f'Export_ECFP{2*radius}.xlsx', index=False)
    # Export MACCS Keys
    df_MACCS = pd.DataFrame([descriptors['MACCS']])
    df_MACCS.to_excel('Export_MACCS.xlsx', index=False)
    # Export MQN Descriptors
    df_MQN = pd.DataFrame([descriptors['MQN']])
    df_MQN.to_excel('Export_MQN.xlsx', index=False)


def export_csv(smiles):
    descriptors = compute_descriptors(smiles)   
    #Export_2D
    df_2D = pd.DataFrame([descriptors['2D']])
    df_2D.to_csv('Export_2D.csv', index=False)
    for radius in [1, 2, 3]:  # Corresponds to ECFP2, ECFP4, ECFP6
        df_ECFP = pd.DataFrame([descriptors[f'ECFP{2*radius}']])
        df_ECFP.to_csv(f'Export_ECFP{2*radius}.csv', index=False)
    # Export MACCS Keys
    df_MACCS = pd.DataFrame([descriptors['MACCS']])
    df_MACCS.to_csv('Export_MACCS.csv', index=False)
    # Export MQN Descriptors
    df_MQN = pd.DataFrame([descriptors['MQN']])
    df_MQN.to_csv('Export_MQN.csv', index=False)


def compute_and_export_descriptors(smiles, export_format, options, export_directory):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")

    descriptors = {}

    if options['2D']:
        descriptors_2D = {desc_name: desc_func(mol) for desc_name, desc_func in Descriptors.descList}
        descriptors['2D'] = descriptors_2D

    if options['ECFP']:
        for radius in [1, 2, 3]:
            morgan_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=2048)
            descriptors[f'ECFP{2*radius}'] = {f'Bit_{i}': morgan_fp.GetBit(i) for i in range(2048)}

    if options['MACCS']:
        maccs_keys = MACCSkeys.GenMACCSKeys(mol)
        descriptors['MACCS'] = {f'Bit_{i}': maccs_keys.GetBit(i) for i in range(167)}

    if options['MQN']:
        mqn = rdMolDescriptors.MQNs_(mol)
        descriptors['MQN'] = {f'MQN_{i+1}': val for i, val in enumerate(mqn)}


    for desc_type, desc_data in descriptors.items():
        df = pd.DataFrame([desc_data])
        file_name = f'{smiles}_{desc_type}.{export_format.lower()}'
        full_path = os.path.join(export_directory, file_name)
        
        if export_format == 'csv':
            df.to_csv(full_path, index=False)
        elif export_format == 'xlsx':
            df.to_excel(full_path, index=False,  engine='openpyxl')

