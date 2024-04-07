# Import necessary libraries and modules
import math
import os
import pickle
import gzip
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Global variable to store fragment scores
fragment_scores = None

# Function to load fragment scores from a compressed file
def load_frag_scores(file_name='fpscores'):
    global fragment_scores
    full_path = os.path.join(os.path.dirname(__file__), file_name) if file_name == "fpscores" else file_name
    with gzip.open(f'{full_path}.pkl.gz', 'rb') as file:
        data = pickle.load(file)
    score_dict = {item[j]: float(item[0]) for item in data for j in range(1, len(item))}
    fragment_scores = score_dict

# Function to calculate the number of bridgehead and spiro atoms
def calc_bridge_spiro_atoms(molecule, ring_info=None):
    spiro_atoms = rdMolDescriptors.CalcNumSpiroAtoms(molecule)
    bridgehead_atoms = rdMolDescriptors.CalcNumBridgeheadAtoms(molecule)
    return bridgehead_atoms, spiro_atoms

# Function to calculate the synthetic accessibility score
# Function to calculate the synthetic accessibility score
def synth_access_score(mol):
    global fragment_scores
    if fragment_scores is None:
        load_frag_scores()

    fingerprint = rdMolDescriptors.GetMorganFingerprint(mol, 2)  # radius of 2
    nonzero_elements = fingerprint.GetNonzeroElements()
    frag_score = sum(fragment_scores.get(bitId, -4) * count for bitId, count in nonzero_elements.items()) / sum(nonzero_elements.values())

    atom_count = mol.GetNumAtoms()
    chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    ring_info = mol.GetRingInfo()
    bridgeheads, spiro = calc_bridge_spiro_atoms(mol, ring_info)
    macrocycles = sum(1 for x in ring_info.AtomRings() if len(x) > 8)

    penalties = sum([atom_count**1.005 - atom_count] + [math.log10(x + 1) for x in [chiral_centers, spiro, bridgeheads]])
    if macrocycles > 0: penalties += math.log10(2)

    if atom_count > len(nonzero_elements):
        penalties += 0.5 * math.log(float(atom_count) / len(nonzero_elements))

    sa_score = 11. - (frag_score + penalties + 3.5) / 6.5 * 9.

    # Adjust sa_score smoothly for values above 8
    if sa_score > 8:
        sa_score = 8. + math.log(sa_score - 7.) if sa_score > 8 else sa_score
        sa_score = min(sa_score, 10.0)  # Ensure sa_score does not exceed 10
    elif sa_score < 1.:
        sa_score = 1.0  # Ensure sa_score does not go below 1

    return sa_score



# Main function to process molecules and print their synthetic accessibility score
def process_molecules(molecules):
    print('SMILES\tName\tSA_Score')
    for idx, molecule in enumerate(molecules):
        if molecule:
            score = synth_access_score(molecule)
            smiles = Chem.MolToSmiles(molecule)
            print(f'{smiles}\t{molecule.GetProp("_Name")}\t{score:.3f}')

if __name__ == '__main__':
    import sys
    import time

    start_time = time.time()
    load_frag_scores("fpscores")
    load_time = time.time()

    molecule_supplier = Chem.SmilesMolSupplier(sys.argv[1])
    process_time = time.time()
    process_molecules(molecule_supplier)
    end_time = time.time()

    print(f'Reading took {load_time - start_time:.2f} seconds. Processing took {end_time - process_time:.2f} seconds', file=sys.stderr)
