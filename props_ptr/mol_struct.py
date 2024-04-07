from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def calculate_molecular_parameters(smiles):
    # Convert the SMILES string to a molecule object
    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)

    # Calculating various molecular parameters
    params = {
        "Atom Count": molecule_h.GetNumAtoms(),
        "Heavy Atom Count": Lipinski.HeavyAtomCount(molecule_h),
        "Asymmetric Atom Count": len(Chem.FindMolChiralCenters(molecule_h, includeUnassigned=True)),
        "Rotatable Bond Count": Lipinski.NumRotatableBonds(molecule_h),
        "Ring Count": Lipinski.RingCount(molecule_h),
        "Aromatic Ring Count": Lipinski.NumAromaticRings(molecule_h),
        "Hetero Ring Count": sum(1 for ring in molecule_h.GetRingInfo().AtomRings() if any(molecule_h.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring)),
        "FSP3": round(rdMolDescriptors.CalcFractionCSP3(molecule_h),3),
        "Hydrogen Bond Donor Count": Lipinski.NumHDonors(molecule_h),
        "Hydrogen Bond Acceptor Count": Lipinski.NumHAcceptors(molecule_h),
        "Formal Charge": Chem.GetFormalCharge(molecule_h),
        "Topological Polar Surface Area (Ångströms^2)": round(rdMolDescriptors.CalcTPSA(molecule_h),3),
        }
    return params

