from rdkit import Chem
from rdkit.Chem import QED
from props_ptr import sa_score_calc

def calculate_sa_score(smiles):
    mol = Chem.MolFromSmiles(smiles)
    sa_score = sa_score_calc.synth_access_score(mol)
    qed_score = QED.qed(mol)
    score = {
        "Synthetic Accessibility Score": (round(sa_score,3)),
        "QED Score" :(round(qed_score,3))
    }
    return score

