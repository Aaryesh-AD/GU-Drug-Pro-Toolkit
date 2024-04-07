import os, joblib
from props_ptr.data_utils import DataUtils
from rdkit import Chem


class PkaPredictor(object):

    def __init__(self, model_dir="models", feature_type="morgan+macc"):
        self._load_model(clf_modelpath=os.path.join(model_dir, "pka_classification.pkl"),
                         acidic_modelpath=os.path.join(model_dir, "pka_acidic_regression.pkl"),
                         basic_modelpath=os.path.join(model_dir, "pka_basic_regression.pkl"))
        self.feature_type = feature_type

    def _load_model(self, clf_modelpath, acidic_modelpath, basic_modelpath):
        self.acidic_reg = joblib.load(acidic_modelpath)
        self.basic_reg = joblib.load(basic_modelpath)
        self.clf = joblib.load(clf_modelpath)

    def predict(self, mols):
        mols_features = DataUtils.get_molecular_features(mols, self.feature_type)
        clf_labels = self.clf.predict(mols_features)
        acidic_scores = self.acidic_reg.predict(mols_features)
        basic_scores = self.basic_reg.predict(mols_features)
        rets = []
        for idx, clf_label in enumerate(clf_labels):
            results = []
            if clf_label[0] == 1: results.append(f"Acidic with pKa: {acidic_scores[idx]:.2f}")
            if clf_label[1] == 1: results.append(f"Basic with pKa: {basic_scores[idx]:.2f}")
            rets.append(", ".join(results))
        return rets

def get_pKa(smiles):
    mol = Chem.MolFromSmiles(smiles)
    pka_predictor = PkaPredictor()
    return pka_predictor.predict([mol])[0]
