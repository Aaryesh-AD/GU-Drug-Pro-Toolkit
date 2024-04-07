import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen
import sys
from props_ptr.ESOL_prop import predict_sol

def graph_prop_1(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    molecule_h = Chem.AddHs(molecule)
    fsp3 = round(rdMolDescriptors.CalcFractionCSP3(molecule_h),3)
    tpsa = round(rdMolDescriptors.CalcTPSA(molecule_h),3)
    rot= Lipinski.NumRotatableBonds(molecule_h)
    mol_wt = round(Descriptors.MolWt(molecule), 3)
    logp = round(Crippen.MolLogP(molecule),3)
    sol = predict_sol(smiles)
    logs = sol['LogS']
    params = np.array([logp, mol_wt, tpsa, logs, fsp3, rot])
    return params
    


def plot_radar_chart(root, compound_values, figsize=(4, 4), font_size=8):
    labels = np.array(['            Lipophilicity (logP)', 'Molecular Weight', 'Polar Surface Area (TPSA)', 'Solubility (logS)      ', 'Saturation (FSP3)', 'Flexibility (Rotatable Bonds)'])
    upper_limits = np.array([5.0, 500, 130, 0, 1, 9])
    lower_limits = np.array([-0.7, 150, 20, -6, 0.25, 0])
    normalized_values = (compound_values - lower_limits) / (upper_limits - lower_limits) * 1
    normalized_upper = np.ones(len(upper_limits))

    normalized_values = np.append(normalized_values, normalized_values[0])
    normalized_upper = np.append(normalized_upper, normalized_upper[0])

    num_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]

    # Clear the existing figure if it exists and create a new one
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, polar=True)

    ax.plot(angles, normalized_upper, color='green', linewidth=1, linestyle='dotted', label='Optimal Range')
    ax.fill(angles, normalized_upper, color='green', alpha=0.25)
    ax.plot(angles, normalized_values, color='red', linewidth=1, linestyle='solid', label='Compound Properties')
    ax.fill(angles, normalized_values, color='red', alpha=0.25)

    ax.set_xticks(angles[:-1])
    ax.set_ylim(0, 2)

    for label, angle in zip(labels, angles[:-1]):
        ax.text(angle, 2.7, label, horizontalalignment='center', color='black', fontsize=font_size)

    ax.spines['polar'].set_color((0, 0, 0, 0.5))
    ax.legend(loc='upper right', bbox_to_anchor=(-0.3, 0.1), fontsize=font_size)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    # Attach the figure to the Tkinter root. The canvas is created here.
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack(side='top', fill='both', expand=True)
    canvas.draw()

    # Return both the canvas and fig to allow external manipulation
    return canvas, fig


def plot_radar_chart_web(compound_values, figsize=(4, 4), font_size=8,filename='temp_img/radar_chart.png'):
    labels = np.array(['            Lipophilicity (logP)', 'Molecular Weight', 'Polar Surface Area (TPSA)', 'Solubility (logS)      ', 'Saturation (FSP3)', 'Flexibility (Rotatable Bonds)'])
    upper_limits = np.array([5.0, 500, 130, 0, 1, 9])
    lower_limits = np.array([-0.7, 150, 20, -6, 0.25, 0])
    normalized_values = (compound_values - lower_limits) / (upper_limits - lower_limits) * 1
    normalized_upper = np.ones(len(upper_limits))

    normalized_values = np.append(normalized_values, normalized_values[0])
    normalized_upper = np.append(normalized_upper, normalized_upper[0])

    num_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]

    # Clear the existing figure if it exists and create a new one
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, polar=True)

    ax.plot(angles, normalized_upper, color='green', linewidth=1, linestyle='dotted', label='Optimal Range')
    ax.fill(angles, normalized_upper, color='green', alpha=0.25)
    ax.plot(angles, normalized_values, color='red', linewidth=1, linestyle='solid', label='Compound Properties')
    ax.fill(angles, normalized_values, color='red', alpha=0.25)

    ax.set_xticks(angles[:-1])
    ax.set_ylim(0, 2)

    for label, angle in zip(labels, angles[:-1]):
        ax.text(angle, 2.7, label, horizontalalignment='center', color='black', fontsize=font_size)

    ax.spines['polar'].set_color((0, 0, 0, 0.5))
    ax.legend(loc='upper right', bbox_to_anchor=(-0.3, 0.1), fontsize=font_size)
    plt.title('Bioavailability Spiderweb', pad=40)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    plt.savefig(filename, bbox_inches='tight')  # Save the figure as an image
    plt.close()
