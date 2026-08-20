import sys
import os
import torch
import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
from torch.utils.data import DataLoader

# Add local repository paths
base_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(base_dir)
sys.path.append(base_dir)
sys.path.append(repo_root)
sys.path.append(os.path.join(repo_root, "model_code"))

from Graph2pickle import graph_save
from model_code.Dataloader.dataloader import LeadOptDataset, collate_fn
from model_code.predict.predict import predict

# Set target data directory
data_dir = os.path.join(repo_root, "case", "toy_data")

# 1. Pocket Extraction
def extract(ligand, pdb):
    parser = PDBParser(QUIET=True)
    if not os.path.exists(pdb):
        return
    structure = parser.get_structure("protein", pdb)
    lp = [l.GetConformer().GetPositions() for l in ligand]
    ligand_positions = np.concatenate(lp)

    class ResidueSelect(Select):
        def accept_residue(self, residue):
            residue_positions = np.array([
                np.array(list(atom.get_vector()))
                for atom in residue.get_atoms() if "H" not in atom.get_id()
            ])
            if len(residue_positions.shape) < 2:
                return 0
            min_dis = np.min(distance_matrix(residue_positions, ligand_positions))
            return 1 if min_dis < 8.0 else 0

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb.replace('protein.pdb', 'pocket.pdb'), ResidueSelect())

def pocket_extract(sdf_files, protein_file):
    ligands = [Chem.MolFromMolFile(a) for a in sdf_files if Chem.MolFromMolFile(a) is not None]
    extract(ligands, protein_file)

sdfs = [os.path.join(data_dir, i) for i in os.listdir(data_dir) if i.endswith('.sdf')]
protein_file = os.path.join(data_dir, 'protein.pdb')

for sdf in sdfs:
    pocket_extract([sdf], protein_file)

# 2. Graph Pickles Generation
pocket_pdb = os.path.join(data_dir, 'pocket.pdb')
for sdf in sdfs:
    graph_save(sdf, pocket_pdb, sdf.replace('.sdf', '.pkl'))

# 3. Pairwise File Generation
pickles = [os.path.join(data_dir, i) for i in os.listdir(data_dir) if i.endswith('.pkl')]
ic50 = [7.47, 8.52, 8.3, 7.77, 7.6, 8.1]  # Dummy or experimental pIC50 values

N1, N2, D1, D2, L = [], [], [], [], []
for i, pkl1 in enumerate(pickles):
    for j, pkl2 in enumerate(pickles):
        N1.append(pkl1.split(os.sep)[-1])
        D1.append(pkl1)
        N2.append(pkl2.split(os.sep)[-1])
        D2.append(pkl2)
        L.append(ic50[i] - ic50[j])

predict_csv_path = os.path.join(data_dir, 'predict.csv')
pd.DataFrame({
    'lig1': N1, 'lig2': N2, 'Label': L, 
    'Label1': L, 'Label2': L, 'dir_1': D1, 'dir_2': D2
}).to_csv(predict_csv_path, index=False)

# 4. Model Inference
# Set to torch.device('cuda') if GPU is available
device = torch.device('cpu') 

model_path = os.path.join(os.path.dirname(base_dir), "PBCNet2.pth")
model = torch.load(model_path, map_location=device, weights_only=False)
model.to(device)

test_dataset = LeadOptDataset(predict_csv_path)
test_dataloader = DataLoader(test_dataset, collate_fn=collate_fn, batch_size=8, shuffle=False)

_, _, _, _, valid_prediction, _, _, _, _ = predict(model, test_dataloader, device)

# Save predictions back to CSV
df_file = pd.read_csv(predict_csv_path)
df_file['pre'] = valid_prediction
df_file.to_csv(predict_csv_path, index=False)

print("Inference completed successfully!")
print(df_file.head())