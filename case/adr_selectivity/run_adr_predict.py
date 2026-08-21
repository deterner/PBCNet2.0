"""
Pairwise PBCNet2.0 relative binding-affinity predictions for epinephrine vs.
norepinephrine across several adrenergic receptor subtypes, each with its own
natively resolved pocket (see README.md in this folder for the receptor/PDB
list and how the pockets + ligand poses here were produced).

Usage (from repo root, with model deps installed):
    python case/adr_selectivity/run_adr_predict.py
"""
import os
import sys
import pandas as pd
import torch
from torch.utils.data import DataLoader

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))  # .../PBCNet2.0

sys.path.append(os.path.join(REPO_ROOT, "model_code"))
from Dataloader.dataloader import LeadOptDataset, collate_fn
from predict.predict import predict

# name -> (pkl path, receptor subtype, ligand)
entries = {
    "alpha1A_epi":           (os.path.join(HERE, "alpha1A_epi/epinephrine.pkl"),            "alpha1A-AR (8THL)", "epinephrine"),
    "alpha1A_norepi":        (os.path.join(HERE, "alpha1A_norepi/norepinephrine.pkl"),       "alpha1A-AR (7YMH)", "norepinephrine"),
    "alpha2A_norepi":        (os.path.join(HERE, "alpha2A_norepi/norepinephrine.pkl"),       "alpha2A-AR (7EJ0)", "norepinephrine"),
    "beta1_epi":             (os.path.join(HERE, "beta1_epi/epinephrine.pkl"),               "beta1-AR (7BTS)",   "epinephrine"),
    "beta1_norepi":          (os.path.join(HERE, "beta1_norepi/norepinephrine.pkl"),         "beta1-AR (7BU6)",   "norepinephrine"),
    "beta2_epi":             (os.path.join(HERE, "beta2_epi/epinephrine.pkl"),               "beta2-AR (4LDO)",   "epinephrine"),
    "beta3_epi":             (os.path.join(HERE, "beta3_epi/epinephrine.pkl"),               "beta3-AR (9IJE)",   "epinephrine"),
    "alpha2A_epi_docked":    (os.path.join(HERE, "alpha2A_docked_ref/epinephrine.pkl"),      "alpha2A-AR (7EJA, docked)", "epinephrine"),
    "alpha2A_norepi_docked": (os.path.join(HERE, "alpha2A_docked_ref/norepinephrine.pkl"),   "alpha2A-AR (7EJA, docked)", "norepinephrine"),
}

names = list(entries.keys())
N1, N2, D1, D2, L = [], [], [], [], []
for n1 in names:
    for n2 in names:
        N1.append(n1)
        D1.append(entries[n1][0])
        N2.append(n2)
        D2.append(entries[n2][0])
        L.append(0.0)

csv_path = os.path.join(HERE, "adr_predict.csv")
pd.DataFrame({
    "lig1": N1, "lig2": N2,
    "Label": L, "Label1": L, "Label2": L,
    "dir_1": D1, "dir_2": D2,
}).to_csv(csv_path, index=False)

device = torch.device("cpu")
model = torch.load(os.path.join(REPO_ROOT, "PBCNet2.pth"), map_location=device, weights_only=False)
model.to(device)

test_dataset = LeadOptDataset(csv_path)
test_dataloader = DataLoader(test_dataset, collate_fn=collate_fn, batch_size=8, shuffle=False)
_, _, _, _, valid_prediction, _, _, _, _ = predict(model, test_dataloader, device)

df = pd.read_csv(csv_path)
df["pred_delta_pIC50"] = valid_prediction
df.to_csv(csv_path, index=False)

print(df[["lig1", "lig2", "pred_delta_pIC50"]].to_string(index=False))
