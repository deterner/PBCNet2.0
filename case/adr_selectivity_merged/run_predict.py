"""
Pairwise PBCNet2.0 relative binding-affinity predictions for epinephrine vs.
norepinephrine at alpha1A-AR and beta1-AR, using ONE shared pocket per
subtype (see README.md and merge_pockets.py for how the shared pocket was
built by superimposing the norepinephrine-bound structure onto the
epinephrine-bound one).

Usage (from repo root, with model deps installed):
    python case/adr_selectivity_merged/run_predict.py
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

entries = {
    "alpha1A_epi":    os.path.join(HERE, "alpha1A/epinephrine.pkl"),
    "alpha1A_norepi": os.path.join(HERE, "alpha1A/norepinephrine.pkl"),
    "beta1_epi":      os.path.join(HERE, "beta1/epinephrine.pkl"),
    "beta1_norepi":   os.path.join(HERE, "beta1/norepinephrine.pkl"),
}

names = list(entries.keys())
N1, N2, D1, D2, L = [], [], [], [], []
for n1 in names:
    for n2 in names:
        N1.append(n1)
        D1.append(entries[n1])
        N2.append(n2)
        D2.append(entries[n2])
        L.append(0.0)

csv_path = os.path.join(HERE, "merged_predict.csv")
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
