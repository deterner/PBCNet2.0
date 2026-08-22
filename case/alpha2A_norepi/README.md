# PBCNet2.0 relative-affinity prediction: alpha2A-adrenergic receptor (7EJ0)

Self-contained example: prepared structures, DGL graph pickles, and prediction output for a
small analysis run with the pretrained `PBCNet2.pth` checkpoint on the alpha2A-adrenergic
receptor structure PDB [7EJ0](https://www.rcsb.org/structure/7EJ0) (alpha2A-AR + norepinephrine,
GoA-coupled complex). Everything needed to inspect or rerun this analysis lives in this folder;
paths in `predict.csv` are relative to it.

## 1. Objective

PBCNet2.0 predicts **relative** binding affinity (Δ pIC50) between two ligand-pocket complexes —
it does not predict an absolute affinity for a single pose. The question asked here: for the
endogenous catecholamines epinephrine and norepinephrine, what delta does the model predict, and
does using the *receptor-specific* alpha2A pocket (vs. a generic adrenergic pocket used elsewhere
in this repo, under `main/`) change the prediction?

## 2. Methodological framework

PBCNet2.0 (`model_code/models/readout.py::PBCNetv2`) is a Cartesian-tensor equivariant GNN
(TensorNet encoder, `model_code/models/tensornet.py`) that:

1. Encodes each ligand+pocket complex as a DGL heterograph (ligand atoms + nearby protein pocket
   atoms as one node type `atom`, bonds/proximity as edges) and embeds it independently.
2. Subtracts the two complex embeddings (`emb1 - emb2`) and passes the difference through a small
   feed-forward head to output a single scalar: the predicted Δ pIC50 between complex 1 and
   complex 2.
3. Is trained to be anti-symmetric (`g(1,2) ≈ -g(2,1)`); in practice, as seen below, this holds
   only approximately for an off-the-shelf checkpoint.

Three complexes were built and compared pairwise (3×3 = 9 predictions):

| label | files in this folder | ligand source | pocket |
|---|---|---|---|
| `norepinephrine_alpha2A` | `norepinephrine_alpha2A.sdf` / `.pkl` | co-crystallized pose, extracted directly from `7EJ0.pdb` (real, not docked) | `pocket.pdb` — alpha2A-specific, extracted from 7EJ0 |
| `norepinephrine_generic` | `norepinephrine_generic.sdf` / `.pkl` | pre-existing docked pose (copied from `main/norepinephrine.sdf`) | `ne_pocket.pdb` — generic adrenergic pocket |
| `epinephrine_generic` | `epinephrine_generic.sdf` / `.pkl` | pre-existing docked pose (copied from `main/epinephrine_docked_v21.sdf`) | `ne_pocket.pdb` — generic adrenergic pocket |

There is no epinephrine-bound alpha2A structure available, so no receptor-specific epinephrine
pose exists — this is the main asymmetry in the comparison set (see Limitations).

## 3. Steps performed

1. **Pocket extraction** (`main/adr/extract_pockets.py`, run against `7EJ0.pdb`): for residue
   `E5E` (norepinephrine), selected all standard amino-acid residues with any non-H atom within
   8.0 Å of any ligand non-H atom, saved as `pocket.pdb` (BioPython `PDBParser`/`PDBIO`).
2. **Native ligand extraction**: pulled the 12 `HETATM` lines for residue `E5E` out of `7EJ0.pdb`
   directly (real crystal coordinates, no re-docking), parsed with
   `Chem.MolFromPDBBlock(sanitize=False)`, then assigned bond orders and formal charges by
   matching against a SMILES template (`NCC(O)c1ccc(O)c(O)c1`) with
   `AllChem.AssignBondOrdersFromTemplate`, added explicit hydrogens with `Chem.AddHs(addCoords=True)`,
   and wrote the result to `norepinephrine_alpha2A.sdf`.
3. **Graph construction** (`Graph2pickle.graph_save`): converted each ligand SDF + pocket PDB pair
   into a DGL heterograph pickle (`*.pkl` in this folder). Note: DGL graph pickles are not
   guaranteed portable across DGL versions — if a `.pkl` here fails to unpickle in your
   environment, regenerate it from the corresponding `.sdf` + pocket PDB with `graph_save`.
4. **Pairwise CSV** (`predict.csv`): all 9 ordered pairs among the three complexes above, with
   dummy `Label`/`Label1`/`Label2` columns (unused at inference — they only matter for
   training/evaluation against known experimental values), and `dir_1`/`dir_2` pointing at the
   `.pkl` files in this folder.
5. **Inference**: loaded `PBCNet2.pth` (repo root) on CPU (`torch.load(..., map_location='cpu')`),
   ran `model_code/predict/predict.predict(model, loader, device)` over the CSV via
   `LeadOptDataset` + `DataLoader(collate_fn=collate_fn)`, wrote predictions back to `predict.csv`
   as `pred_delta_pIC50`.

## 4. Results

`pred_delta_pIC50` = predicted (pIC50 of `lig1`) − (pIC50 of `lig2`); higher = `lig1` binds
tighter.

| lig1 | lig2 | pred Δ pIC50 |
|---|---|---:|
| norepinephrine_alpha2A | norepinephrine_alpha2A | −0.080 |
| norepinephrine_alpha2A | norepinephrine_generic | +0.531 |
| norepinephrine_alpha2A | epinephrine_generic | −0.908 |
| norepinephrine_generic | norepinephrine_alpha2A | −0.756 |
| norepinephrine_generic | norepinephrine_generic | −0.080 |
| norepinephrine_generic | epinephrine_generic | −1.052 |
| epinephrine_generic | norepinephrine_alpha2A | +0.895 |
| epinephrine_generic | norepinephrine_generic | +1.060 |
| epinephrine_generic | epinephrine_generic | −0.080 |

Full output: [`predict.csv`](./predict.csv).

## 5. Interpretation

- **Epinephrine vs. norepinephrine**: the model consistently predicts epinephrine binds tighter
  than norepinephrine, both in the generic pocket (Δ ≈ +1.06 / −1.05) and when norepinephrine is
  represented by its native alpha2A pose (Δ ≈ +0.89 to −0.91). Direction is stable across all
  four cross-ligand comparisons.
- **Pocket sensitivity**: the same norepinephrine ligand scores differently depending on which
  pocket structure it's paired with (alpha2A-native vs. generic: Δ ≈ +0.53 / −0.76, not
  symmetric — see below). This indicates the model's output is sensitive to the specific pocket
  geometry, as intended by design, though the two pocket structures used here also differ in
  ligand pose source (crystal vs. docked), so this delta conflates pocket identity with pose
  quality.
- **Self-comparison ("diagonal") is not exactly zero**: comparing a complex against itself gives
  −0.080 for all three, rather than the mathematically expected 0. This is the FNN head's fixed
  bias term surviving the `emb1 - emb2 = 0` subtraction, not a property of the ligand or pocket —
  treat it as a constant offset baked into this checkpoint, not signal.
- **Approximate, not exact, anti-symmetry**: `g(A,B)` and `-g(B,A)` are close but not identical
  (e.g. +0.531 vs. −0.756, or +0.895 vs. −0.908). The gap (~0.1–0.2 pIC50) is a rough indicator of
  the model's internal noise floor on this input distribution.

## 6. Limitations

- **No receptor-specific epinephrine pose.** No epinephrine-bound alpha2A crystal structure was
  used; epinephrine is only represented via a generic-pocket docked pose. The alpha2A-specific
  results here are therefore only receptor-accurate for the norepinephrine side of each
  comparison.
- **"Generic" pocket provenance is not receptor-subtype-specific.** `ne_pocket.pdb` was built for
  an earlier, unrelated demo (see `main/run_predict.py`) and is not itself an alpha2A structure —
  comparisons involving it mix receptor identity with pocket/pose source and should not be read
  as "alpha2A vs. alpha2A."
- **Cross-pocket comparison is outside the model's intended use case.** PBCNet2.0 is designed and
  benchmarked for relative binding affinity *within a single target's pocket* (e.g. lead
  optimization series, FEP benchmarks). Comparing complexes built from two different protein
  structures (as done here for `*_alpha2A` vs. `*_generic`) is not part of its validated use
  case; treat those specific deltas as exploratory, not benchmark-grade.
- **No experimental ground truth.** No pIC50/Ki values were supplied for norepinephrine or
  epinephrine at alpha2A — dummy zero labels were used since they only affect ranking-metric
  computation, not the model's Δ prediction. Predicted deltas are not validated against any
  measured affinity in this analysis.
- **Ligand protonation/bond-order assumption.** The alpha2A-native norepinephrine pose was
  reconstructed from heavy-atom crystal coordinates using a SMILES template match
  (`AssignBondOrdersFromTemplate`); the crystal's `FORMUL` record indicates a protonated amine
  (`C8 H12 N O3 1+`), which the template preserves, but this is inferred, not directly observed
  from the electron density.
- **Arbitrary pocket cutoff.** The 8.0 Å residue-inclusion cutoff (`extract_pockets.py`) is a
  fixed heuristic, not tuned per-target; a different cutoff could change predictions.
- **Single checkpoint, no fine-tuning, no uncertainty estimate.** This used the pretrained
  `PBCNet2.pth` as-is (no fine-tuning on adrenergic-receptor or catecholamine data, see
  `model_code/Finetune.py`/`run_fintune.sh` for that workflow), and PBCNet2.0 does not output a
  confidence interval — only the point-estimate deltas above.

## 7. Replication

Requires the `pbcnet` conda environment (see repo-root `README.md` / `CLAUDE.md` for install
steps) with `PBCNet2.pth` present at the repo root, on the `alpha2a_cpu` branch.

```bash
# from the repo root
git checkout alpha2a_cpu

# 1. (re)extract the alpha2A pocket from the raw structure, if pocket.pdb needs regenerating
python main/adr/extract_pockets.py   # writes main/adr/alpha2A_norepi/pocket.pdb from 7EJ0.pdb
```

```python
# 2. extract the native norepinephrine ligand from 7EJ0.pdb (skip if norepinephrine_alpha2A.sdf
#    already exists in this folder)
from rdkit import Chem
from rdkit.Chem import AllChem

lines = [l for l in open("case/alpha2A_norepi/7EJ0.pdb") if l.startswith("HETATM") and "E5E" in l]
mol = Chem.MolFromPDBBlock("".join(lines) + "END\n", sanitize=False, removeHs=False)
template = Chem.MolFromSmiles("NCC(O)c1ccc(O)c(O)c1")
mol2 = Chem.AddHs(AllChem.AssignBondOrdersFromTemplate(template, mol), addCoords=True)
Chem.SDWriter("case/alpha2A_norepi/norepinephrine_alpha2A.sdf").write(mol2)

# 3. build (or rebuild) the DGL graph pickles
import sys; sys.path.append(".")
from Graph2pickle import graph_save
base = "case/alpha2A_norepi/"
graph_save(base+"norepinephrine_alpha2A.sdf", base+"pocket.pdb",     base+"norepinephrine_alpha2A.pkl")
graph_save(base+"norepinephrine_generic.sdf", base+"ne_pocket.pdb",  base+"norepinephrine_generic.pkl")
graph_save(base+"epinephrine_generic.sdf",    base+"ne_pocket.pdb",  base+"epinephrine_generic.pkl")

# 4. build the pairwise CSV and run inference
import os, torch, pandas as pd
from torch.utils.data import DataLoader
sys.path.append("model_code")
from model_code.Dataloader.dataloader import LeadOptDataset, collate_fn
from model_code.predict.predict import predict

pkl_paths = {
    "norepinephrine_alpha2A": base+"norepinephrine_alpha2A.pkl",
    "norepinephrine_generic": base+"norepinephrine_generic.pkl",
    "epinephrine_generic":    base+"epinephrine_generic.pkl",
}
names = list(pkl_paths)
rows = [(n1, n2, pkl_paths[n1], pkl_paths[n2]) for n1 in names for n2 in names]
df = pd.DataFrame(rows, columns=["lig1", "lig2", "dir_1", "dir_2"])
df["Label"] = df["Label1"] = df["Label2"] = 0.0
df.to_csv(base+"predict.csv", index=False)

device = torch.device("cpu")
model = torch.load("PBCNet2.pth", map_location=device, weights_only=False).to(device)
loader = DataLoader(LeadOptDataset(base+"predict.csv"), collate_fn=collate_fn, batch_size=8, shuffle=False)
_, _, _, _, pred, _, _, _, _ = predict(model, loader, device)
df["pred_delta_pIC50"] = pred
df.to_csv(base+"predict.csv", index=False)
```

## Contents of this folder

- `7EJ0.pdb` — raw source structure (alpha2A-AR + norepinephrine, GoA complex)
- `pocket.pdb` — alpha2A pocket, extracted from `7EJ0.pdb` (8 Å cutoff around the ligand)
- `norepinephrine_alpha2A.sdf` / `.pkl` — native ligand pose extracted from `7EJ0.pdb`, paired
  with `pocket.pdb`
- `ne_pocket.pdb` — generic adrenergic pocket (copied from `main/`), used for the two
  `*_generic` complexes below
- `norepinephrine_generic.sdf` / `.pkl`, `epinephrine_generic.sdf` / `.pkl` — pre-existing docked
  poses (copied from `main/`), paired with `ne_pocket.pdb`
- `predict.csv` — the 3×3 pairwise prediction results (§4)
