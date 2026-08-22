# PBCNet2.0 relative-affinity prediction: alpha2A-adrenergic receptor + epinephrine (9CBL)

Self-contained example, in the same style as [`case/alpha2A_norepi/`](../alpha2A_norepi/README.md)
(read that first — this folder follows its exact methodology and does not modify it). This folder
adds a **receptor-specific epinephrine pose** using PDB
[9CBL](https://www.rcsb.org/structure/9CBL) ("Cryo-EM structure of epinephrine-bound
alpha-2A-adrenergic receptor in complex with heterotrimeric Gi-protein", Lou et al. 2024,
*Exp. Mol. Med.* 56:1952-1966,
[doi:10.1038/s12276-024-01296-x](https://doi.org/10.1038/s12276-024-01296-x)) — closing the
"no epinephrine-bound alpha2A structure available" limitation noted in `alpha2A_norepi/README.md`.
Everything needed to inspect or rerun this analysis lives in this folder; paths in `predict.csv`
are relative to it. Nothing in `case/alpha2A_norepi/`, `main/adr/`, or any other existing file was
changed to produce this — the reference complexes used below are unmodified copies.

## 1. Objective

Same question as `alpha2A_norepi/`, extended: with a *receptor-specific* epinephrine pose now
available (not just a generic-pocket docked pose), does PBCNet2.0's predicted Δ pIC50 between
epinephrine and norepinephrine still favor epinephrine, and how sensitive is the prediction to
pocket/pose provenance across four complexes instead of three?

## 2. Methodological framework

Same model and pairwise-subtraction mechanism as described in `alpha2A_norepi/README.md` §2.
Four complexes were built and compared pairwise (4×4 = 16 predictions):

| label | files in this folder | ligand source | pocket |
|---|---|---|---|
| `epinephrine_alpha2A` | `epinephrine_alpha2A.sdf` / `.pkl` | co-crystallized pose, extracted directly from `9CBL.pdb` (real, not docked) | `pocket.pdb` — alpha2A-specific, extracted from 9CBL (Gi complex) |
| `norepinephrine_alpha2A` | `norepinephrine_alpha2A.sdf` / `.pkl` | unmodified copy from `case/alpha2A_norepi/` — native pose extracted from `7EJ0.pdb` | `pocket_alpha2A_norepi.pdb` — unmodified copy of `case/alpha2A_norepi/pocket.pdb` (alpha2A-specific, extracted from 7EJ0, GoA complex) |
| `norepinephrine_generic` | `norepinephrine_generic.sdf` / `.pkl` | unmodified copy from `case/alpha2A_norepi/` — pre-existing docked pose | `ne_pocket.pdb` — unmodified copy; generic adrenergic pocket |
| `epinephrine_generic` | `epinephrine_generic.sdf` / `.pkl` | unmodified copy from `case/alpha2A_norepi/` — pre-existing docked pose | `ne_pocket.pdb` — unmodified copy; generic adrenergic pocket |

Note `epinephrine_alpha2A` and `norepinephrine_alpha2A` come from two different structures
(9CBL/Gi vs. 7EJ0/GoA) — same receptor and orthosteric pocket, but resolved in complex with a
different G protein and therefore a potentially different overall conformational state; see
Limitations.

## 3. Steps performed

1. **Pocket extraction**: called `extract_pocket()` from `main/adr/extract_pockets.py`
   (imported, unmodified — that file's `TARGETS` list was *not* edited) directly against
   `9CBL.pdb` for residue `ALE` (epinephrine), same 8.0 Å non-H-atom cutoff as every other
   pocket in this repo, saved as `pocket.pdb` in this folder.
2. **Native ligand extraction**: pulled the 13 `HETATM` lines for residue `ALE` out of
   `9CBL.pdb` directly (real crystal coordinates, no re-docking), parsed with
   `Chem.MolFromPDBBlock(sanitize=False)`, assigned bond orders by matching against a SMILES
   template (`CNCC(O)c1ccc(O)c(O)c1` — epinephrine differs from norepinephrine only by the
   N-methyl group) with `AllChem.AssignBondOrdersFromTemplate`, added explicit hydrogens with
   `Chem.AddHs(addCoords=True)`, wrote to `epinephrine_alpha2A.sdf`.
3. **Reference complexes**: copied `norepinephrine_alpha2A.{sdf,pkl}`, `pocket.pdb` (renamed here
   to `pocket_alpha2A_norepi.pdb` to avoid clashing with this folder's own `pocket.pdb`),
   `norepinephrine_generic.{sdf,pkl}`, `epinephrine_generic.{sdf,pkl}`, and `ne_pocket.pdb`
   unmodified from `case/alpha2A_norepi/`.
4. **Graph construction** (`Graph2pickle.graph_save`): built `epinephrine_alpha2A.pkl` from
   `epinephrine_alpha2A.sdf` + this folder's `pocket.pdb`. Note: DGL graph pickles are not
   guaranteed portable across DGL versions — if a `.pkl` here fails to unpickle in your
   environment, regenerate it from the corresponding `.sdf` + pocket PDB with `graph_save`.
5. **Pairwise CSV** (`predict.csv`): all 16 ordered pairs among the four complexes above, with
   dummy `Label`/`Label1`/`Label2` columns (unused at inference), and `dir_1`/`dir_2` pointing at
   the `.pkl` files in this folder.
6. **Inference**: loaded `PBCNet2.pth` (repo root) on CPU (`torch.load(..., map_location='cpu')`),
   ran `model_code/predict/predict.predict(model, loader, device)` over the CSV via
   `LeadOptDataset` + `DataLoader(collate_fn=collate_fn)`, wrote predictions back to `predict.csv`
   as `pred_delta_pIC50`.

## 4. Results

`pred_delta_pIC50` = predicted (pIC50 of `lig1`) − (pIC50 of `lig2`); higher = `lig1` binds
tighter. (The three columns/rows not involving `epinephrine_alpha2A` reproduce
`case/alpha2A_norepi/predict.csv` exactly — same unmodified input pickles, same checkpoint.)

| lig1 | lig2 | pred Δ pIC50 |
|---|---|---:|
| norepinephrine_alpha2A | norepinephrine_alpha2A | −0.080 |
| norepinephrine_alpha2A | norepinephrine_generic | +0.531 |
| norepinephrine_alpha2A | epinephrine_alpha2A | −1.120 |
| norepinephrine_alpha2A | epinephrine_generic | −0.908 |
| norepinephrine_generic | norepinephrine_alpha2A | −0.756 |
| norepinephrine_generic | norepinephrine_generic | −0.080 |
| norepinephrine_generic | epinephrine_alpha2A | −1.376 |
| norepinephrine_generic | epinephrine_generic | −1.052 |
| epinephrine_alpha2A | norepinephrine_alpha2A | +1.184 |
| epinephrine_alpha2A | norepinephrine_generic | +1.311 |
| epinephrine_alpha2A | epinephrine_alpha2A | −0.080 |
| epinephrine_alpha2A | epinephrine_generic | −0.319 |
| epinephrine_generic | norepinephrine_alpha2A | +0.895 |
| epinephrine_generic | norepinephrine_generic | +1.060 |
| epinephrine_generic | epinephrine_alpha2A | +0.260 |
| epinephrine_generic | epinephrine_generic | −0.080 |

Full output: [`predict.csv`](./predict.csv).

## 5. Interpretation

- **Epinephrine vs. norepinephrine**: the model consistently predicts epinephrine binds tighter
  than norepinephrine across all four cross-ligand comparisons — generic pocket (Δ ≈ +1.06 /
  −1.05), native alpha2A poses on both sides (Δ ≈ +1.18 / −1.12), and the two mixed
  native/generic pairings (Δ ≈ +1.31 / −1.38 and +0.89 / −0.91). Direction and rough magnitude
  are stable whether or not a receptor-specific pose is used for either ligand — this is the main
  new finding this folder adds over `alpha2A_norepi/`, where the epinephrine side was
  generic-pocket-only.
- **Pocket sensitivity**: the same ligand scores differently depending on which pocket structure
  it's paired with — e.g. norepinephrine, alpha2A-native vs. generic: Δ ≈ +0.53 / −0.76;
  epinephrine, alpha2A-native vs. generic: Δ ≈ −0.32 / +0.26 (not symmetric — see below). This
  indicates the model's output is sensitive to specific pocket geometry, as intended by design,
  though the pockets being compared also differ in ligand pose source (crystal vs. docked) and,
  for the two alpha2A pockets, in which G protein the receptor was captured with (see
  Limitations), so these deltas conflate pocket identity with pose quality and receptor
  conformational state.
- **Self-comparison ("diagonal") is not exactly zero**: comparing a complex against itself gives
  −0.080 for all four, rather than the mathematically expected 0 — the FNN head's fixed bias term
  surviving the `emb1 - emb2 = 0` subtraction, not a property of the ligand or pocket. Same as
  observed in `alpha2A_norepi/`.
- **Approximate, not exact, anti-symmetry**: `g(A,B)` and `-g(B,A)` are close but not identical
  (e.g. +1.184 vs. −1.120, or +0.260 vs. −0.319). The gap (~0.06–0.5 pIC50) is a rough indicator
  of the model's internal noise floor on this input distribution.

## 6. Limitations

- **The two alpha2A-native poses come from structures in different G-protein-coupled states.**
  `epinephrine_alpha2A` is from 9CBL (alpha2A-AR-Gi complex); `norepinephrine_alpha2A` is from
  7EJ0 (alpha2A-AR-GoA complex). Both are active-state agonist-bound structures of the same
  receptor and the same orthosteric pocket, but resolved with different intracellular binding
  partners, which can produce small differences in overall receptor conformation. Comparisons
  between the two `*_alpha2A` complexes therefore are not a perfectly matched same-structure
  comparison, though they are a closer match than either pose against the unrelated `*_generic`
  pocket.
- **"Generic" pocket provenance is not receptor-subtype-specific.** `ne_pocket.pdb` was built for
  an earlier, unrelated demo (see `main/run_predict.py`) and is not itself an alpha2A structure —
  comparisons involving it mix receptor identity with pocket/pose source and should not be read
  as "alpha2A vs. alpha2A."
- **Cross-pocket comparison is outside the model's intended use case.** PBCNet2.0 is designed and
  benchmarked for relative binding affinity *within a single target's pocket* (e.g. lead
  optimization series, FEP benchmarks). Comparing complexes built from two different protein
  structures (as done here for `*_alpha2A` vs. `*_generic`, and between the two `*_alpha2A`
  complexes) is not part of its validated use case; treat those specific deltas as exploratory,
  not benchmark-grade.
- **No experimental ground truth.** No pIC50/Ki values were supplied for norepinephrine or
  epinephrine at alpha2A — dummy zero labels were used since they only affect ranking-metric
  computation, not the model's Δ prediction. Predicted deltas are not validated against any
  measured affinity in this analysis.
- **Ligand protonation/bond-order assumption.** `epinephrine_alpha2A` was reconstructed from
  heavy-atom crystal coordinates using a SMILES template match
  (`AssignBondOrdersFromTemplate`); 9CBL's `FORMUL` record for epinephrine (`C9 H13 N O3`) is
  neutral, which the template preserves, but this is inferred from the deposited formula, not
  directly observed from the electron density. (`norepinephrine_alpha2A`, copied unmodified from
  `alpha2A_norepi/`, was built the same way against 7EJ0's `FORMUL` record, which indicates a
  protonated amine, `C8 H12 N O3 1+`.)
- **Arbitrary pocket cutoff.** The 8.0 Å residue-inclusion cutoff (from
  `main/adr/extract_pockets.py::extract_pocket`) is a fixed heuristic, not tuned per-target; a
  different cutoff could change predictions.
- **Single checkpoint, no fine-tuning, no uncertainty estimate.** This used the pretrained
  `PBCNet2.pth` as-is (no fine-tuning on adrenergic-receptor or catecholamine data, see
  `model_code/Finetune.py`/`run_fintune.sh` for that workflow), and PBCNet2.0 does not output a
  confidence interval — only the point-estimate deltas above.

## 7. Replication

Requires the `pbcnet` conda environment (see repo-root `README.md` / `CLAUDE.md` for install
steps) with `PBCNet2.pth` present at the repo root, on the `alpha2a_cpu` branch. Does not require
or modify anything in `main/adr/extract_pockets.py` beyond importing its `extract_pocket()`
function, nor anything in `case/alpha2A_norepi/`.

```bash
# from the repo root
git checkout alpha2a_cpu

# 1. download the raw structure (skip if 9CBL.pdb already exists in this folder)
curl -s -o case/alpha2A_epi/9CBL.pdb https://files.rcsb.org/download/9CBL.pdb
```

```python
# 2. extract the alpha2A pocket from 9CBL.pdb, reusing extract_pocket() as-is (no edits to
#    main/adr/extract_pockets.py or its TARGETS list)
import sys; sys.path.append("main/adr")
from extract_pockets import extract_pocket
extract_pocket("case/alpha2A_epi/9CBL.pdb", "ALE", "case/alpha2A_epi/pocket.pdb")

# 3. extract the native epinephrine ligand from 9CBL.pdb
from rdkit import Chem
from rdkit.Chem import AllChem

lines = [l for l in open("case/alpha2A_epi/9CBL.pdb") if l.startswith("HETATM") and "ALE" in l]
mol = Chem.MolFromPDBBlock("".join(lines) + "END\n", sanitize=False, removeHs=False)
template = Chem.MolFromSmiles("CNCC(O)c1ccc(O)c(O)c1")
mol2 = Chem.AddHs(AllChem.AssignBondOrdersFromTemplate(template, mol), addCoords=True)
Chem.SDWriter("case/alpha2A_epi/epinephrine_alpha2A.sdf").write(mol2)

# 4. build the DGL graph pickle for the new complex
sys.path.append(".")
from Graph2pickle import graph_save
base = "case/alpha2A_epi/"
graph_save(base+"epinephrine_alpha2A.sdf", base+"pocket.pdb", base+"epinephrine_alpha2A.pkl")

# 5. copy the unmodified reference complexes from case/alpha2A_norepi/ (skip if already present)
import shutil
src = "case/alpha2A_norepi/"
for fname, dest in [
    ("norepinephrine_alpha2A.sdf", "norepinephrine_alpha2A.sdf"),
    ("norepinephrine_alpha2A.pkl", "norepinephrine_alpha2A.pkl"),
    ("pocket.pdb",                 "pocket_alpha2A_norepi.pdb"),
    ("norepinephrine_generic.sdf", "norepinephrine_generic.sdf"),
    ("norepinephrine_generic.pkl", "norepinephrine_generic.pkl"),
    ("epinephrine_generic.sdf",    "epinephrine_generic.sdf"),
    ("epinephrine_generic.pkl",    "epinephrine_generic.pkl"),
    ("ne_pocket.pdb",               "ne_pocket.pdb"),
]:
    shutil.copy(src+fname, base+dest)

# 6. build the pairwise CSV and run inference
import os, torch, pandas as pd
from torch.utils.data import DataLoader
sys.path.append("model_code")
from model_code.Dataloader.dataloader import LeadOptDataset, collate_fn
from model_code.predict.predict import predict

pkl_paths = {
    "norepinephrine_alpha2A": base+"norepinephrine_alpha2A.pkl",
    "norepinephrine_generic": base+"norepinephrine_generic.pkl",
    "epinephrine_alpha2A":    base+"epinephrine_alpha2A.pkl",
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

- `9CBL.pdb` — raw source structure (alpha2A-AR + epinephrine, Gi complex)
- `pocket.pdb` — alpha2A pocket, extracted from `9CBL.pdb` (8 Å cutoff around the ligand)
- `epinephrine_alpha2A.sdf` / `.pkl` — native ligand pose extracted from `9CBL.pdb`, paired with
  `pocket.pdb`
- `norepinephrine_alpha2A.sdf` / `.pkl`, `pocket_alpha2A_norepi.pdb` — unmodified copies from
  `case/alpha2A_norepi/` (native norepinephrine pose + its alpha2A pocket from 7EJ0)
- `norepinephrine_generic.sdf` / `.pkl`, `epinephrine_generic.sdf` / `.pkl`, `ne_pocket.pdb` —
  unmodified copies from `case/alpha2A_norepi/` (generic-pocket docked poses)
- `predict.csv` — the 4×4 pairwise prediction results (§4)
