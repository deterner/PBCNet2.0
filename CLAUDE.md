# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PBCNet2.0 is a research codebase (paper: "Advancing Ligand Binding Affinity Prediction with
Cartesian Tensor-Based Deep Learning") implementing a Cartesian-tensor GNN (built on TensorNet)
that predicts **relative** protein-ligand binding affinity: given two ligand poses (docked into the
same or different pockets), the model predicts the delta in binding affinity (pIC50/pKi) between
them, not an absolute value. `PBCNet2.pth` at the repo root is the pretrained model checkpoint.

There is no build system, package manifest, linter, or test suite in this repo — it is a set of
scripts run directly with `python`, structured for training, fine-tuning, and running inference.

## Environment setup

```bash
conda create --name pbcnet python=3.8
pip3 install torch torchvision torchaudio
pip install pandas packaging PyYAML pydantic scipy matplotlib rdkit networkx psutil tqdm dataloader scikit-learn bio
pip install dgl==1.0.2 -f https://data.dgl.ai/wheels/cu113/repo.html --no-deps
```

Tested on Ubuntu 18.04, CUDA 11.4. GPU is expected by default (`cuda:N` device strings throughout),
but the current branch (`alpha2a_cpu`) contains CPU-only adaptations in `main/` — those scripts
explicitly use `torch.device('cpu')` / `map_location=torch.device('cpu')`.

## Common commands

There are no test/lint/build commands — "running" this project means executing one of the pipeline
scripts below with `python`.

**Train from scratch:**
```bash
./model_code/run_train.sh
```
Wraps `model_code/train.py`. Key args: `--train_path` (CSV of pairwise training examples),
`--hidden_dim`, `--radius` (number of TensorNet message-passing layers), `--retrain` (continue
from an existing checkpoint), `--device` (GPU index).

**Fine-tune on the FEP benchmark:**
```bash
./model_code/run_fintune.sh
```
Wraps `model_code/Finetune.py`. Iterates over FEP1/FEP2 benchmark systems and reference-ligand
counts, fine-tuning a copy of `PBCNet2.pth` per system and logging Spearman/Pearson/Kendall
correlation and RMSE before/after fine-tuning to `results/finetune/<system>/ref<N>_results/`.

**Preprocess structures into graphs (required before training/fine-tuning/prediction):**
```bash
python Graph2pickle.py   # or main/Graph2pickle.py — identical copies, kept in sync manually
```
`graph_save(ligand_file, pock_file, pickle_save)` converts one ligand SDF + one pocket/protein PDB
into a DGL heterograph pickle (`.pkl`). Training/eval CSVs reference these `.pkl` paths directly
(columns `dir_1`/`dir_2` or `Ligand1`/`Ligand2` depending on the loader — see below).

**Run a prediction (recommended entry point):** `case/try.ipynb` — step-by-step notebook example.
For a scripted equivalent see `main/run_predict.py` or `main/run_local.py`, which both:
1. call `graph_save` on ligand SDF(s) + pocket PDB to produce `.pkl` graphs,
2. build a pairwise CSV (`lig1`, `lig2`, `Label`, `Label1`, `Label2`, `dir_1`, `dir_2`),
3. `torch.load` `PBCNet2.pth`, wrap the CSV in `LeadOptDataset` + `DataLoader(collate_fn=collate_fn)`,
4. call `model_code/predict/predict.predict(model, loader, device)` and write predictions back to
   the CSV as delta pIC50 (`lig1` value minus `lig2` value).

**Reproduce paper results:** notebooks/scripts under `Result_in_paper/<section>/`. Paths inside
them are hardcoded to the original authors' environment and need updating before use.

## Architecture

**Pairwise, not pointwise.** The model (`model_code/models/readout.py::PBCNetv2`) never scores a
single pose — `forward(g1, g2)` embeds two ligand-pocket graphs independently through a shared
`TensorNet` encoder, subtracts their embeddings (`emb1 - emb2`), and feeds the difference through
an FNN readout to predict the affinity delta. It also computes the same thing in reverse
(`emb2 - emb1`) as `out_put_neg`, and training losses are applied to both directions
symmetrically (see `train.py`/`Finetune.py`: `loss_func(logits, label) + loss_func(logits_neg, -label)`).
This anti-symmetry is central to the design — any change to the readout must preserve it.

**Graph encoder is TensorNet** (`model_code/models/tensornet.py`), a Cartesian-tensor equivariant
GNN (rank-0/1/2 tensor messages, O(3) or SO(3) equivariance group configurable). Ligand and pocket
atoms live in one DGL heterograph with node type `atom` and edge type `int`; protein atoms are
masked out of the final readout via `g.nodes['atom'].data['type']` (see `_readout` in
`readout.py`) — only ligand-atom embeddings are pooled, the protein pocket only contributes through
message passing.

**Graph construction** (`Graph2pickle.py`) builds this heterograph from a ligand SDF + a pocket/
protein PDB using RDKit (ligand featurization: chirality, degree, hybridization, ring membership,
etc. via `allowable_features`) and BioPython's `PDBParser` (pocket residues). Output is a pickled
DGL graph consumed by the dataloader — this preprocessing step is a hard prerequisite for every
downstream script; nothing reads SDF/PDB directly at train/predict time.

**Data flow / loaders** (`model_code/Dataloader/dataloader.py`): several near-duplicate
`collate_fn*` variants exist because the two example CSV column conventions used across the repo
are *not* interchangeable:
- root-level / `main/` scripts use `dir_1`/`dir_2` + `Label`/`Label1`/`Label2` and `collate_fn`
  (paths point directly at `.pkl` files).
- `model_code/train.py`/`Finetune.py`/paper-reproduction code use `Ligand1`/`Ligand2` +
  `Lable`/`Lable1`/`Lable2` (note the misspelling — it's intentional/consistent, not a typo to
  fix) and derive `.pkl` paths by string manipulation relative to `code_path` (e.g.
  `collate_fn_fep_ft` rewrites `.../pose/<system>/<lig>.sdf` style paths into
  `data/FEP/pose_graph/<system>/<lig>_dgl_group.pkl`). When adding a new data source, match the
  convention of whichever pipeline (root/`main` vs `model_code`) you're extending rather than
  introducing a fifth variant.
- `Label`/`Lable` is the pairwise delta (what the model is trained/scored on); `Label1`/`Lable1`
  and `Label2`/`Lable2` are the two ligands' absolute experimental values, used only for
  reconstructing absolute predictions and ranking metrics (Spearman/Pearson/Kendall) after
  inference, not as model inputs.

**`code_path` bootstrap pattern**: nearly every script under `model_code/` computes
`code_path = os.path.dirname(os.path.abspath(__file__))` then does `sys.path.append` and
`.rsplit("/", 1)[0]` gymnastics to locate sibling packages (`Dataloader`, `models`, `predict`,
`utilis`) and the repo root, because there's no installed package / `PYTHONPATH` setup. Scripts
under `model_code/` must be run with `model_code/` as the effective import root (as the
`run_*.sh` wrappers do via `$script_dir`); scripts under `main/` instead `sys.path.append` both
the repo root and `model_code/` explicitly and import via `model_code.Dataloader...` /
`model_code.predict...`. Don't mix the two import styles in one script.

**Directory map:**
- `model_code/models/` — `tensornet.py` (encoder), `readout.py` (`PBCNetv2`, the pairwise head),
  `utils.py`/`wrappers.py` (RBF, activations, distance helpers).
- `model_code/Dataloader/dataloader.py` — CSV-driven `Dataset` + `collate_fn` variants (see above).
- `model_code/predict/predict.py` — inference loops, including FEP-benchmark-specific
  `test_fep`/`test_fep_nobond` helpers that also compute ranking correlations.
- `model_code/utilis/` — loss functions (`function.py`), Noam-style LR scheduler
  (`scheduler.py`), weight init (`initial.py`), pickle load helpers (`utilis.py`), misc/logging
  (`trick.py`, includes the `Writer` used for all training/fine-tune logs).
- `model_code/Selection.py` — standalone active-learning-style ligand selection/fine-tuning driver
  (separate from `train.py`/`Finetune.py`).
- `data/` — benchmark datasets: `FEP/` (FEP1/FEP2 + fine-tune splits), `Mutation/`, `SAR-Diff/`,
  `Selection/`, `F-Opt/`. Training data (`training_clip_862W.csv`/`.zip`) is not checked in —
  download from the Zenodo record linked in `README.md`.
- `Result_in_paper/` — per-paper-section notebooks/scripts and their saved logs
  (`Finetune_results/<target>/ref*/N_record.txt`) reproducing published figures/tables.
- `case/` — `try.ipynb`, the recommended end-to-end usage example, plus `toy_data/` fixtures.
- `main/` — CPU-oriented, adrenergic-receptor-specific scripts on this branch (`alpha2a_cpu`):
  pocket extraction, graph generation, and prediction for epinephrine/norepinephrine-type ligands.

## Agent skills

### Issue tracker

Issues live as GitHub issues on `deterner/PBCNet2.0` (this fork, `origin`; not `upstream`). See `docs/agents/issue-tracker.md`.

### Domain docs

Single-context: `CONTEXT.md` + `docs/adr/` at the repo root (once created). See `docs/agents/domain.md`.
