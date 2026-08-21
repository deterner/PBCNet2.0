# Adrenergic receptor subtype selectivity: epinephrine vs. norepinephrine

Uses PBCNet2.0 to predict relative binding affinity (Δ pIC50) between
epinephrine and norepinephrine across several adrenergic receptor (AR)
subtypes, each with its own natively resolved binding pocket, instead of a
single shared pocket (as in `case/toy_data`).

## Source structures

Each subtype folder holds the receptor pocket and the natively
co-crystallized catecholamine ligand extracted from a real PDB structure
(source PDBs and ideal-ligand templates live on the `alpha2a_cpu` branch,
not duplicated here — see `extract_native_ligands.py`):

| Folder | PDB ID | Subtype | Ligand | Complex |
|---|---|---|---|---|
| `alpha1A_epi` | 8THL | alpha1A-AR | epinephrine (ALE) | Gq complex |
| `alpha1A_norepi` | 7YMH | alpha1A-AR | norepinephrine (E5E) | Nb29-miniGsq complex |
| `alpha2A_norepi` | 7EJ0 | alpha2A-AR | norepinephrine (E5E) | GoA complex |
| `beta1_epi` | 7BTS | beta1-AR | epinephrine (ALE) | nanobody complex |
| `beta1_norepi` | 7BU6 | beta1-AR | norepinephrine (E5E) | nanobody complex |
| `beta2_epi` | 4LDO | beta2-AR | epinephrine (ALE) | nanobody complex |
| `beta3_epi` | 9IJE | beta3-AR | epinephrine (ALE) | — |
| `alpha2A_docked_ref` | 7EJA | alpha2A-AR | both, docked into one shared pocket | reference from the original `main/run_predict.py` run |

Each subtype's `pocket.pdb` was extracted (8 Å around the ligand,
biopython `Bio.PDB`) directly from its own source structure — so, unlike
`alpha2A_docked_ref`, epinephrine and norepinephrine at the same subtype
come from *different* PDB entries (different intracellular binding
partners), not the same shared pocket.

**Coverage gap:** only alpha1A and beta1 have both ligands natively
resolved. alpha2A has no native epinephrine structure, and beta2/beta3
have no native norepinephrine structure — the `alpha2A_docked_ref` entry
covers alpha2A via docking instead of native structures for both ligands
in the same pocket.

## Pipeline

1. `extract_native_ligands.py` — pulls the HETATM block for the ALE/E5E
   residue out of each source PDB and assigns bond orders via RDKit's
   `AllChem.AssignBondOrdersFromTemplate` against `ALE_ideal.sdf` /
   `LNR_ideal.sdf` (both molecules are epinephrine/norepinephrine; `E5E`
   is simply a newer PDB chemical-component ID for the same
   norepinephrine molecule as the `LNR_ideal.sdf` template).
2. `Graph2pickle.graph_save` (`model_code`) builds each `<ligand>.pkl`
   graph from the ligand SDF + `pocket.pdb`.
3. `run_adr_predict.py` builds the full pairwise `adr_predict.csv` (every
   entry vs. every entry) and runs `PBCNet2.pth` over it.

## Headline results (`adr_predict.csv`)

Same-subtype epinephrine vs. norepinephrine (Δ pIC50 = epi − norepi):

| Subtype | Δ pIC50 | Reading |
|---|---|---|
| alpha1A-AR | +0.836 | epinephrine favored (~7x) |
| beta1-AR | +1.743 | epinephrine favored (~55x) |
| alpha2A-AR (docked, shared pocket) | +1.060 | epinephrine favored (~11x) |

The model predicts epinephrine binds more tightly than norepinephrine at
every subtype where a same-receptor comparison is possible, strongest at
beta1-AR.

**Caveats:** these are relative-ranking predictions from a pairwise model,
not absolute affinities. The alpha1A/beta1 same-subtype comparisons use
different PDB entries (different complex partners) for each ligand, so
pocket geometry may not be identical between the epi- and norepi-bound
structures of the same subtype. alpha2A/beta2/beta3 lack a native
same-pocket comparison for both ligands.
