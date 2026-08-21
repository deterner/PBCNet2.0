# Adrenergic receptor selectivity, cleaner version: shared pocket per subtype

Follow-up to `case/adr_selectivity/`. There, alpha1A-AR and beta1-AR each
had epinephrine and norepinephrine bound in *different* PDB entries
(different complex partners), so the "same subtype" comparison wasn't
strictly the same structure. This folder fixes that: for each subtype, both
ligands are now placed in one shared pocket.

## Method (`merge_pockets.py`)

1. Superimpose the norepinephrine-bound receptor chain onto the
   epinephrine-bound receptor chain — biopython `Superimposer`, rigid-body
   fit on matched Cα atoms (by residue number) of the receptor chain only.
2. Apply that same rotation/translation to the already bond-order-correct
   norepinephrine ligand coordinates (from `case/adr_selectivity/`), moving
   norepinephrine into the epinephrine structure's reference frame.
3. Re-extract a single 8 Å pocket from the epinephrine (reference)
   structure, using the union of both ligands' positions — so both ligands
   now share the exact same `pocket.pdb`.

| Subtype | Reference (epi) | Superposed (norepi) | Matched Cα atoms | Superposition RMSD |
|---|---|---|---|---|
| alpha1A-AR | 8THL | 7YMH | 243 | 0.981 Å |
| beta1-AR | 7BTS | 7BU6 | 449 | 0.029 Å |

beta1-AR's two structures are nearly identical conformations (0.029 Å),
so that comparison is about as clean as it gets. alpha1A-AR's fit is
looser (0.981 Å, different complex partners: Gq vs. Nb29-miniGsq) but
still a good rigid-body superposition.

## Results (`merged_predict.csv`)

Same-pocket epinephrine vs. norepinephrine (Δ pIC50 = epi − norepi):

| Subtype | Mixed-structure (adr_selectivity) | Shared-pocket (this folder) |
|---|---|---|
| alpha1A-AR | +0.836 | **+1.047** |
| beta1-AR | +1.743 | **+1.756** |

Both subtypes still predict epinephrine binding more tightly than
norepinephrine, and the shared-pocket result is close to the original
mixed-structure estimate in both cases (beta1-AR barely moves, consistent
with its near-zero superposition RMSD; alpha1A-AR shifts up somewhat).
This is a useful sanity check: the earlier mixed-structure comparison
wasn't badly confounded by using different PDB entries per ligand.

Cross-subtype: beta1-AR is still predicted as the tighter-binding pocket
for both ligands relative to alpha1A-AR (beta1_epi − alpha1A_epi = +0.891;
beta1_norepi − alpha1A_norepi = +0.661).

## Caveat

The pocket-merge is a rigid-body superposition of the receptor structure,
not a re-docking — norepinephrine's pose is carried over from its own
(different) PDB entry and may not perfectly reflect how it would actually
sit in the epinephrine structure's exact side-chain rotamers, especially
at alpha1A-AR where the fit RMSD is higher.
