# Adrenergic receptor selectivity, cleaner version: shared pocket per subtype

Follow-up to `case/adr_selectivity/`. There, each subtype's epinephrine and
norepinephrine comparison used pockets from *different* PDB entries (or, for
beta2/beta3, was missing one ligand entirely), so the "same subtype"
comparison wasn't strictly the same structure. This folder fixes that: for
each of alpha1A-AR, beta1-AR, beta2-AR, and beta3-AR, both ligands are now
placed in one shared pocket.

Two different methods were needed, depending on whether a native
norepinephrine structure of that exact subtype exists.

## Method 1 — same subtype, two PDB entries (`merge_pockets.py`)

Used for alpha1A-AR (8THL epi / 7YMH norepi) and beta1-AR (7BTS epi / 7BU6
norepi), where a real norepinephrine-bound structure of the *same* receptor
exists.

1. Superimpose the norepinephrine-bound receptor chain onto the
   epinephrine-bound receptor chain — biopython `Superimposer`, rigid-body
   fit on matched Cα atoms (by residue number, since both PDB entries use
   the same receptor numbering) of the receptor chain only.
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

beta1-AR's two structures are nearly identical conformations (0.029 Å), so
that comparison is about as clean as it gets. alpha1A-AR's fit is looser
(0.981 Å, different complex partners: Gq vs. Nb29-miniGsq) but still a good
rigid-body superposition.

## Method 2 — cross-subtype pose transfer (`merge_pockets_cross_subtype.py`)

Used for beta2-AR (4LDO) and beta3-AR (9IJE), neither of which has a native
norepinephrine-bound structure at all — only epinephrine. There's no
same-subtype norepinephrine structure to superimpose, so instead this
borrows norepinephrine's pose from **beta1-AR's native norepinephrine
structure (7BU6)** — the closest available source, same beta-adrenergic
family — and cross-superimposes it onto each of beta2 and beta3's own
epinephrine-bound structure.

Because beta1/beta2/beta3 are different genes (different residue
numbering), the Cα correspondence can't be found by residue number like
Method 1. Instead this uses biopython's `CEAligner` (structure-based
alignment, independent of residue numbering) to align the receptor chains,
then recovers the exact rigid transform CEAligner applied (via a
before/after Cα correspondence trick) and applies it to the norepinephrine
ligand coordinates.

| Subtype | Norepinephrine source | Whole-chain CEAligner RMSD | Epi/norepi pocket-centroid distance |
|---|---|---|---|
| beta2-AR | beta1-AR (7BU6) | 4.714 Å | 2.505 Å |
| beta3-AR | beta1-AR (7BU6) | 3.681 Å | 3.921 Å |

The whole-chain RMSD is inflated by flexible loops/termini that don't
correspond well between different receptor genes — it is **not** a measure
of pocket-region accuracy. As a local sanity check, the epinephrine/
norepinephrine centroid distance after transfer (2.5 Å / 3.9 Å) is
consistent with two structurally similar catecholamines sharing one
orthosteric site, which suggests the transferred pose lands in a
physically reasonable place, but this is a rougher approximation than
Method 1 and should be read with more caution.

## Results (`merged_predict.csv`)

Same-pocket epinephrine vs. norepinephrine (Δ pIC50 = epi − norepi):

| Subtype | Mixed-structure (`adr_selectivity`) | Shared-pocket (this folder) | Method |
|---|---|---|---|
| alpha1A-AR | +0.836 | **+1.047** | same-subtype merge |
| beta1-AR | +1.743 | **+1.756** | same-subtype merge |
| beta2-AR | n/a (no native norepi) | **+0.291** | cross-subtype transfer |
| beta3-AR | n/a (no native norepi) | **+0.723** | cross-subtype transfer |

All four subtypes predict epinephrine binding more tightly than
norepinephrine at their own pocket. The same-subtype comparisons
(alpha1A, beta1) are close to their earlier mixed-structure estimates —
a useful sanity check that those weren't badly confounded. The
cross-subtype comparisons (beta2, beta3) show a smaller epi/norepi gap,
consistent with the extra approximation involved in transferring
norepinephrine's pose from a different receptor.

Cross-subtype affinity ranking (relative to alpha1A_epi, all values are
target − alpha1A_epi): beta1_epi (+0.891) > beta2_epi (+0.699) >
beta3_epi (−0.194) > alpha1A_epi (0) — beta1-AR is predicted as the
tightest-binding pocket for both ligands throughout.

## Caveats

- The pocket-merge is a rigid-body superposition of the receptor
  structure, not a re-docking — the transferred ligand's pose reflects
  the source structure's binding mode, not necessarily the exact
  side-chain rotamers of the target structure.
- Method 1 (same subtype) is the more reliable of the two; Method 2
  (cross-subtype, beta2/beta3) carries the added uncertainty of
  transferring a pose across different receptor genes and should be
  treated as a rougher estimate.
