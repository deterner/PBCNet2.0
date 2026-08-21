"""
In case/adr_selectivity/, the alpha1A and beta1 epinephrine/norepinephrine
comparisons used pockets from two *different* PDB entries per subtype
(different complex partners), so they weren't a strictly apples-to-apples
comparison. This script merges each pair into ONE shared structure:

1. Superimpose the norepinephrine-bound receptor chain onto the
   epinephrine-bound receptor chain (biopython Superimposer, matched-CA
   rigid-body fit on the receptor chain only, matched by residue number).
2. Apply that same rotation/translation to the already-extracted,
   correctly bonded norepinephrine ligand coordinates -- norepinephrine
   now lives in the epinephrine structure's reference frame.
3. Re-extract a single unified 8 A pocket from the epinephrine (reference)
   structure, using the union of both ligands' positions -- so both
   ligands are now paired with the exact same pocket.pdb.

Source PDBs (8THL, 7YMH, 7BTS, 7BU6) and their extracted native ligand SDFs
live on the alpha2a_cpu branch / case/adr_selectivity/ -- not duplicated
here. This script is included for reproducibility/reference; its outputs
(pocket.pdb, epinephrine.sdf, norepinephrine.sdf per subtype) are already
committed in this folder.

Usage:
    python merge_pockets.py --wt /path/to/alpha2a_cpu/main --out /path/to/this/folder
"""
import argparse
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from scipy.spatial import distance_matrix
from rdkit import Chem

# subtype: (ref_folder, ref_pdb, ref_chain, mobile_folder, mobile_pdb, mobile_chain)
PAIRS = {
    "alpha1A": ("adr/alpha1A_epi", "8THL", "R", "adr/alpha1A_norepi", "7YMH", "A"),
    "beta1":   ("adr/beta1_epi", "7BTS", "A", "adr/beta1_norepi", "7BU6", "A"),
}

CUTOFF = 8.0


def ca_map(structure, chain_id):
    d = {}
    chain = structure[0][chain_id]
    for res in chain:
        if res.id[0] == " " and "CA" in res:
            d[res.id[1]] = res["CA"]
    return d


def main():
    parser_args = argparse.ArgumentParser()
    parser_args.add_argument("--wt", required=True, help="path to alpha2a_cpu branch's main/ folder (source PDBs + native ligand SDFs)")
    parser_args.add_argument("--out", required=True, help="path to this adr_selectivity_merged folder")
    args = parser_args.parse_args()

    parser = PDBParser(QUIET=True)

    for subtype, (ref_folder, ref_pdb, ref_chain, mob_folder, mob_pdb, mob_chain) in PAIRS.items():
        out_dir = os.path.join(args.out, subtype)
        os.makedirs(out_dir, exist_ok=True)

        ref_struct = parser.get_structure("ref", os.path.join(args.wt, ref_folder, f"{ref_pdb}.pdb"))
        mob_struct = parser.get_structure("mob", os.path.join(args.wt, mob_folder, f"{mob_pdb}.pdb"))

        ref_ca = ca_map(ref_struct, ref_chain)
        mob_ca = ca_map(mob_struct, mob_chain)
        common = sorted(set(ref_ca) & set(mob_ca))

        fixed_atoms = [ref_ca[i] for i in common]
        moving_atoms = [mob_ca[i] for i in common]

        sup = Superimposer()
        sup.set_atoms(fixed_atoms, moving_atoms)
        rot, tran = sup.rotran
        print(f"[{subtype}] matched {len(common)} CA atoms, superposition RMSD = {sup.rms:.3f} A")

        # transform the already-extracted norepinephrine ligand into the reference frame
        norepi_sdf_in = os.path.join(args.wt, mob_folder, "norepinephrine.sdf")
        mol = Chem.MolFromMolFile(norepi_sdf_in, removeHs=False)
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = np.array(conf.GetAtomPosition(i))
            new_pos = np.dot(pos, rot) + tran
            conf.SetAtomPosition(i, new_pos.tolist())

        norepi_sdf_out = os.path.join(out_dir, "norepinephrine.sdf")
        Chem.MolToMolFile(mol, norepi_sdf_out)

        # epinephrine is already native to the reference structure -- just copy it
        epi_sdf_in = os.path.join(args.wt, ref_folder, "epinephrine.sdf")
        epi_mol = Chem.MolFromMolFile(epi_sdf_in, removeHs=False)
        epi_sdf_out = os.path.join(out_dir, "epinephrine.sdf")
        Chem.MolToMolFile(epi_mol, epi_sdf_out)

        # re-extract ONE unified pocket from the reference structure, union of both ligand positions
        epi_positions = epi_mol.GetConformer().GetPositions()
        norepi_positions = mol.GetConformer().GetPositions()
        ligand_positions = np.concatenate([epi_positions, norepi_positions])

        class ResidueSelect(Select):
            def accept_residue(self, residue):
                if residue.id[0] != " ":
                    return 0
                residue_positions = np.array([
                    atom.get_coord() for atom in residue.get_atoms() if "H" not in atom.get_id()
                ])
                if residue_positions.ndim < 2 or residue_positions.shape[0] == 0:
                    return 0
                min_dis = np.min(distance_matrix(residue_positions, ligand_positions))
                return 1 if min_dis < CUTOFF else 0

        io = PDBIO()
        io.set_structure(ref_struct)
        pocket_out = os.path.join(out_dir, "pocket.pdb")
        io.save(pocket_out, ResidueSelect())

        print(f"[{subtype}] wrote {epi_sdf_out}, {norepi_sdf_out}, {pocket_out}")


if __name__ == "__main__":
    main()
