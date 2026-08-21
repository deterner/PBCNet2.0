"""
beta2-AR (4LDO) and beta3-AR (9IJE) have no native norepinephrine-bound PDB
structure (only epinephrine). Unlike alpha1A/beta1 (see merge_pockets.py --
same subtype, two PDB entries, CA matched by residue number), there is no
norepinephrine structure of beta2 or beta3 itself to superimpose.

This script borrows norepinephrine's pose from beta1-AR's native
norepinephrine structure (7BU6) -- the closest available source, same
beta-adrenergic receptor family -- and cross-superimposes it onto each of
beta2 and beta3's own epinephrine-bound structure. Because beta1/beta2/beta3
are different genes (not the same residue numbering), the CA correspondence
can't be found by residue number like the alpha1A/beta1 case; instead this
uses biopython's CEAligner (structure-based alignment, sequence-number
independent) to align the two receptor chains.

Because this borrows a pose from a genuinely different receptor subtype,
treat the resulting norepinephrine placement as a rougher approximation
than the same-subtype alpha1A/beta1 merge -- see README.md for the
whole-chain CEAligner RMSD (4.7 A for beta2, 3.7 A for beta3 -- inflated by
flexible loops/termini in the global chain alignment) and the local
pocket-region sanity check (epinephrine/norepinephrine centroid distance:
2.5 A / 3.9 A, consistent with two similar catecholamines sharing one
orthosteric site).

Source PDBs and native ligand SDFs live on the alpha2a_cpu branch /
case/adr_selectivity/ -- not duplicated here.

Usage:
    python merge_pockets_cross_subtype.py --wt /path/to/alpha2a_cpu/main --out /path/to/this/folder
"""
import argparse
import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.cealign import CEAligner
from scipy.spatial import distance_matrix
from rdkit import Chem

CUTOFF = 8.0

NOREPI_SOURCE_FOLDER = "adr/beta1_norepi"
NOREPI_SOURCE_PDB = "7BU6"
NOREPI_SOURCE_CHAIN = "A"

# subtype: (ref_folder, ref_pdb, ref_chain)
TARGETS = {
    "beta2": ("adr/beta2_epi", "4LDO", "A"),
    "beta3": ("adr/beta3_epi", "9IJE", "R"),
}


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def write_single_chain(structure, chain_id, out_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ChainSelect(chain_id))


def ca_atoms_in_order(structure):
    atoms = []
    for chain in structure[0]:
        for res in chain:
            if res.id[0] == " " and "CA" in res:
                atoms.append(res["CA"])
    return atoms


def main():
    parser_args = argparse.ArgumentParser()
    parser_args.add_argument("--wt", required=True, help="path to alpha2a_cpu branch's main/ folder (source PDBs + native ligand SDFs)")
    parser_args.add_argument("--out", required=True, help="path to this adr_selectivity_merged folder")
    args = parser_args.parse_args()

    parser = PDBParser(QUIET=True)

    norepi_full = parser.get_structure("norepi_full", os.path.join(args.wt, NOREPI_SOURCE_FOLDER, f"{NOREPI_SOURCE_PDB}.pdb"))
    norepi_chain_pdb = os.path.join(args.out, "beta1_norepi_chain_only.pdb")
    write_single_chain(norepi_full, NOREPI_SOURCE_CHAIN, norepi_chain_pdb)

    for subtype, (ref_folder, ref_pdb, ref_chain) in TARGETS.items():
        out_dir = os.path.join(args.out, subtype)
        os.makedirs(out_dir, exist_ok=True)

        ref_full = parser.get_structure("ref_full", os.path.join(args.wt, ref_folder, f"{ref_pdb}.pdb"))
        ref_chain_pdb = os.path.join(out_dir, f"{ref_pdb}_chain_only.pdb")
        write_single_chain(ref_full, ref_chain, ref_chain_pdb)

        ref_iso = parser.get_structure("ref_iso", ref_chain_pdb)
        mob_before = parser.get_structure("mob_before", norepi_chain_pdb)
        mob_after = parser.get_structure("mob_after", norepi_chain_pdb)

        before_ca = ca_atoms_in_order(mob_before)

        aligner = CEAligner()
        aligner.set_reference(ref_iso)
        aligner.align(mob_after, transform=True)
        print(f"[{subtype}] CEAligner RMSD (beta1_norepi chain -> {ref_pdb} chain) = {aligner.rms:.3f} A")

        after_ca = ca_atoms_in_order(mob_after)
        assert len(before_ca) == len(after_ca)

        # recover the rigid transform CEAligner applied, via before/after CA correspondence
        sup = Superimposer()
        sup.set_atoms(after_ca, before_ca)
        rot, tran = sup.rotran

        # apply that transform to the already bond-order-correct norepinephrine ligand coordinates
        norepi_sdf_in = os.path.join(args.wt, NOREPI_SOURCE_FOLDER, "norepinephrine.sdf")
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
        io.set_structure(ref_full)
        pocket_out = os.path.join(out_dir, "pocket.pdb")
        io.save(pocket_out, ResidueSelect())

        print(f"[{subtype}] wrote {epi_sdf_out}, {norepi_sdf_out}, {pocket_out}")


if __name__ == "__main__":
    main()
