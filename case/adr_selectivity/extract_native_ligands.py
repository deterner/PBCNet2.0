"""
Extract the natively co-crystallized epinephrine (ALE) / norepinephrine (E5E)
ligand from each adrenergic-receptor source PDB structure, using RDKit
template-based bond-order assignment (against ALE_ideal.sdf / LNR_ideal.sdf)
so the extracted ligand is a valid RDKit mol with correct bonds -- required
before Graph2pickle.graph_save can build a graph.

The source PDB files (8THL.pdb, 7YMH.pdb, 7EJ0.pdb, 7BTS.pdb, 7BU6.pdb,
4LDO.pdb, 9IJE.pdb) and the ideal-ligand templates (ALE_ideal.sdf,
LNR_ideal.sdf) live on the `alpha2a_cpu` branch under main/ and main/adr/ --
they are not duplicated into this folder. This script is included for
reproducibility/reference; the outputs it produces (<ligand>.sdf per
subtype) are already committed alongside each subtype's pocket.pdb here.

Usage (from a checkout of the alpha2a_cpu branch):
    python extract_native_ligands.py --base /path/to/alpha2a_cpu/main --out /path/to/this/folder
"""
import argparse
import os
from rdkit import Chem
from rdkit.Chem import AllChem

# (folder, pdb_id, chain, resnum, resname, template_sdf, out_name)
TARGETS = [
    ("adr/alpha1A_epi",    "8THL", "R", 401, "ALE", "ALE_ideal.sdf", "epinephrine"),
    ("adr/alpha1A_norepi", "7YMH", "A", 401, "E5E", "LNR_ideal.sdf", "norepinephrine"),
    ("adr/alpha2A_norepi", "7EJ0", None, None, "E5E", "LNR_ideal.sdf", "norepinephrine"),
    ("adr/beta1_epi",      "7BTS", None, None, "ALE", "ALE_ideal.sdf", "epinephrine"),
    ("adr/beta1_norepi",   "7BU6", None, None, "E5E", "LNR_ideal.sdf", "norepinephrine"),
    ("adr/beta2_epi",      "4LDO", None, None, "ALE", "ALE_ideal.sdf", "epinephrine"),
    ("adr/beta3_epi",      "9IJE", None, None, "ALE", "ALE_ideal.sdf", "epinephrine"),
]


def extract_ligand_block(pdb_path, resname, chain=None, resnum=None):
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            rn = line[17:20].strip()
            ch = line[21].strip()
            rnum = line[22:26].strip()
            if rn != resname:
                continue
            if chain is not None and ch != chain:
                continue
            if resnum is not None and rnum != str(resnum):
                continue
            lines.append(line)
    lines.append("END\n")
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, help="path to alpha2a_cpu branch's main/ folder (source PDBs + templates)")
    parser.add_argument("--out", required=True, help="path to this adr_selectivity folder (destination for <ligand>.sdf files)")
    args = parser.parse_args()

    for folder, pdbid, chain, resnum, resname, template_name, out_name in TARGETS:
        subtype = folder.split("/", 1)[1]  # e.g. "alpha1A_epi"
        pdb_path = os.path.join(args.base, folder, f"{pdbid}.pdb")
        template_path = os.path.join(args.base, template_name)
        out_path = os.path.join(args.out, subtype, f"{out_name}.sdf")

        block = extract_ligand_block(pdb_path, resname, chain, resnum)
        raw_mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=True)
        if raw_mol is None:
            print(f"[FAIL] {subtype}: could not parse HETATM block for {resname}")
            continue

        template = Chem.MolFromMolFile(template_path)
        try:
            fixed_mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
        except Exception as e:
            print(f"[FAIL] {subtype}: template assignment failed: {e}")
            continue

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        Chem.MolToMolFile(fixed_mol, out_path)
        print(f"[OK] {subtype}: {resname} ({pdbid}) -> {out_path}  atoms={fixed_mol.GetNumAtoms()}")


if __name__ == "__main__":
    main()
