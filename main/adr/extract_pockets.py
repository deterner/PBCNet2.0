"""
다른 아드레날린 수용체 서브타입들의 포켓 추출 (epinephrine/norepinephrine 실제 결합 구조 사용).

RCSB PDB에서 받은 구조들은 리간드(epinephrine=ALE, norepinephrine=E5E)가
이미 결합된 상태라서 별도 도킹 없이 구조 안의 리간드 좌표를 그대로 사용한다.
포켓 추출 로직은 case/try.ipynb, main/run_local.py의 extract() 함수와 동일
(리간드 non-H 원자 기준 8.0 Å 이내 residue를 선택, PDBIO로 저장).

실행: conda activate pbcnet && python main/adr/extract_pockets.py
"""
import os
import numpy as np
from scipy.spatial import distance_matrix
from Bio.PDB import PDBParser, PDBIO, Select

# (폴더명, PDB ID, 리간드 residue 코드, 설명)
TARGETS = [
    ("alpha1A_epi",    "8THL", "ALE", "alpha1A-AR + epinephrine (Gq complex)"),
    ("alpha1A_norepi", "7YMH", "E5E", "alpha1A-AR + norepinephrine (Nb29-miniGsq)"),
    ("alpha2A_norepi", "7EJ0", "E5E", "alpha2A-AR + norepinephrine (GoA complex)"),
    ("beta1_epi",      "7BTS", "ALE", "beta1-AR + epinephrine (nanobody)"),
    ("beta1_norepi",   "7BU6", "E5E", "beta1-AR + norepinephrine (nanobody)"),
    ("beta2_epi",      "4LDO", "ALE", "beta2-AR + epinephrine (nanobody)"),
    ("beta3_epi",      "9IJE", "ALE", "beta3-AR + epinephrine"),
]

base_dir = os.path.dirname(os.path.abspath(__file__))
CUTOFF = 8.0


def extract_pocket(pdb_path, ligand_resname, out_path, cutoff=CUTOFF):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)

    ligand_positions = np.array([
        atom.get_coord()
        for res in structure.get_residues()
        if res.get_resname() == ligand_resname
        for atom in res
        if atom.element != "H"
    ])
    if len(ligand_positions) == 0:
        raise ValueError(f"ligand '{ligand_resname}' not found in {pdb_path}")

    class ResidueSelect(Select):
        def accept_residue(self, residue):
            # 표준 아미노산(ATOM) residue만 후보로 (리간드/이온/지질/water/버퍼 등 HETATM 제외)
            if residue.id[0] != " ":
                return 0
            residue_positions = np.array([
                atom.get_coord() for atom in residue.get_atoms()
                if "H" not in atom.get_id()
            ])
            if residue_positions.ndim < 2 or residue_positions.shape[0] == 0:
                return 0
            min_dis = np.min(distance_matrix(residue_positions, ligand_positions))
            return 1 if min_dis < cutoff else 0

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, ResidueSelect())


if __name__ == "__main__":
    for name, pdbid, ligcode, desc in TARGETS:
        folder = os.path.join(base_dir, name)
        pdb_path = os.path.join(folder, f"{pdbid}.pdb")
        out_path = os.path.join(folder, "pocket.pdb")
        extract_pocket(pdb_path, ligcode, out_path)
        print(f"[{name}] {pdbid} ({desc}) ligand={ligcode} -> {out_path}")
