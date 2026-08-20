import os
from Graph2pickle import graph_save

print("현재 작업 디렉토리:", os.getcwd())
print("리간드 파일 존재?", os.path.exists("epinephrine_final1.sdf"))
print("포켓 파일 존재?", os.path.exists("pocket.pdb"))

graph_save(
    ligand_file="epinephrine_final1.sdf",
    pock_file="pocket.pdb",
    pickle_save="epinephrine.pkl"
)

print("pkl 생성됨?", os.path.exists("epinephrine.pkl"))