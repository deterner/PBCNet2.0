import os
from Graph2pickle import graph_save

epi_sdf = "epinephrine_docked_v21.sdf"
nor_sdf = "norepinephrine.sdf"
pocket = "ne_pocket.pdb"

print("에피네프린 SDF 존재?", os.path.exists(epi_sdf))
print("노르에피네프린 SDF 존재?", os.path.exists(nor_sdf))
print("포켓 존재?", os.path.exists(pocket))

graph_save(epi_sdf, pocket, "epinephrine.pkl")
graph_save(nor_sdf, pocket, "norepinephrine.pkl")

print("epinephrine.pkl 생성됨?", os.path.exists("epinephrine.pkl"))
print("norepinephrine.pkl 생성됨?", os.path.exists("norepinephrine.pkl"))