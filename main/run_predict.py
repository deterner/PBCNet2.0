"""
PBCNet2.0 CPU 예측 스크립트 (adrenaline receptor: epinephrine vs norepinephrine)

준비된 파일 (main/ 폴더):
  - epinephrine_docked_v21.sdf, norepinephrine.sdf : 리간드 (도킹 포즈)
  - ne_pocket.pdb                                    : 포켓 (두 리간드 공통)
  - ../PBCNet2.pth                                    : 모델 가중치 (repo 루트)

README.md 4번 항목, case/try.ipynb 예측 파트를 참고해 작성.
실행: conda activate pbcnet && python main/run_predict.py
"""
import os
import sys
import torch
import pandas as pd
from torch.utils.data import DataLoader

# ---- 경로 설정 ----
base_dir = os.path.dirname(os.path.abspath(__file__))   # .../PBCNet2.0/main
repo_root = os.path.dirname(base_dir)                     # .../PBCNet2.0
sys.path.append(repo_root)
sys.path.append(os.path.join(repo_root, "model_code"))

from Graph2pickle import graph_save
from model_code.Dataloader.dataloader import LeadOptDataset, collate_fn
from model_code.predict.predict import predict

# ---- 1. 리간드 -> 그래프(.pkl) 생성 ----
# (name, sdf 파일, 포켓 파일)
ligands = [
    ("epinephrine", "epinephrine_docked_v21.sdf", "ne_pocket.pdb"),
    ("norepinephrine", "norepinephrine.sdf", "ne_pocket.pdb"),
]

pkl_paths = {}
for name, sdf, pocket in ligands:
    sdf_path = os.path.join(base_dir, sdf)
    pocket_path = os.path.join(base_dir, pocket)
    pkl_path = os.path.join(base_dir, f"{name}.pkl")
    graph_save(sdf_path, pocket_path, pkl_path)
    pkl_paths[name] = pkl_path
    print(f"[graph] {name}: {pkl_path}")

# ---- 2. pairwise 입력 csv 생성 ----
# 실험값(pIC50)을 모르므로 Label은 0으로 둔다 (예측값에는 영향 없음, 있으면 채워서 비교용으로 사용 가능)
names = list(pkl_paths.keys())
N1, N2, D1, D2, L = [], [], [], [], []
for n1 in names:
    for n2 in names:
        N1.append(n1)
        D1.append(pkl_paths[n1])
        N2.append(n2)
        D2.append(pkl_paths[n2])
        L.append(0.0)

predict_csv_path = os.path.join(base_dir, "predict.csv")
pd.DataFrame({
    "lig1": N1, "lig2": N2,
    "Label": L, "Label1": L, "Label2": L,
    "dir_1": D1, "dir_2": D2,
}).to_csv(predict_csv_path, index=False)

# ---- 3. CPU로 모델 로드 및 추론 ----
device = torch.device("cpu")
model_path = os.path.join(repo_root, "PBCNet2.pth")
model = torch.load(model_path, map_location=device, weights_only=False)
model.to(device)

test_dataset = LeadOptDataset(predict_csv_path)
test_dataloader = DataLoader(test_dataset, collate_fn=collate_fn, batch_size=8, shuffle=False)

_, _, _, _, valid_prediction, _, _, _, _ = predict(model, test_dataloader, device)

df = pd.read_csv(predict_csv_path)
df["pred_delta_pIC50"] = valid_prediction  # lig1 기준 - lig2 기준 상대 결합친화도(pIC50) 예측값
df.to_csv(predict_csv_path, index=False)

print("\n=== 예측 결과 (pred_delta_pIC50 = lig1의 pIC50 - lig2의 pIC50) ===")
print(df[["lig1", "lig2", "pred_delta_pIC50"]].to_string(index=False))
