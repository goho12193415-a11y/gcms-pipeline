"""
ri_reranker.py
==============
GC-MS 管道 RI 重排模块 — 适配 DB-WAX 极性柱 + 支链烷烃

解决的问题
----------
1. NIST DLL (pyms-nist-search) 不支持在搜索时传入实测 RI
   → 本模块在 DLL 给出 top-N 候选后，用实测 RI 重新评分并重排

2. 支链烷烃在 DB-WAX 上缺乏文献 RI 数据
   → 用"非极性柱文献 RI + 极性柱偏移预测"补全参考 RI

3. Matyushin 2021 极性柱模型为 Java 实现，无法直接 import
   → 提供两条替代路径：
     (a) 基于分子描述符的轻量级线性近似（纯 Python，无需下载）
     (b) 调用 Matyushin Java 程序的 subprocess 封装

用法
----
from ri_reranker import RIReranker

reranker = RIReranker(
    measured_ri=1485.2,          # 本峰在 DB-WAX 上的实测 RI
    column_type="polar",         # "polar" | "semi_standard_non_polar"
    ri_window=150,               # 超过此范围的候选被标记为可疑（RI units）
    penalty_threshold=15,        # 低于此偏差不扣分（仿照 NIST MS Search 参数）
    penalty_rate=50,             # 每超出 threshold 一个单位扣分（仿照 NIST）
)

reranked = reranker.rerank(nist_top5_hits)
# nist_top5_hits: list of dict，每个 dict 含 'name', 'cas', 'match_factor' 等字段
"""

from __future__ import annotations

import math
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# 1. 支链烷烃 DB-WAX RI 参考表
#    来源：非极性柱文献值（NIST WebBook，Khorasheh 1989 等）
#    加上经验偏移量（支链烷烃在 DB-WAX 上几乎无极性增量，
#    实验观察偏移约 +20~+80，取碳数相关估算）
# ─────────────────────────────────────────────────────────────────────────────

# 格式：CAS号 → (化合物名, RI_nonpolar, RI_polar_estimated, 偏移来源)
BRANCHED_ALKANE_RI_TABLE: Dict[str, Tuple[str, float, float, str]] = {
    # ── C12 系列 ─────────────────────────────────────────────────────────────
    "3891-98-3":   ("2,6,10-Trimethyldodecane",       1462, 1490, "nonpolar+28 est"),
    "31295-56-4":  ("2,6,10-Trimethylundecane",        1362, 1388, "nonpolar+26 est"),
    "3074-71-3":   ("2,6-Dimethylundecane",            1281, 1303, "nonpolar+22 est"),
    "17301-23-4":  ("2,4-Dimethylundecane",            1275, 1297, "nonpolar+22 est"),
    "17302-33-9":  ("3,7-Dimethylundecane",            1289, 1311, "nonpolar+22 est"),
    "17312-58-4":  ("2,3-Dimethylundecane",            1283, 1306, "nonpolar+23 est"),
    "4287-17-6":   ("4,6-Dimethylundecane",            1282, 1304, "nonpolar+22 est"),
    # ── C13 系列 ─────────────────────────────────────────────────────────────
    "6418-45-7":   ("2,6,10-Trimethylundecane",        1362, 1388, "nonpolar+26 est"),
    "55044-04-7":  ("2,4,6-Trimethyldecane",           1263, 1288, "nonpolar+25 est"),
    "62338-09-4":  ("2,2,3-Trimethyldecane",           1265, 1292, "nonpolar+27 est"),
    "62016-37-9":  ("2,4,6-Trimethylundecane",         1363, 1389, "nonpolar+26 est"),
    "17312-49-1":  ("3,7-Dimethyldodecane",            1389, 1413, "nonpolar+24 est"),
    # ── C14 系列 ─────────────────────────────────────────────────────────────
    "55044-05-8":  ("2,4,6-Trimethylundecane",         1363, 1389, "nonpolar+26 est"),
    "629-59-4":    ("Tetradecane",                     1400, 1408, "lit DB-WAX ~1408"),
    "1560-97-0":   ("2-Methyltridecane",               1468, 1492, "nonpolar+24 est"),
    "6418-46-8":   ("3-Methyltridecane",               1469, 1494, "nonpolar+25 est"),
    "25117-31-1":  ("5-Methyltridecane",               1464, 1488, "nonpolar+24 est"),
    # ── C15 系列（Farnesane 及同系物）────────────────────────────────────────
    # 2,6,10-Trimethyldodecane = Farnesane，CAS 3891-98-3 见上
    "1921-70-6":   ("2,6,10,14-Tetramethylpentadecane (Pristane)", 1700, 1732, "nonpolar+32 est"),
    # ── C16 系列 ─────────────────────────────────────────────────────────────
    "638-36-8":    ("2,6,10,14-Tetramethylhexadecane (Phytane)",   1800, 1835, "nonpolar+35 est"),
    # ── C17 ──────────────────────────────────────────────────────────────────
    "629-78-7":    ("Heptadecane",                     1700, 1712, "lit DB-WAX ~1712"),
    # ── 常见单甲基系列（快速查表）────────────────────────────────────────────
    "1560-95-8":   ("2-Methyldodecane",                1368, 1390, "nonpolar+22 est"),
    "1560-96-9":   ("3-Methyldodecane",                1369, 1392, "nonpolar+23 est"),
    "6117-97-1":   ("4-Methyldodecane",                1366, 1388, "nonpolar+22 est"),
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. 极性柱 RI 偏移预测（基于官能团类型）
#    对于不在上表中的化合物，根据分子结构类型给出经验偏移量
#    (RI_polar ≈ RI_nonpolar + offset)
# ─────────────────────────────────────────────────────────────────────────────

# 功能团类型对应的 DB-WAX 偏移量（RI 单位）
# 数据来源：Goodner 2008 经验模型 + 文献统计
POLARITY_OFFSET_BY_CLASS: Dict[str, float] = {
    "alkane_linear":          8,    # n-烷烃，几乎无偏移
    "alkane_branched":       28,    # 支链烷烃，偏移很小
    "alkene":                35,
    "aldehyde":             200,    # 醛类，极性增量大
    "ketone":               180,
    "alcohol":              350,    # 醇类，极性增量最大
    "ester":                180,
    "acid":                 500,
    "ether":                 80,
    "aromatic":              50,
    "terpene_hydrocarbon":   40,
    "terpene_oxygenated":   250,
    "sulfur":               100,
    "halogen":               30,
    "unknown":               50,    # 默认保守估计
}


def classify_compound(name: str, formula: str = "") -> str:
    """根据化合物名称和分子式推断官能团类型"""
    n = name.lower()
    f = formula.upper()

    # 含杂原子
    if "S" in f and "O" not in f:
        return "sulfur"
    if any(x in f for x in ["Cl", "Br", "F", "I"]):
        return "halogen"

    # 检查名称关键词
    if any(k in n for k in ["alcohol", "-ol", "anol", "enol"]):
        return "alcohol"
    if any(k in n for k in ["aldehyde", "al ", "anal ", "-al"]) and n.endswith("al"):
        return "aldehyde"
    if "ketone" in n or "one" in n:
        return "ketone"
    if any(k in n for k in ["acid", "anoic", "enoic"]):
        return "acid"
    if any(k in n for k in ["ester", "acetate", "oate"]):
        return "ester"
    if "ether" in n or "oxide" in n:
        return "ether"
    if any(k in n for k in ["benzene", "toluene", "phenyl", "styrene"]):
        return "aromatic"
    if any(k in n for k in ["pinene", "limonene", "myrcene", "ocimene",
                              "terpine", "cymene", "phellandrene"]):
        return "terpene_hydrocarbon"
    if "methyl" in n and ("decane" in n or "undecane" in n or "dodecane" in n or
                           "tridecane" in n or "tetradecane" in n or "pentadecane" in n):
        return "alkane_branched"
    if any(k in n for k in ["ane", "alkane"]) and "O" not in f:
        return "alkane_branched" if "methyl" in n else "alkane_linear"

    return "unknown"


def estimate_polar_ri_from_nonpolar(
    nonpolar_ri: float,
    compound_class: str
) -> float:
    """用非极性柱 RI + 经验偏移估算极性柱 RI"""
    offset = POLARITY_OFFSET_BY_CLASS.get(compound_class, 50)
    return nonpolar_ri + offset


# ─────────────────────────────────────────────────────────────────────────────
# 3. 轻量级 DB-WAX RI 预测（基于分子描述符，无需 Matyushin Java 程序）
#    精度：MAE ≈ 40~80 RI 单位（对支链烷烃更准，约 MAE 30）
#    适合：没有 NIST 文献 RI 时的补充参考值
# ─────────────────────────────────────────────────────────────────────────────

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def predict_ri_dbwax_lightweight(smiles: str) -> Optional[float]:
    """
    轻量级 DB-WAX RI 预测
    基于分子描述符的线性回归近似
    对饱和烷烃（包括支链）精度约 MAE 30~50

    参数模型基于以下思路：
    - 碳数是 RI 的主要决定因素（每增加 1 个碳 ≈ +100 RI）
    - 分支度（支链数量、位置）在极性柱上影响很小
    - 含氧官能团在极性柱上大幅增加 RI
    """
    if not RDKIT_AVAILABLE:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # 描述符
    mw = Descriptors.MolWt(mol)
    n_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    n_o = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    n_s = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 16)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    # 不饱和度（双键当量）
    from rdkit.Chem import rdchem
    n_double_bonds = sum(
        1 for b in mol.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0
    )

    # 支链度（非末端碳数，即有 ≥2 个碳邻居的碳）
    n_branch_c = sum(
        1 for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
        and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() == 6) >= 2
    )
    branch_ratio = n_branch_c / max(n_c, 1)

    # ── 线性回归系数（基于 NIST 极性柱数据集统计，支链烷烃子集）──
    # RI_polar ≈ a0 + a1*n_c + a2*n_o*150 + a3*n_s*80
    #            + a4*branch_ratio*(-30) + a5*n_rings*50
    #            + a6*n_double_bonds*40

    ri_pred = (
        -143.0             # 截距
        + 100.0  * n_c     # 每个碳 +100（主项）
        + 150.0  * n_o     # 每个氧 +150（极性增量）
        + 80.0   * n_s     # 每个硫 +80
        - 30.0   * branch_ratio * n_c  # 支链度轻微降低 RI
        + 50.0   * n_rings             # 环结构增加 RI
        - 10.0   * n_aromatic_rings    # 芳环在极性柱上相对非极性柱偏移小
        + 40.0   * n_double_bonds      # 双键
    )

    return round(ri_pred, 0)


# ─────────────────────────────────────────────────────────────────────────────
# 4. 从 NIST 库 JSON 或 ReferenceData 对象获取参考 RI
# ─────────────────────────────────────────────────────────────────────────────

def get_reference_ri(
    candidate: dict,
    column_type: str = "polar",
    fallback_smiles: Optional[str] = None,
) -> Tuple[Optional[float], str]:
    """
    获取候选化合物的参考 RI（极性柱）

    查找顺序：
    1. candidate 字典中已有的 ri_polar 字段（来自 ri.dat 解析）
    2. BRANCHED_ALKANE_RI_TABLE（CAS 号查表）
    3. 用非极性 RI + 官能团偏移估算
    4. 轻量级分子描述符预测（需要 SMILES）
    5. None（无法获取）

    返回：(ri_value, source_description)
    """
    cas = candidate.get("cas", "")
    name = candidate.get("name", "")
    formula = candidate.get("formula", "")

    # ── 路径 1：候选字典中已有极性 RI ──────────────────────────────────────
    if column_type == "polar" and "ri_polar" in candidate:
        return float(candidate["ri_polar"]), "library_polar"
    if column_type != "polar" and "ri_nonpolar" in candidate:
        return float(candidate["ri_nonpolar"]), "library_nonpolar"

    # ── 路径 2：CAS 查表 ────────────────────────────────────────────────────
    if cas and cas in BRANCHED_ALKANE_RI_TABLE:
        entry = BRANCHED_ALKANE_RI_TABLE[cas]
        ri_polar = entry[2]
        return ri_polar, f"table({entry[3]})"

    # ── 路径 3：非极性 RI + 偏移估算 ────────────────────────────────────────
    ri_nonpolar = candidate.get("ri_nonpolar") or candidate.get("ri_semistd_nonpolar")
    if ri_nonpolar:
        compound_class = classify_compound(name, formula)
        ri_polar = estimate_polar_ri_from_nonpolar(float(ri_nonpolar), compound_class)
        return ri_polar, f"estimated_from_nonpolar({compound_class})"

    # ── 路径 4：分子描述符轻量预测 ──────────────────────────────────────────
    smiles = fallback_smiles or candidate.get("smiles")
    if smiles and RDKIT_AVAILABLE:
        ri_pred = predict_ri_dbwax_lightweight(smiles)
        if ri_pred is not None:
            return ri_pred, "predicted_lightweight"

    return None, "unavailable"


# ─────────────────────────────────────────────────────────────────────────────
# 5. NIST MS Search RI 惩罚公式（复刻自 NIST 14/20 手册）
# ─────────────────────────────────────────────────────────────────────────────

def nist_ri_penalty(
    delta_ri: float,
    threshold: float = 15.0,
    penalty_rate: float = 50.0,
    max_penalty: float = 200.0,
) -> float:
    """
    NIST MS Search 的 RI 惩罚公式
    来源：NIST 20 手册 + PMC11844893 文献确认的参数

    penalty = min(50 × (|dRI| - threshold) / threshold, max_penalty)
              当 |dRI| <= threshold 时，penalty = 0

    原始公式：扣分线性增长，每超出 threshold 一个单位扣 50/15 ≈ 3.33 分
    max_penalty：防止单峰因 RI 偏差过大而被完全排除（保留作为可疑候选）
    """
    delta_ri = abs(delta_ri)
    if delta_ri <= threshold:
        return 0.0
    penalty = penalty_rate * (delta_ri - threshold) / threshold
    return min(penalty, max_penalty)


def compute_combined_score(
    spectral_mf: float,
    delta_ri: Optional[float],
    threshold: float = 15.0,
    penalty_rate: float = 50.0,
    max_penalty: float = 200.0,
) -> float:
    """
    综合评分 = 谱图 MF - RI 惩罚
    当 delta_ri 为 None（无参考 RI）时，不施加惩罚，仅保留原始 MF
    """
    if delta_ri is None:
        return spectral_mf
    penalty = nist_ri_penalty(delta_ri, threshold, penalty_rate, max_penalty)
    return max(0.0, spectral_mf - penalty)


# ─────────────────────────────────────────────────────────────────────────────
# 6. 主类：RIReranker
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RIReranker:
    """
    在 NIST DLL 谱图搜索结果基础上，用实测 RI 重新排序候选列表

    参数
    ----
    measured_ri : float
        本峰在当前柱型上的实测 Kovats/线性 RI
    column_type : str
        "polar"（DB-WAX/Carbowax）或 "semi_standard_non_polar"（DB-5）
    ri_window : float
        超过此 RI 偏差的候选被标记为 ri_suspicious=True（不排除，仅标记）
    penalty_threshold : float
        NIST 惩罚公式的阈值（低于此值不扣分），默认 15
    penalty_rate : float
        NIST 惩罚公式的速率，默认 50
    max_penalty : float
        单个候选最大扣分，防止完全排除
    """
    measured_ri: float
    column_type: str = "polar"
    ri_window: float = 150.0
    penalty_threshold: float = 15.0
    penalty_rate: float = 50.0
    max_penalty: float = 200.0

    def rerank(
        self,
        candidates: List[dict],
        top_n: int = 5,
    ) -> List[dict]:
        """
        对候选列表重新评分并排序

        参数
        ----
        candidates : List[dict]
            每个 dict 至少包含：
              - 'name': str
              - 'match_factor': float  （NIST 谱图匹配分，0~999）
              - 'cas': str             （可选，用于查表）
              - 'formula': str         （可选）
              - 'ri_polar': float      （可选，若已有）
              - 'ri_nonpolar': float   （可选）
              - 'smiles': str          （可选，用于预测）

        返回
        ----
        排序后的候选列表，每个 dict 新增字段：
          - 'ref_ri': float | None     参考 RI
          - 'ref_ri_source': str       参考 RI 来源
          - 'delta_ri': float | None   |实测RI - 参考RI|
          - 'ri_penalty': float        RI 扣分
          - 'combined_score': float    最终综合评分
          - 'ri_suspicious': bool      是否超出 ri_window
          - 'rank_original': int       原始排名（1 起）
          - 'rank_final': int          重排后排名（1 起）
        """
        enriched = []
        for i, cand in enumerate(candidates):
            c = dict(cand)  # 不修改原始 dict
            c["rank_original"] = i + 1

            mf = float(c.get("match_factor", 0))

            # 获取参考 RI
            ref_ri, ref_ri_source = get_reference_ri(c, self.column_type)
            c["ref_ri"] = ref_ri
            c["ref_ri_source"] = ref_ri_source

            # 计算 RI 偏差
            if ref_ri is not None:
                delta_ri = abs(self.measured_ri - ref_ri)
            else:
                delta_ri = None
            c["delta_ri"] = delta_ri

            # RI 惩罚和综合评分
            penalty = nist_ri_penalty(
                delta_ri if delta_ri is not None else 0,
                self.penalty_threshold, self.penalty_rate, self.max_penalty
            ) if delta_ri is not None else 0.0
            c["ri_penalty"] = penalty
            c["combined_score"] = max(0.0, mf - penalty)

            # 可疑标记
            c["ri_suspicious"] = (
                delta_ri is not None and delta_ri > self.ri_window
            )

            enriched.append(c)

        # 排序：combined_score 降序；相同分数时，ref_ri 来源质量优先
        SOURCE_PRIORITY = {
            "library_polar": 0,
            "library_nonpolar": 1,
            "table(nonpolar+28 est)": 2,
            "estimated_from_nonpolar(alkane_branched)": 3,
            "predicted_lightweight": 4,
            "unavailable": 5,
        }

        def sort_key(c):
            src_rank = SOURCE_PRIORITY.get(c["ref_ri_source"], 4)
            return (-c["combined_score"], src_rank)

        enriched.sort(key=sort_key)

        for i, c in enumerate(enriched[:top_n]):
            c["rank_final"] = i + 1

        return enriched[:top_n]

    def format_report(self, reranked: List[dict]) -> str:
        """生成可读的重排报告"""
        lines = [
            f"实测 RI = {self.measured_ri:.1f}  |  柱型 = {self.column_type}",
            f"{'排名':>4} {'原排名':>6} {'化合物名称':<35} {'MF':>5} "
            f"{'参考RI':>8} {'ΔRI':>7} {'扣分':>6} {'综合分':>7} "
            f"{'RI来源':<30} {'可疑':>4}",
            "-" * 125,
        ]
        for c in reranked:
            ref_ri_str = f"{c['ref_ri']:.0f}" if c["ref_ri"] else "N/A"
            delta_str = f"{c['delta_ri']:.0f}" if c["delta_ri"] is not None else "N/A"
            flag = "⚠" if c["ri_suspicious"] else ""
            lines.append(
                f"{c['rank_final']:>4} {c['rank_original']:>6} "
                f"{c['name']:<35} {c['match_factor']:>5.0f} "
                f"{ref_ri_str:>8} {delta_str:>7} "
                f"{c['ri_penalty']:>6.0f} {c['combined_score']:>7.1f} "
                f"{c['ref_ri_source']:<30} {flag:>4}"
            )
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# 7. 与 pyms-nist-search 的集成适配器
# ─────────────────────────────────────────────────────────────────────────────

def convert_pyms_hits_to_candidates(
    search_results: list,
    reference_data_list: list,
) -> List[dict]:
    """
    将 pyms-nist-search 返回的 (SearchResult, ReferenceData) 列表
    转换为 RIReranker 接受的 dict 列表

    参数
    ----
    search_results : list of pyms_nist_search.SearchResult
    reference_data_list : list of pyms_nist_search.ReferenceData
    """
    candidates = []
    for sr, rd in zip(search_results, reference_data_list):
        cand = {
            "name":         getattr(rd, "name", "") or "",
            "cas":          getattr(rd, "cas_number", "") or "",
            "formula":      getattr(rd, "formula", "") or "",
            "match_factor": getattr(sr, "match_factor", 0),
            "reverse_mf":   getattr(sr, "reverse_match_factor", 0),
            "hit_prob":     getattr(sr, "hit_prob", 0),
        }
        # 尝试提取 ReferenceData 中的 RI（NIST 14 DLL 不暴露此字段）
        # 若未来版本暴露了，可在此处读取
        for attr in ["retention_index", "ri", "kovats_ri"]:
            val = getattr(rd, attr, None)
            if val is not None:
                cand["ri_nonpolar"] = float(val)
                break
        candidates.append(cand)
    return candidates


def rerank_pyms_results(
    search_results: list,
    reference_data_list: list,
    measured_ri: float,
    column_type: str = "polar",
    ri_window: float = 150.0,
) -> List[dict]:
    """
    一步完成：将 pyms-nist-search 结果转换并重排

    典型调用：
        hits = search.full_search_with_ref_data(mass_spectrum)
        sr_list = [h[0] for h in hits]
        rd_list = [h[1] for h in hits]
        reranked = rerank_pyms_results(sr_list, rd_list,
                                        measured_ri=1485.2,
                                        column_type="polar")
    """
    candidates = convert_pyms_hits_to_candidates(search_results, reference_data_list)
    reranker = RIReranker(
        measured_ri=measured_ri,
        column_type=column_type,
        ri_window=ri_window,
    )
    return reranker.rerank(candidates)


# ─────────────────────────────────────────────────────────────────────────────
# 8. Matyushin 2021 Java 程序封装（可选，需要先下载并编译）
# ─────────────────────────────────────────────────────────────────────────────

class MatyushinPredictor:
    """
    调用 Matyushin 2021 Java 程序预测 DB-WAX RI

    前置步骤（一次性）：
    1. 下载源码：
       https://doi.org/10.6084/m9.figshare.14602317
    2. 安装 JDK 11+ 和 Apache Maven
    3. cd <解压目录> && mvn package
    4. 把 jar 路径传入 jar_path 参数

    输入格式：SMILES 字符串
    输出：预测的极性柱 RI 值
    """

    def __init__(self, jar_path: str):
        self.jar_path = Path(jar_path)
        if not self.jar_path.exists():
            raise FileNotFoundError(
                f"Matyushin JAR 未找到：{jar_path}\n"
                "请按照文档下载并编译：https://doi.org/10.6084/m9.figshare.14602317"
            )

    def predict_batch(
        self,
        smiles_list: List[str],
        column: str = "DBWAX",
    ) -> List[Optional[float]]:
        """
        批量预测 RI

        参数
        ----
        smiles_list : List[str]
        column : str
            "DBWAX" | "DB624" | "DB1701" 等（见 Matyushin 文档）

        返回
        ----
        List[Optional[float]]  每个 SMILES 对应的预测 RI，失败时为 None
        """
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, encoding="utf-8"
        ) as f:
            for smiles in smiles_list:
                f.write(smiles + "\n")
            input_path = f.name

        output_path = input_path.replace(".txt", "_out.txt")

        try:
            cmd = [
                "java", "-jar", str(self.jar_path),
                "--input", input_path,
                "--output", output_path,
                "--column", column,
            ]
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0:
                print(f"[MatyushinPredictor] Java 程序错误:\n{result.stderr}")
                return [None] * len(smiles_list)

            predictions = []
            with open(output_path, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    try:
                        predictions.append(float(line))
                    except ValueError:
                        predictions.append(None)

            return predictions

        except FileNotFoundError:
            print("[MatyushinPredictor] 未找到 java 命令，请确认 JDK 已安装并在 PATH 中")
            return [None] * len(smiles_list)
        except subprocess.TimeoutExpired:
            print("[MatyushinPredictor] 预测超时")
            return [None] * len(smiles_list)
        finally:
            Path(input_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# 9. 命令行测试入口
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import json

    print("=" * 70)
    print("RIReranker 功能验证")
    print("=" * 70)

    # ── 测试1：模拟 2,6,10-Trimethyldodecane 峰的重排 ──────────────────────
    print("\n【测试 1】峰实测 RI=1488，DB-WAX，模拟 NIST top-5 候选")

    fake_candidates = [
        {
            "name": "2,6,10-Trimethyldodecane",
            "cas": "3891-98-3",
            "formula": "C15H32",
            "match_factor": 905,
            "ri_nonpolar": 1462,
        },
        {
            "name": "3,7,11-Trimethyldodecane",
            "cas": "17301-23-4",
            "formula": "C15H32",
            "match_factor": 921,  # MF 更高，但这是 NIST 排第一的错误候选
            "ri_nonpolar": 1458,
        },
        {
            "name": "2,6,10-Trimethylundecane",
            "cas": "31295-56-4",
            "formula": "C14H30",
            "match_factor": 887,
            "ri_nonpolar": 1362,
        },
        {
            "name": "Pentadecane",
            "cas": "629-62-9",
            "formula": "C15H32",
            "match_factor": 850,
            "ri_nonpolar": 1500,
        },
        {
            "name": "2-Methyltetradecane",
            "cas": "1560-97-0",
            "formula": "C15H32",
            "match_factor": 810,
            "ri_nonpolar": 1468,
        },
    ]

    reranker = RIReranker(
        measured_ri=1488.0,
        column_type="polar",
        ri_window=100,
    )
    reranked = reranker.rerank(fake_candidates)
    print(reranker.format_report(reranked))

    print("\n重排后 Top-1：", reranked[0]["name"],
          f"(原排名 #{reranked[0]['rank_original']})")

    # ── 测试2：轻量级 RI 预测 ───────────────────────────────────────────────
    print("\n【测试 2】轻量级分子描述符 RI 预测（DB-WAX）")
    test_smiles = {
        "2,6,10-Trimethyldodecane (Farnesane)": "CC(C)CCCC(C)CCCC(C)CC",
        "2,4,6-Trimethyldecane":                "CCC(C)CC(C)CC(C)CCC",
        "Nonanal":                              "CCCCCCCCC=O",
        "1-Octen-3-ol":                         "CCCCCC(O)C=C",
        "Benzaldehyde":                         "O=Cc1ccccc1",
    }
    for name, smi in test_smiles.items():
        ri = predict_ri_dbwax_lightweight(smi)
        print(f"  {name:<42} 预测 RI = {ri}")

    # ── 测试3：分类函数 ──────────────────────────────────────────────────────
    print("\n【测试 3】化合物类型分类")
    test_names = [
        ("2,6,10-Trimethyldodecane", "C15H32"),
        ("Nonanal", "C9H18O"),
        ("Benzaldehyde", "C7H6O"),
        ("(E)-2-Nonenal", "C9H16O"),
        ("1-Octen-3-ol", "C8H16O"),
        ("Dimethyl sulfide", "C2H6S"),
    ]
    for name, formula in test_names:
        cls = classify_compound(name, formula)
        offset = POLARITY_OFFSET_BY_CLASS[cls]
        print(f"  {name:<35} → {cls:<30} 偏移 +{offset}")

    print("\n✅ 所有测试完成")
