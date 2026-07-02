# GC-MS 自动数据处理系统 v2.1

**原始数据 → 峰检测 → NIST 搜索 → RI 标定 → Excel 审核报告，一键完成。**

---

## 作者

go ho

## 简介

本系统用于 GC-MS（气相色谱-质谱联用）数据的全自动处理。**跨仪器支持**：Thermo Fisher `.RAW`（原生直读）、`.mzML`、岛津 GCMSsolution `.qgd`，自动识别格式。自动完成基线校正、色谱峰检测、空气背景离子剔除、NIST 质谱库搜索、保留指数（RI）标定与核对，最终输出分级的 Excel 审核报告。

核心技术：SNIP 基线校正（向量化实现）、ICIS 峰检测（赛默飞 Genesis 同等逻辑）、空气背景离子剔除（解决低丰度峰被水峰主导的问题）、通过 pyms-nist-search 调用 NIST 官方引擎（mainlib + 复本库 replib）、**按样品 RI 标定**（用本仪器跑的正构烷烃标品自动建立 RT→RI，并对每个鉴定做 MS+RI 双证据核对）。

保留指数（RI）数据库含约 53,000 条 DB-5 与 16,600 条 DB-WAX 文献值。

## 快速开始

```bash
# 安装依赖（版本已锁定；olefile 读 .qgd，pythonnet 原生读 .RAW）
pip install -r requirements.txt

# 图形界面（推荐）：选样品(.RAW/.mzML/.qgd) + RI 标品 + 溶剂延迟 → START
双击 软件\启动.bat

# 命令行
python pipeline.py --files sample.qgd --standard alkanes.qgd --nist --solvent-delay 2.0
```

## 输出

每个样品输出一个独立的 Excel（文件名=样品名，多样品互不合并）。每个 Excel 含三张表（按此顺序看最省力）：
- **汇总**：峰总数、可信免审、待复核、微量峰、污染物、低置信的数量
- **待复核**：只列需人工核对的峰（黄色无 RI 佐证 + RI 可疑），按面积排序；
  可信(绿+RI双证据)/微量/污染物/低置信已自动归类，无需逐个翻
- **Results**：全部峰的前 8 名候选，含 FMF/RMF、双柱 RI、`RI_Check`
  （OK/suspect/noRI/noCal，MS+RI 双证据核对）、`Source`（replib 标记）、
  审核状态（绿/黄/红/灰）、人工确认列

## 性能（跨仪器验证）

- **Thermo**（海带挥发物，对照 49 个人工鉴定，`eval/score.py`）：top-1 41%、top-5 47%、top-8 53%。
- **岛津 .qgd**（对照 GCMSsolution 结果，61 峰）：top-1 50%、top-5 72%。
- 这是"与厂商/人工结果的一致度"，非绝对准确率（双方都是 NIST 搜索结果）；
  剩余分歧主要是 EI 物理极限（支链烷烃同分异构体）+ 弱峰 + 厂商低置信标注。
  详见《项目技术档案.md》，回归守卫见 `eval/regression.py`。
- 想要**绝对准确率**：跑一个成分已知的标准品，用 `eval/truth_check.py`
  （`--candidates <样品>_candidates.json --truth <清单>`）测 top-1/5/8 真实命中率。
  真值独立于人工与管道，是唯一真正的准确率数字。

## 系统要求

- Python 3.11+
- NIST 质谱库
- Windows 操作系统

## 在自己电脑上复现

本项目为作者自用而建，`config.py` / `nist_engine.py` / `step1_parse.py` 等文件中的
库路径、DLL 路径是作者本机的绝对路径。若要在自己的电脑上复现，需要：

1. 安装依赖：`pip install -r requirements.txt`（Python 3.11+，Windows）
2. 自备 NIST 质谱库（版权原因不随仓库分发），并把上述文件里的 mainlib / replib、
   ThermoRawFileParser DLL 等路径改成你本机的实际路径
3. 原生读 Thermo `.RAW` 需 `thermo_lib`（RawFileReader.dll，同样不随仓库分发）

如果你真的想在自己电脑上复现这个项目，遇到问题欢迎联系：**g0ho@qq.com**

## 许可证

MIT License · Copyright (c) 2026 go ho
