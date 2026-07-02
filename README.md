# GC-MS 自动鉴定分析工具

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

## 验证与定位

**独立、可复现、结构级验证**（`eval/bench_massbank.py`）：把 MassBank **独立采集**的
EI 谱送入本项目的 NIST 搜索，按 **CAS / InChIKey 结构**判定命中（不是名字比对——杜绝命名
歧义、不循环）。在类挥发物子集（未衍生化、MW≤200，n=300）上：

- **只要目标化合物在 NIST 库中可检索到（占 71%），引擎就以 top-1 ≈ 89%、top-5 ≈ 100%
  把它排到最前** —— 排序能力已逼近 EI 质谱的可分辨极限。
- 整体 top-1 ≈ 60%、top-5 ≈ 70%；剩下约 29% 的差距来自**库未收录**或 **EI 物理上无法
  区分的同分异构体**（数据/物理极限，**非算法可解**）。

而这只是**质谱单证据**的成绩——真实运行时本工具还叠加**保留指数（RI）双证据核对**，
正好破解上面 MS 单独分不开的异构体，因此实际场景的可靠性高于此基准。这也说明：早期
"与人工 top-1 仅 41% 一致"主要是**人工端**（疲劳、只选 top-1）拉低，而非引擎问题。

> 严格起见：以上为**单一化合物库谱**（比真实色谱峰干净，属乐观上限）且**未用 RI**。
> 想要贴合自己样品的绝对准确率，用成分已知的标准品跑 `eval/truth_check.py`。
> 历史一致度数据见《项目技术档案.md》，回归守卫见 `eval/regression.py`。

**定位**：面向审核提效的**辅助**工具，不替代人工鉴定与标准品确认——核心价值是把
"逐个翻几百个峰"压到"只看几十个待复核峰"，且每个鉴定都带 MS+RI 双证据分级。

## 系统要求

- Python 3.11+
- NIST 质谱库
- Windows 操作系统

## 本地部署与复现

本项目为作者自用而建，`config.py` / `nist_engine.py` / `step1_parse.py` 等文件中的
库路径、DLL 路径为作者本机的绝对路径。在其他环境部署时，需要：

1. 安装依赖：`pip install -r requirements.txt`（Python 3.11+，Windows）
2. 自备 NIST 质谱库（版权原因不随仓库分发），并将上述文件中的 mainlib / replib、
   ThermoRawFileParser DLL 等路径改为本机的实际路径
3. 原生读取 Thermo `.RAW` 需 `thermo_lib`（RawFileReader.dll，同样不随仓库分发）

如需在本地部署或复现本项目，欢迎来信交流：**g0ho@qq.com**

## 许可证

MIT License · Copyright (c) 2026 go ho
