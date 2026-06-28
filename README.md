# GC-MS 自动数据处理系统 v2.1

**Thermo .RAW → 峰检测 → NIST 搜索 → Excel 报告，一键完成。**

---

## 作者

go ho

## 简介

本系统用于 GC-MS（气相色谱-质谱联用）数据的全自动处理。输入 Thermo Fisher 仪器的 `.RAW` 原始文件或 `.mzML` 开放格式文件，自动完成基线校正、色谱峰检测、谱图增强、NIST 质谱库搜索和保留指数（RI）计算，最终输出带颜色标注的 Excel 审核报告。

核心技术：SNIP 基线校正（向量化实现）、ICIS 峰检测算法（赛默飞 Genesis 同等逻辑）、Xcalibur 风格的谱图增强（峰顶多扫描平均减去背景扣除），以及通过 pyms-nist-search 调用的 NIST 官方质谱搜索引擎。

保留指数（RI）数据库包含 46,000 条 DB-5 文献值和 16,000 条 DB-WAX 文献值，支持双柱交叉验证鉴定。

## 性能

在本课题样品（海带挥发性化合物，DB-WAX 柱）上，系统自动鉴定与人工鉴定的符合度约为三分之一。剩余三分之二的分歧源自：（1）支链烷烃同分异构体在 EI 质谱中无法区分且缺少 DB-WAX 文献 RI 值；（2）E/Z 异构体的四极杆质谱物理极限；（3）低信号化合物的信噪比不足。详细性能说明见《使用教程.txt》。

## 快速开始

```bash
# 安装依赖
pip install numpy scipy pandas openpyxl pymzml pyms-nist-search

# 图形界面（推荐）
双击 软件\启动.bat

# 命令行
python pipeline.py --mzml-files sample.mzML --nist
```

## 输出

Excel 报告包含：
- 峰面积与峰高
- 每个峰的前三名 NIST 候选鉴定
- 匹配因子（FMF/RMF）
- DB-WAX 和 DB-5 双柱文献保留指数
- 自动审核状态：绿（自动确认）· 黄（待审核）· 红（污染物）· 灰（低置信）
- 人工确认列（供用户填写）

## 系统要求

- Python 3.11+
- NIST 质谱库
- Windows 操作系统

## 许可证

MIT License · Copyright (c) 2026 go ho
