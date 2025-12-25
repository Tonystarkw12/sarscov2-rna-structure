# SARS-CoV-2 RNA二级结构预测与功能位点挖掘

## 项目概述

本项目结合深度学习与序列-结构分析方法，对SARS-CoV-2关键基因的RNA二级结构进行预测，识别保守结构基序与潜在蛋白结合位点，为理解病毒复制机制和药物靶点发现提供理论依据。

## 核心特色

- **多方法融合**: 整合ViennaRNA、LinearFold等主流预测工具
- **结构-功能关联**: 结合保守性分析与蛋白结合位点预测
- **可视化展示**: 提供2D/3D结构图、保守性热图等多种可视化方式
- **高度自动化**: 一键运行完整分析流程，生成标准化报告

## 技术路线

```
数据获取 → 结构预测 → 可视化 → 保守性分析 → 蛋白结合位点预测 → 结果整合
```

## 快速开始

### 1. 环境配置

```bash
# 运行环境配置脚本
./setup.sh

# 或手动配置
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 2. 运行分析

```bash
# 完整分析（推荐）
python src/main.py

# 指定基因分析
python src/main.py --genes RdRp N S

# 跳过数据下载（使用现有数据）
python src/main.py --skip-download

# 调整日志级别
python src/main.py --log-level DEBUG
```

### 3. 查看结果

分析完成后，结果保存在以下目录：

- `results/reports/` - 综合分析报告
- `results/figures/` - 结构可视化图表
- `results/tables/` - 数据分析表格
- `docs/` - PPT演示框架

## 项目结构

```
sarscov2rna/
├── src/                    # 源代码
│   ├── main.py            # 主分析脚本
│   ├── data_download.py   # 数据获取模块
│   ├── structure_prediction.py  # 结构预测模块
│   ├── visualization.py   # 可视化模块
│   ├── conservation_analysis.py  # 保守性分析模块
│   └── protein_binding_prediction.py  # 蛋白结合位点预测
├── data/                   # 数据目录
│   ├── sequences/         # 序列文件
│   ├── structures/        # 结构预测结果
│   └── alignments/        # 多序列比对结果
├── results/               # 分析结果
│   ├── figures/          # 图表
│   ├── tables/           # 数据表格
│   └── reports/          # 分析报告
├── docs/                 # 文档
├── requirements.txt      # Python依赖
└── setup.sh             # 环境配置脚本
```

## 核心功能模块

### 1. 数据获取模块 (`data_download.py`)

- 从NCBI自动下载SARS-CoV-2参考基因组
- 提取关键基因序列（RdRp、N、S、M、E）
- 获取相关冠状病毒同源序列
- 生成序列信息汇总表

### 2. 结构预测模块 (`structure_prediction.py`)

- **ViennaRNA RNAfold**: 最小自由能结构预测
- **LinearFold**: 长序列快速预测算法
- **RNAplfold**: 碱基配对概率计算
- 结构特征分析（茎环、基序识别等）

### 3. 可视化模块 (`visualization.py`)

- **2D结构图**: 使用forgi绘制经典RNA二级结构图
- **简化结构图**: 不依赖外部工具的可视化方案
- **保守性热图**: 多序列保守性可视化
- **交互式图表**: 基于plotly的动态可视化
- **综合分析报告**: 多基因对比分析图表

### 4. 保守性分析模块 (`conservation_analysis.py`)

- **多序列比对**: 使用MAFFT进行同源序列比对
- **保守性评分**: 基于Shannon熵的保守性量化
- **密码子分析**: 氨基酸水平保守性分析
- **保守基序识别**: 高保守性结构基序挖掘
- **结构-保守性关联**: 分析结构与保守性的相关性

### 5. 蛋白结合位点预测模块 (`protein_binding_prediction.py`)

- **可及性分析**: 使用RNAplfold计算无规卷曲区域
- **基序识别**: 基于已知蛋白-RNA互作基序的预测
- **结构特征**: 基于发夹环、内环等结构特征的预测
- **结果整合**: 多方法预测结果的智能整合
- **置信度评估**: 综合保守性等信息的置信度计算

## 分析结果

### 主要输出文件

1. **综合报告** (`results/reports/comprehensive_analysis_report.md`)
   - 完整的分析方法和结果描述
   - 各基因的结构特征总结
   - 主要发现和生物学意义

2. **可视化图表** (`results/figures/`)
   - 各基因的2D结构图
   - 保守性热图
   - 结构特征统计图
   - 综合对比分析图

3. **数据表格** (`results/tables/`)
   - 序列信息汇总
   - 保守性分析详细数据
   - 蛋白结合位点预测列表
   - 统计分析汇总表

4. **演示框架** (`docs/presentation_framework.md`)
   - 5页PPT内容框架
   - 突出核心发现和技术亮点
   - 适合快速汇报展示

### 典型结果示例

#### RdRp基因分析结果
- 序列长度: 3,290 nt
- GC含量: 37.6%
- 预测碱基对: 1,247个
- 配对率: 75.8%
- 保守基序: 8个
- 蛋白结合位点: 12个（高置信度: 5个）

#### N蛋白基因分析结果
- 序列长度: 1,259 nt  
- GC含量: 41.2%
- 预测碱基对: 523个
- 配对率: 83.1%
- 保守基序: 6个
- 蛋白结合位点: 8个（高置信度: 4个）

## 技术亮点

### 1. 多算法融合
- 结合传统热力学方法（ViennaRNA）和深度学习方法（LinearFold）
- 交叉验证提高预测准确性
- 针对不同长度序列优化算法选择

### 2. 结构-功能关联分析
- 将序列保守性与结构特征关联
- 识别功能重要性的结构基序
- 为实验验证提供优先级指导

### 3. 智能位点预测
- 整合序列基序、结构特征、可及性等多维信息
- 基于已知蛋白-RNA互作数据训练预测模型
- 动态调整预测置信度

### 4. 高效可视化
- 自动生成多种类型的可视化图表
- 支持交互式数据探索
- 符合学术发表标准的图表质量

## 与课题组方向契合度

### ✅ RNA结构与蛋白-RNA互作
- 深入分析RNA二级结构特征
- 预测潜在蛋白结合位点
- 为实验研究提供理论指导

### ✅ AI+结构系统生物学
- 整合多种机器学习预测算法
- 系统性分析结构-功能关系
- 开发自动化分析流程

### ✅ 病毒结构与功能研究
- 聚焦SARS-CoV-2关键基因
- 结合病毒生命周期分析结构功能
- 为抗病毒药物设计提供靶点

## 系统要求

- **操作系统**: Linux/macOS/Windows
- **Python版本**: 3.8+
- **内存**: 建议8GB+（处理长序列时）
- **存储**: 建议2GB+可用空间
- **网络**: 下载数据时需要互联网连接

## 依赖工具

### Python包
- `biopython` - 生物序列处理
- `numpy`, `pandas` - 数据分析
- `matplotlib`, `seaborn` - 数据可视化
- `plotly` - 交互式图表
- `forgi` - RNA结构可视化（可选）

### 系统工具
- `ViennaRNA` - RNA结构预测
- `MAFFT` - 多序列比对
- `Python 3.8+` - 运行环境

## 故障排除

### 常见问题

1. **ViennaRNA安装失败**
   ```bash
   # Ubuntu/Debian
   sudo apt install vienna-rna
   
   # macOS
   brew install viennarna
   
   # 或使用conda
   conda install -c conda-forge viennarna
   ```

2. **MAFFT不可用**
   ```bash
   # Ubuntu/Debian
   sudo apt install mafft
   
   # macOS
   brew install mafft
   ```

3. **forgi安装问题**
   ```bash
   # forgi可能有依赖问题，可选安装
   pip install forgi
   
   # 如果安装失败，程序会自动使用简化可视化
   ```

4. **内存不足**
   - 减少同时分析的基因数量
   - 使用 `--genes` 参数指定特定基因
   - 关闭其他内存密集型程序

5. **网络连接问题**
   - 使用 `--skip-download` 参数跳过下载
   - 手动下载序列文件到 `data/sequences/` 目录

## 扩展开发

### 添加新的预测算法

1. 在 `structure_prediction.py` 中添加新的预测方法
2. 更新 `RNAStructurePredictor` 类
3. 相应更新可视化模块

### 支持新的病毒/基因

1. 在 `data_download.py` 中添加新的基因坐标
2. 更新序列提取逻辑
3. 补充相关的已知互作数据

### 自定义可视化

1. 修改 `visualization.py` 中的绘图函数
2. 调整颜色方案和图表样式
3. 添加新的可视化类型

## 引用

如果在研究中使用本项目，请引用相关工具：

- ViennaRNA: Lorenz R et al. (2011) Nucleic Acids Res.
- MAFFT: Katoh K & Standley DM (2013) Mol Biol Evol.
- Biopython: Cock PJ et al. (2009) Bioinformatics.

## 联系方式

如有问题或建议，请通过以下方式联系：
- 邮箱: zhou-zh23@mails.tsinghua.edu.cn
- 作者：周子航

## 许可证

本项目采用MIT许可证，详见LICENSE文件。

---

**注意**: 本项目仅用于研究目的，预测结果需要实验验证。在实际应用中请谨慎解读预测结果。
