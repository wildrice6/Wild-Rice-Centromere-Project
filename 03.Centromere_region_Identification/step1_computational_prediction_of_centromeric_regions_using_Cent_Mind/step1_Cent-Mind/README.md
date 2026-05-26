# 着丝粒区域预测 - Centromere Area Prediction

基于Transformer架构的染色体着丝粒区域预测深度学习模型。该模型使用多尺度k-mer特征，能够准确预测基因组序列中的着丝粒区域。

**✨ Includes pretrained models - Ready to use out of the box!**

## 项目简介

着丝粒（Centromere）是染色体的重要结构区域，在细胞分裂过程中起关键作用。本项目使用深度学习技术，基于序列特征自动识别和定位着丝粒区域。

### 主要特点

- **🎁 Pretrained Models**: Includes ready-to-use models, no training needed
- **Transformer架构**: 使用纯Transformer Encoder捕捉序列的长程依赖关系
- **多尺度特征**: 整合多个k-mer（64, 128, 256, 512）的统计特征
- **端到端训练**: 从原始特征直接预测着丝粒位置
- **高度可配置**: 灵活的模型配置和训练参数
- **完整工具链**: 包含训练、推理、评估和可视化等完整流程

### 模型架构

```
输入特征 (8维)
    ↓
线性投影 → d_model维
    ↓
位置编码
    ↓
Transformer Encoder (多层)
    ↓
    ├─→ 逐位置分类头 → 着丝粒概率
    └─→ 多尺度卷积 → 区间预测分数
```

## 环境要求

### 系统要求
- Python 3.8+
- CUDA 11.0+ (可选，用于GPU加速)

### 依赖包

```bash
pip install -r requirements.txt
```

主要依赖：
- PyTorch >= 1.10.0
- numpy >= 1.21.0
- pandas >= 1.3.0
- scikit-learn >= 1.0.0
- matplotlib >= 3.4.0
- tensorboard >= 2.8.0

## 快速开始

### Option 1: Inference from FASTA (Recommended for End Users)

If you have a genome FASTA file and a pretrained model:

```bash
# Install dependencies
pip install -r requirements.txt
conda install -c bioconda jellyfish

# One-command pipeline: FASTA → k-mer analysis → feature extraction → inference → BED output
chmod +x scripts/predict_from_fasta.sh
./scripts/predict_from_fasta.sh genome.fasta checkpoints/best_model.pt

# View results
cat predictions_output/predictions/predictions_summary.csv
```

**Detailed Guide**: 📖 [From FASTA to Predictions](docs/FROM_FASTA_TO_PREDICTION.md)

### Option 2: Inference from Feature CSV

If you already have feature CSV files:

```bash
cd src/training
python inference.py \
    --checkpoint checkpoints/best_model.pt \
    --input /path/to/features.csv \
    --output ./predictions \
    --threshold 0.3
```

**Output files:**
- `predictions.json`: Detailed predictions with probabilities
- `predictions_summary.csv`: Summary table
- `centromeres.bed`: BED format for genome browsers (IGV, UCSC)

### Option 3: Train Your Own Model

#### 1. Prepare Training Data

Input data should be CSV format with the following columns:
- `start`, `end`: Genomic coordinates
- `has_cen`: Label (0 or 1)
- Multi-scale k-mer features (8 columns for k=64,128,256,512)

File naming: `*_multi_k_summary.csv`

See [Data Format](docs/DATA_FORMAT.md) for details.

#### 2. Train Model

```bash
cd src/training
python train.py --data_dir /path/to/your/data --device cuda
```

Training parameters:
```bash
python train.py \
    --data_dir /path/to/data \
    --epochs 100 \
    --lr 5e-4 \
    --pos_weight 50.0 \
    --device cuda
```

#### 3. Monitor Training

```bash
tensorboard --logdir=training/logs
```

Open `http://localhost:6006` in your browser.

## 项目结构

```
centromere_prediction_github/
├── src/
│   ├── training/           # 训练模块
│   │   ├── config.py      # 配置文件
│   │   ├── model.py       # Transformer模型定义
│   │   ├── dataset.py     # 数据加载
│   │   ├── train.py       # 训练脚本
│   │   └── inference.py   # 推理脚本
│   └── evaluation/        # 评估模块
│       ├── evaluate_top5_prediction.py
│       ├── generate_summary_report.py
│       ├── compare_predictions.py
│       ├── calculate_average_metrics.py
│       └── process_coverage_statistics.py
├── examples/              # 示例脚本
│   └── run_training.sh   # 训练示例脚本
├── docs/                 # 文档
├── requirements.txt      # 依赖包列表
├── .gitignore           # Git忽略文件
└── README.md            # 本文件
```

## 模型配置

在 `src/training/config.py` 中可以调整以下参数：

### 模型参数
- `d_model`: Transformer特征维度（默认：128）
- `nhead`: 注意力头数（默认：8）
- `num_layers`: Transformer层数（默认：4）
- `dim_feedforward`: 前馈层维度（默认：512）
- `dropout`: Dropout率（默认：0.2）

### 训练参数
- `batch_size`: 批量大小（默认：1）
- `learning_rate`: 学习率（默认：5e-4）
- `num_epochs`: 训练轮数（默认：100）
- `pos_weight`: 正样本权重（默认：50.0，用于处理类别不平衡）
- `patience`: 早停耐心值（默认：20）

### 推理参数
- `threshold`: 分类阈值（默认：0.3）
- `min_region_bins`: 最小区间长度（默认：3）
- `top_n`: 返回top-N预测区间（默认：5）

## 评估指标

模型使用以下指标进行评估：

- **Precision（精确率）**: TP / (TP + FP)
- **Recall（召回率）**: TP / (TP + FN)
- **F1 Score**: 2 × (Precision × Recall) / (Precision + Recall)
- **IoU**: 预测区域与真实区域的交并比
- **AUC**: ROC曲线下面积

## 使用示例

### 训练自定义模型

```python
from src.training.config import Config
from src.training.model import create_model
from src.training.train import train

# 创建配置
config = Config()
config.training.num_epochs = 50
config.training.learning_rate = 1e-4

# 训练
model, metrics = train(config, data_dir="/path/to/data")
```

### 推理单条染色体

```python
from src.training.inference import load_model, predict_single_chromosome

# 加载模型
model, feature_stats, config = load_model("checkpoints/best_model.pt")

# 预测
result = predict_single_chromosome(
    model, 
    "path/to/chromosome.csv", 
    feature_stats, 
    config
)

print(f"预测区域: {result['predicted_regions']}")
```

## 性能优化建议

1. **类别不平衡**: 调整 `pos_weight` 参数（通常设置为正负样本比例）
2. **过拟合**: 增加 `dropout` 或使用早停机制
3. **欠拟合**: 增加模型容量（`d_model`, `num_layers`）或训练轮数
4. **内存不足**: 减小 `max_seq_len` 或使用梯度累积
5. **训练加速**: 使用GPU（`--device cuda`）

## 监控训练

使用TensorBoard监控训练过程：

```bash
tensorboard --logdir=training/logs
```

然后在浏览器访问 `http://localhost:6006`

## 常见问题

### Q: 训练时出现NaN损失？
A: 可能是学习率过高或数据归一化问题。尝试降低学习率或检查数据质量。

### Q: 模型预测全为0或全为1？
A: 调整 `pos_weight` 参数和分类阈值 `threshold`。

### Q: 如何处理超长序列？
A: 可以分段处理或增加 `max_seq_len` 参数（需要更多内存）。

### Q: 支持多GPU训练吗？
A: 当前版本使用单GPU，多GPU支持可通过 `torch.nn.DataParallel` 实现。

## 技术细节

### 位置编码
使用标准的正弦/余弦位置编码：

```python
PE(pos, 2i) = sin(pos / 10000^(2i/d_model))
PE(pos, 2i+1) = cos(pos / 10000^(2i/d_model))
```

### 损失函数
使用加权二元交叉熵损失（Weighted BCE Loss）处理类别不平衡：

```python
Loss = -[w_pos × y × log(ŷ) + w_neg × (1-y) × log(1-ŷ)]
```

### 数据归一化
特征使用Z-score标准化：

```python
X_norm = (X - mean) / std
```

## 引用

如果您在研究中使用了本项目，请引用：

```bibtex
@software{centromere_prediction,
  title = {Centromere Area Prediction with Transformer},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/centromere_prediction}
}
```

## 贡献指南

欢迎贡献！请遵循以下步骤：

1. Fork本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启Pull Request

## 许可证

本项目采用 MIT 许可证。详见 `LICENSE` 文件。

## 联系方式

- 项目主页: [GitHub Repository]
- 问题反馈: [Issues]

## 更新日志

### v1.0.0 (2024-12)
- 初始版本发布
- 实现基于Transformer的着丝粒预测模型
- 支持多尺度k-mer特征
- 提供完整的训练、推理和评估工具链

## 致谢

感谢所有为本项目提供帮助和建议的研究者。


