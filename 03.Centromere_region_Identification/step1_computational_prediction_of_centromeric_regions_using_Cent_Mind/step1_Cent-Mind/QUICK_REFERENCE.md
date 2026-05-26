# 快速参考 - Quick Reference

## 三种使用方式

### 🚀 方式1: 从FASTA一键推理（最简单）

**适用场景**: 您有基因组FASTA文件，想直接得到着丝粒预测结果

```bash
# 一键运行完整流程
./scripts/predict_from_fasta.sh genome.fasta checkpoints/best_model.pt

# 查看结果
cat predictions_output/predictions/predictions_summary.csv

# 在IGV中可视化
# 加载: predictions_output/predictions/centromeres.bed
```

**需要**: 
- ✅ 基因组FASTA文件
- ✅ 预训练模型文件
- ✅ Jellyfish工具: `conda install -c bioconda jellyfish`

**详细教程**: [docs/FROM_FASTA_TO_PREDICTION.md](docs/FROM_FASTA_TO_PREDICTION.md)

---

### 📊 方式2: 从特征CSV推理

**适用场景**: 您已经有处理好的特征CSV文件

```bash
cd src/training
python inference.py \
    --checkpoint ../../checkpoints/best_model.pt \
    --input your_features.csv \
    --output ./predictions
```

**需要**:
- ✅ 特征CSV文件（包含8个k-mer特征列）
- ✅ 预训练模型文件

**CSV格式**: [docs/DATA_FORMAT.md](docs/DATA_FORMAT.md)

---

### 🎓 方式3: 训练自己的模型

**适用场景**: 您有标注好的训练数据，想训练新模型

```bash
cd src/training
python train.py --data_dir /path/to/data --device cuda

# 监控训练
tensorboard --logdir=logs
```

**需要**:
- ✅ 标注的训练数据（多个CSV文件）
- ✅ GPU（可选，但强烈推荐）

**详细教程**: [docs/QUICKSTART_CN.md](docs/QUICKSTART_CN.md)

---

## 常用命令速查

### 安装依赖

```bash
# Python依赖
pip install -r requirements.txt

# k-mer分析工具（用于从FASTA推理）
conda install -c bioconda jellyfish

# 安装为Python包
pip install -e .
```

### 从FASTA生成特征CSV

```bash
# 步骤1: k-mer统计
for k in 64 128 256 512; do
    jellyfish count -m $k -s 1G -t 8 -C -o ${k}mer.jf genome.fasta
    jellyfish dump ${k}mer.jf > ${k}mer_counts.txt
done

# 步骤2: 生成特征
python src/preprocessing/generate_features.py \
    --genome genome.fasta \
    --kmer-dir . \
    --output features.csv
```

### 模型推理

```bash
# 基本推理
python src/training/inference.py \
    --checkpoint checkpoints/best_model.pt \
    --input features.csv \
    --output predictions

# 调整阈值（更敏感）
python src/training/inference.py \
    --checkpoint checkpoints/best_model.pt \
    --input features.csv \
    --output predictions \
    --threshold 0.2

# 使用CPU
python src/training/inference.py \
    --checkpoint checkpoints/best_model.pt \
    --input features.csv \
    --output predictions \
    --device cpu
```

### 转换为BED格式

```bash
# 基本转换
python src/postprocessing/predictions_to_bed.py \
    predictions/predictions.json \
    centromeres.bed

# 只保留高置信度区域
python src/postprocessing/predictions_to_bed.py \
    predictions/predictions.json \
    centromeres.bed \
    --min-prob 0.7

# 只保留top 3区域
python src/postprocessing/predictions_to_bed.py \
    predictions/predictions.json \
    centromeres.bed \
    --top-n 3

# 生成详细BED
python src/postprocessing/predictions_to_bed.py \
    predictions/predictions.json \
    centromeres_detail.bed \
    --detailed
```

### 评估结果

```bash
# 生成汇总报告
python src/evaluation/generate_summary_report.py predictions/

# 计算平均指标
python src/evaluation/calculate_average_metrics.py predictions/

# Top-5评估
python src/evaluation/evaluate_top5_prediction.py predictions/
```

### 模型训练

```bash
# 基本训练
python src/training/train.py --data_dir /path/to/data

# 自定义参数
python src/training/train.py \
    --data_dir /path/to/data \
    --epochs 100 \
    --lr 5e-4 \
    --pos_weight 50.0 \
    --device cuda

# 监控训练
tensorboard --logdir=src/training/logs
```

## 输出文件说明

### 推理输出

| 文件 | 说明 |
|------|------|
| `predictions.json` | 详细预测结果，包含每个位置的概率 |
| `predictions_summary.csv` | 表格汇总，包含预测区域信息 |
| `centromeres.bed` | BED格式，可在IGV等工具中查看 |
| `centromeres_detailed.bed` | 详细BED，包含概率等信息 |

### 训练输出

| 文件/目录 | 说明 |
|----------|------|
| `checkpoints/best_model.pt` | 验证集上最佳模型 |
| `checkpoints/final_model.pt` | 最终训练模型 |
| `outputs/training_results_*.json` | 训练结果统计 |
| `logs/` | TensorBoard日志 |

## 参数速查表

### inference.py参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| --checkpoint | 必需 | 模型文件路径 |
| --input | 必需 | 输入CSV文件或目录 |
| --output | ./predictions | 输出目录 |
| --threshold | 0.3 | 分类阈值（0.1-0.5） |
| --device | 自动检测 | cuda或cpu |

### train.py参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| --data_dir | 必需 | 训练数据目录 |
| --epochs | 100 | 训练轮数 |
| --lr | 5e-4 | 学习率 |
| --pos_weight | 50.0 | 正样本权重 |
| --device | 自动检测 | cuda或cpu |

### generate_features.py参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| --genome | 必需 | FASTA文件 |
| --kmer-dir | 必需 | k-mer计数目录 |
| --output | 必需 | 输出CSV文件 |
| --bin-size | 10000 | Bin大小(bp) |
| --chromosome | 无 | 只处理指定染色体 |

### predict_from_fasta.sh参数

```bash
./scripts/predict_from_fasta.sh <genome.fasta> <model.pt> [output_dir] [bin_size] [threads] [threshold]
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| genome.fasta | 必需 | 基因组文件 |
| model.pt | 必需 | 模型文件 |
| output_dir | predictions_output | 输出目录 |
| bin_size | 10000 | Bin大小 |
| threads | 8 | 线程数 |
| threshold | 0.3 | 预测阈值 |

## 常见问题快速解答

### Q: 我只有FASTA文件，怎么办？
A: 使用方式1 - `./scripts/predict_from_fasta.sh genome.fasta model.pt`

### Q: 推理很慢怎么办？
A: 1) 使用GPU: `--device cuda`; 2) 增加线程数

### Q: 预测结果不准确？
A: 1) 调整阈值 `--threshold`; 2) 如果是新物种，需要重新训练

### Q: 如何批量处理多个基因组？
```bash
for fasta in *.fasta; do
    ./scripts/predict_from_fasta.sh $fasta model.pt ${fasta%.fasta}_results
done
```

### Q: 如何在IGV中查看？
A: File → Load from File → 选择 `centromeres.bed`

### Q: 内存不够？
A: 1) 分染色体处理; 2) 使用 `--chromosome` 参数

### Q: 没有Jellyfish？
A: `conda install -c bioconda jellyfish` 或使用KMC

## 资源链接

- 📖 [完整文档](README.md)
- 🚀 [快速开始](docs/QUICKSTART_CN.md)  
- 🧬 [从FASTA到预测](docs/FROM_FASTA_TO_PREDICTION.md)
- 📊 [数据格式](docs/DATA_FORMAT.md)
- 🏗️ [模型架构](docs/MODEL_ARCHITECTURE.md)
- 🤝 [贡献指南](CONTRIBUTING.md)
- 📝 [更新日志](CHANGELOG.md)

## 联系支持

- GitHub Issues: [提交问题](https://github.com/yourusername/centromere_prediction/issues)
- Email: your.email@example.com

---

**提示**: 如果这是您第一次使用，强烈推荐从[快速开始指南](docs/QUICKSTART_CN.md)开始！

