# 从FASTA到预测结果 - 完整流程指南

本文档介绍如何从原始的FASTA基因组文件出发，一步步得到着丝粒预测结果。

## 流程概览

```
FASTA基因组文件
    ↓
步骤1: k-mer分析
    ↓
步骤2: 生成特征CSV
    ↓
步骤3: 模型推理
    ↓
步骤4: 查看结果
```

## 准备工作

### 所需文件
- ✅ 基因组FASTA文件（例如：`genome.fasta`）
- ✅ 训练好的模型（`best_model.pt`）

### 所需工具
```bash
# 安装本项目
pip install -r requirements.txt

# 安装k-mer分析工具（二选一）
# 选项1: Jellyfish (推荐)
conda install -c bioconda jellyfish

# 选项2: KMC
conda install -c bioconda kmc
```

## 完整流程

### 步骤1: k-mer频率统计

使用Jellyfish计算不同k值的k-mer频率：

```bash
# 为您的基因组文件设置变量
GENOME="genome.fasta"
OUTPUT_DIR="kmer_analysis"
mkdir -p $OUTPUT_DIR

# 计算4个k值的k-mer频率
for k in 64 128 256 512; do
    echo "Processing k=${k}..."
    
    # 计数k-mers
    jellyfish count \
        -m $k \
        -s 1G \
        -t 8 \
        -C \
        -o ${OUTPUT_DIR}/${k}mer.jf \
        $GENOME
    
    # 导出为文本格式
    jellyfish dump \
        ${OUTPUT_DIR}/${k}mer.jf \
        > ${OUTPUT_DIR}/${k}mer_counts.txt
    
    echo "k=${k} completed"
done
```

### 步骤2: 生成特征CSV文件

创建预处理脚本（或使用项目提供的工具）：

```python
# generate_features.py
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import argparse

def load_kmer_counts(kmer_file):
    """加载k-mer计数文件"""
    kmer_counts = {}
    with open(kmer_file, 'r') as f:
        while True:
            kmer = f.readline().strip()
            if not kmer:
                break
            count = int(f.readline().strip())
            kmer_counts[kmer] = count
    return kmer_counts

def calculate_bin_statistics(genome_file, kmer_counts_dict, bin_size=10000):
    """
    计算每个bin的统计特征
    
    Args:
        genome_file: FASTA文件路径
        kmer_counts_dict: {k_value: kmer_counts}的字典
        bin_size: bin大小（默认10kb）
    
    Returns:
        DataFrame with features
    """
    results = []
    
    # 读取基因组序列
    for record in SeqIO.parse(genome_file, "fasta"):
        seq = str(record.seq).upper()
        seq_len = len(seq)
        
        # 按bin分割
        for start in range(0, seq_len, bin_size):
            end = min(start + bin_size, seq_len)
            bin_seq = seq[start:end]
            
            if len(bin_seq) < bin_size // 2:  # 跳过太短的bin
                continue
            
            bin_features = {
                'chromosome': record.id,
                'start': start,
                'end': end,
                'has_cen': 0  # 推理时设为0，如果有真实标签可修改
            }
            
            # 计算每个k值的特征
            for k, kmer_counts in kmer_counts_dict.items():
                if len(bin_seq) < k:
                    bin_features[f'{k}_highlighted_percent'] = 0.0
                    bin_features[f'{k}_coverage_depth_avg'] = 0.0
                    continue
                
                # 提取bin中的所有k-mers
                bin_kmers = []
                for i in range(len(bin_seq) - k + 1):
                    kmer = bin_seq[i:i+k]
                    if 'N' not in kmer:  # 跳过含N的k-mer
                        bin_kmers.append(kmer)
                
                if not bin_kmers:
                    bin_features[f'{k}_highlighted_percent'] = 0.0
                    bin_features[f'{k}_coverage_depth_avg'] = 0.0
                    continue
                
                # 计算统计量
                counts = [kmer_counts.get(kmer, 0) for kmer in bin_kmers]
                
                # highlighted_percent: 高频k-mer的比例
                # 定义"高频"为大于中位数的2倍
                median_count = np.median(counts) if counts else 0
                threshold = median_count * 2
                highlighted = sum(1 for c in counts if c > threshold)
                highlighted_percent = highlighted / len(counts) if counts else 0
                
                # coverage_depth_avg: 平均覆盖深度
                avg_depth = np.mean(counts) if counts else 0
                
                bin_features[f'{k}_highlighted_percent'] = highlighted_percent
                bin_features[f'{k}_coverage_depth_avg'] = avg_depth
            
            results.append(bin_features)
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description='Generate feature CSV from FASTA and k-mer counts')
    parser.add_argument('--genome', required=True, help='Input FASTA file')
    parser.add_argument('--kmer-dir', required=True, help='Directory containing k-mer count files')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--bin-size', type=int, default=10000, help='Bin size (default: 10000)')
    
    args = parser.parse_args()
    
    print("Loading k-mer counts...")
    kmer_counts_dict = {}
    for k in [64, 128, 256, 512]:
        kmer_file = f"{args.kmer_dir}/{k}mer_counts.txt"
        print(f"  Loading k={k}...")
        kmer_counts_dict[k] = load_kmer_counts(kmer_file)
    
    print("Calculating bin statistics...")
    df = calculate_bin_statistics(args.genome, kmer_counts_dict, args.bin_size)
    
    print(f"Saving to {args.output}...")
    df.to_csv(args.output, index=False)
    
    print(f"Done! Generated {len(df)} bins")
    print(f"\nOutput columns:")
    print(df.columns.tolist())
    print(f"\nFirst few rows:")
    print(df.head())

if __name__ == '__main__':
    main()
```

运行特征提取：

```bash
python generate_features.py \
    --genome genome.fasta \
    --kmer-dir kmer_analysis \
    --output genome_multi_k_summary.csv \
    --bin-size 10000
```

### 步骤3: 模型推理

使用训练好的模型进行预测：

```bash
cd src/training

python inference.py \
    --checkpoint ../../checkpoints/best_model.pt \
    --input ../../genome_multi_k_summary.csv \
    --output ../../predictions \
    --threshold 0.3 \
    --device cuda
```

**参数说明**：
- `--checkpoint`: 模型文件路径
- `--input`: 步骤2生成的CSV文件
- `--output`: 结果保存目录
- `--threshold`: 分类阈值（0.1-0.5，越小召回率越高）
- `--device`: 使用cuda或cpu

### 步骤4: 查看结果

#### 4.1 查看JSON详细结果

```bash
# 查看预测的详细信息
cat predictions/predictions.json
```

JSON格式示例：
```json
{
  "csv_file": "genome_multi_k_summary.csv",
  "seq_len": 1000,
  "predictions": [0.05, 0.08, 0.12, 0.89, 0.92, 0.88, 0.15, ...],
  "predicted_regions": [
    {
      "start_bin": 350,
      "end_bin": 450,
      "start_pos": 3500000,
      "end_pos": 4500000,
      "length_bins": 100,
      "length_bp": 1000000,
      "avg_prob": 0.87,
      "max_prob": 0.95
    }
  ]
}
```

#### 4.2 查看CSV汇总结果

```bash
# 表格形式查看预测区域
cat predictions/predictions_summary.csv
```

CSV格式示例：
```csv
file,seq_len,num_regions,top_region_start,top_region_end,top_region_prob
genome_multi_k_summary.csv,1000,3,3500000,4500000,0.8700
```

#### 4.3 提取预测区域为BED格式

创建转换脚本：

```python
# predictions_to_bed.py
import json
import sys

def json_to_bed(json_file, bed_file, min_prob=0.5):
    """将预测结果转换为BED格式"""
    with open(json_file, 'r') as f:
        results = json.load(f)
    
    with open(bed_file, 'w') as f:
        # 如果是单个染色体的结果
        if isinstance(results, dict):
            results = [results]
        
        for result in results:
            for i, region in enumerate(result['predicted_regions']):
                if region['avg_prob'] >= min_prob:
                    # BED格式: chr start end name score
                    chrom = result.get('chromosome', 'chr1')
                    start = region['start_pos']
                    end = region['end_pos']
                    name = f"centromere_{i+1}"
                    score = int(region['avg_prob'] * 1000)
                    
                    f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\n")

if __name__ == '__main__':
    json_to_bed(
        'predictions/predictions.json',
        'predictions/centromeres.bed',
        min_prob=0.5
    )
    print("BED file created: predictions/centromeres.bed")
```

运行：
```bash
python predictions_to_bed.py
```

## 简化流程（一键脚本）

将所有步骤整合到一个脚本中：

```bash
#!/bin/bash
# predict_from_fasta.sh - 从FASTA到预测结果的一键脚本

set -e

# 参数
GENOME=$1
MODEL=$2
OUTPUT_DIR=${3:-"predictions_output"}
BIN_SIZE=${4:-10000}
THREADS=${5:-8}

if [ -z "$GENOME" ] || [ -z "$MODEL" ]; then
    echo "Usage: $0 <genome.fasta> <model.pt> [output_dir] [bin_size] [threads]"
    echo "Example: $0 genome.fasta best_model.pt results 10000 8"
    exit 1
fi

echo "========================================="
echo "Centromere Prediction Pipeline"
echo "========================================="
echo "Genome: $GENOME"
echo "Model: $MODEL"
echo "Output: $OUTPUT_DIR"
echo "Bin size: $BIN_SIZE"
echo "Threads: $THREADS"
echo "========================================="

# 创建输出目录
mkdir -p $OUTPUT_DIR
KMER_DIR="$OUTPUT_DIR/kmer_analysis"
mkdir -p $KMER_DIR

# 步骤1: k-mer统计
echo ""
echo "[Step 1/4] Computing k-mer frequencies..."
for k in 64 128 256 512; do
    echo "  Processing k=$k..."
    jellyfish count -m $k -s 1G -t $THREADS -C \
        -o ${KMER_DIR}/${k}mer.jf $GENOME
    jellyfish dump ${KMER_DIR}/${k}mer.jf \
        > ${KMER_DIR}/${k}mer_counts.txt
    rm ${KMER_DIR}/${k}mer.jf  # 删除中间文件节省空间
done
echo "  k-mer analysis completed!"

# 步骤2: 生成特征
echo ""
echo "[Step 2/4] Generating feature CSV..."
python generate_features.py \
    --genome $GENOME \
    --kmer-dir $KMER_DIR \
    --output ${OUTPUT_DIR}/features.csv \
    --bin-size $BIN_SIZE
echo "  Feature CSV created!"

# 步骤3: 模型推理
echo ""
echo "[Step 3/4] Running model inference..."
python src/training/inference.py \
    --checkpoint $MODEL \
    --input ${OUTPUT_DIR}/features.csv \
    --output ${OUTPUT_DIR}/predictions \
    --threshold 0.3
echo "  Inference completed!"

# 步骤4: 生成BED文件
echo ""
echo "[Step 4/4] Generating BED file..."
python predictions_to_bed.py
echo "  BED file created!"

# 总结
echo ""
echo "========================================="
echo "Pipeline completed successfully!"
echo "========================================="
echo "Results:"
echo "  - Feature CSV: ${OUTPUT_DIR}/features.csv"
echo "  - Predictions JSON: ${OUTPUT_DIR}/predictions/predictions.json"
echo "  - Summary CSV: ${OUTPUT_DIR}/predictions/predictions_summary.csv"
echo "  - BED file: ${OUTPUT_DIR}/predictions/centromeres.bed"
echo "========================================="
```

使用一键脚本：

```bash
chmod +x predict_from_fasta.sh
./predict_from_fasta.sh genome.fasta checkpoints/best_model.pt
```

## 结果解读

### 预测区域信息

| 字段 | 说明 |
|------|------|
| start_pos / end_pos | 着丝粒在基因组上的位置（bp） |
| length_bp | 着丝粒长度 |
| avg_prob | 平均预测概率（0-1） |
| max_prob | 最高预测概率（0-1） |

### 可信度评估

- **avg_prob > 0.8**: 高置信度，很可能是着丝粒
- **avg_prob 0.5-0.8**: 中等置信度，较可能是着丝粒
- **avg_prob < 0.5**: 低置信度，可能是假阳性

### 可视化结果

使用IGV或其他基因组浏览器查看：

```bash
# 在IGV中加载
# 1. 加载参考基因组
# 2. 加载预测的BED文件: predictions/centromeres.bed
# 3. 查看预测的着丝粒区域
```

## 完整示例

假设您有一个拟南芥基因组：

```bash
# 1. 下载或准备基因组
# genome.fasta (已有)

# 2. 下载预训练模型
# best_model.pt (已有)

# 3. 运行完整流程
./predict_from_fasta.sh \
    genome.fasta \
    checkpoints/best_model.pt \
    arabidopsis_results

# 4. 查看结果
cat arabidopsis_results/predictions/predictions_summary.csv

# 5. 在IGV中可视化
# File -> Load from File -> arabidopsis_results/predictions/centromeres.bed
```

## 性能和资源

### 计算时间估算

对于一个200Mb的基因组：
- k-mer统计: ~10-30分钟
- 特征生成: ~5-10分钟
- 模型推理: ~1-5分钟（GPU）或~10-30分钟（CPU）
- **总计**: ~20-60分钟

### 内存需求

- k-mer统计: ~2-8GB（取决于基因组大小）
- 特征生成: ~1-4GB
- 模型推理: ~2-4GB（GPU）或~1-2GB（CPU）

### 磁盘空间

- k-mer中间文件: ~1-5GB（可在生成CSV后删除）
- 特征CSV: ~10-100MB
- 预测结果: ~1-10MB

## 常见问题

### Q: 没有Jellyfish怎么办？

A: 可以使用KMC或其他k-mer计数工具，但需要修改脚本的相应部分。

### Q: 基因组很大（>1GB）怎么办？

A: 
1. 增加Jellyfish的内存限制（-s参数）
2. 分染色体处理
3. 使用更多线程加速

### Q: 如何批量处理多个基因组？

A:
```bash
for genome in *.fasta; do
    name=$(basename $genome .fasta)
    ./predict_from_fasta.sh $genome checkpoints/best_model.pt ${name}_results
done
```

### Q: 预测结果不理想怎么办？

A:
1. 调整阈值参数（--threshold）
2. 检查特征CSV的数值范围是否合理
3. 如果是新物种，可能需要重新训练模型

### Q: 能否跳过k-mer统计步骤？

A: 如果您已有k-mer统计结果，可以直接从步骤2开始，但必须确保格式兼容。

## 进阶使用

### 自定义bin大小

```bash
# 使用5kb的bin
python generate_features.py \
    --genome genome.fasta \
    --kmer-dir kmer_analysis \
    --output features_5kb.csv \
    --bin-size 5000
```

### 调整预测阈值

```bash
# 更敏感（更高召回率，可能更多假阳性）
python inference.py --threshold 0.2 ...

# 更保守（更高精确率，可能漏掉部分区域）
python inference.py --threshold 0.5 ...
```

### 分染色体处理

```bash
# 如果基因组很大，可以分染色体处理
for chr in chr*.fasta; do
    ./predict_from_fasta.sh $chr model.pt ${chr}_results
done

# 合并结果
cat */predictions/centromeres.bed > all_centromeres.bed
```

## 获取帮助

- 📖 查看 [快速开始指南](QUICKSTART_CN.md)
- 🔧 查看 [数据格式文档](DATA_FORMAT.md)
- 💬 提交 [GitHub Issue](https://github.com/yourusername/centromere_prediction/issues)

---

**最后更新**: 2024-12-19

