# 数据格式说明

本文档详细说明模型所需的输入数据格式。

## 概述

模型接收CSV格式的输入文件，每个文件代表一条染色体或基因组片段。文件名必须以 `_multi_k_summary.csv` 结尾。

## CSV文件格式

### 必需列

| 列名 | 数据类型 | 说明 | 示例值 |
|------|---------|------|--------|
| `start` | int | 区间起始位置（碱基对） | 0 |
| `end` | int | 区间结束位置（碱基对） | 10000 |
| `has_cen` | int | 是否为着丝粒（0或1） | 1 |
| `64_highlighted_percent` | float | k=64时高亮百分比 | 0.85 |
| `64_coverage_depth_avg` | float | k=64时平均覆盖深度 | 15.2 |
| `128_highlighted_percent` | float | k=128时高亮百分比 | 0.82 |
| `128_coverage_depth_avg` | float | k=128时平均覆盖深度 | 14.5 |
| `256_highlighted_percent` | float | k=256时高亮百分比 | 0.78 |
| `256_coverage_depth_avg` | float | k=256时平均覆盖深度 | 13.8 |
| `512_highlighted_percent` | float | k=512时高亮百分比 | 0.75 |
| `512_coverage_depth_avg` | float | k=512时平均覆盖深度 | 12.9 |

### 列说明

#### 位置列
- **start/end**: 定义基因组上的区间
  - 通常为等长的bins（如10kb）
  - 必须连续且不重叠
  - 单位：碱基对（bp）

#### 标签列
- **has_cen**: 着丝粒标签
  - 0: 非着丝粒区域
  - 1: 着丝粒区域
  - 训练时必需，推理时可选

#### 特征列
对于每个k值（64, 128, 256, 512），有两个特征：

1. **highlighted_percent**: 高亮百分比
   - 范围: 0.0 到 1.0
   - 表示该区间内被标记为"高亮"的k-mer比例
   - 着丝粒区域通常有较高的值

2. **coverage_depth_avg**: 平均覆盖深度
   - 范围: 通常 0 到 100+
   - 表示该区间的平均测序深度或k-mer计数
   - 着丝粒区域通常有较高的值

### 示例CSV文件

```csv
start,end,has_cen,64_highlighted_percent,64_coverage_depth_avg,128_highlighted_percent,128_coverage_depth_avg,256_highlighted_percent,256_coverage_depth_avg,512_highlighted_percent,512_coverage_depth_avg
0,10000,0,0.15,2.3,0.12,1.8,0.10,1.5,0.08,1.2
10000,20000,0,0.18,2.5,0.15,2.1,0.13,1.9,0.10,1.5
20000,30000,1,0.85,15.2,0.82,14.5,0.78,13.8,0.75,12.9
30000,40000,1,0.90,16.5,0.88,15.8,0.85,15.1,0.82,14.3
40000,50000,1,0.87,15.8,0.84,15.2,0.80,14.5,0.77,13.7
50000,60000,0,0.16,2.4,0.13,2.0,0.11,1.7,0.09,1.3
```

## 文件命名规范

### 命名格式

```
{sample_id}_{chromosome_id}_multi_k_summary.csv
```

### 示例文件名

- `sample1_chr1_multi_k_summary.csv`
- `HG002_chr2_multi_k_summary.csv`
- `human_genome_hap1_chr3_multi_k_summary.csv`
- `arabidopsis_chr5_multi_k_summary.csv`

### 命名要求

1. 必须以 `_multi_k_summary.csv` 结尾
2. 建议包含样本ID和染色体ID便于识别
3. 避免使用空格和特殊字符
4. 使用小写字母和下划线

## 数据组织

### 目录结构示例

```
data/
├── training_data/
│   ├── sample1_chr1_multi_k_summary.csv
│   ├── sample1_chr2_multi_k_summary.csv
│   ├── sample2_chr1_multi_k_summary.csv
│   └── sample2_chr2_multi_k_summary.csv
├── test_data/
│   ├── sample3_chr1_multi_k_summary.csv
│   └── sample3_chr2_multi_k_summary.csv
└── validation_data/
    ├── sample4_chr1_multi_k_summary.csv
    └── sample4_chr2_multi_k_summary.csv
```

### 数据划分建议

- **按样本划分**: 将不同个体的数据分到不同集合
- **按染色体划分**: 每个染色体作为独立样本
- **避免数据泄露**: 同一样本的不同染色体不应同时出现在训练和测试集

## 数据质量要求

### 基本要求

1. **完整性**: 所有必需列都存在
2. **一致性**: 同一文件内的bin大小应相同
3. **连续性**: start和end应该连续不重叠
4. **有效性**: 数值在合理范围内

### 数据检查清单

- [ ] CSV格式正确，可被pandas读取
- [ ] 所有必需列都存在
- [ ] start < end 对所有行成立
- [ ] has_cen只包含0和1
- [ ] 特征值在合理范围内（无NaN或Inf）
- [ ] 文件名以_multi_k_summary.csv结尾
- [ ] 每个文件至少有10行数据

### 数据验证脚本

```python
import pandas as pd
import numpy as np

def validate_csv(filepath):
    """验证CSV文件格式"""
    required_cols = [
        'start', 'end', 'has_cen',
        '64_highlighted_percent', '64_coverage_depth_avg',
        '128_highlighted_percent', '128_coverage_depth_avg',
        '256_highlighted_percent', '256_coverage_depth_avg',
        '512_highlighted_percent', '512_coverage_depth_avg'
    ]
    
    try:
        df = pd.read_csv(filepath)
        
        # 检查列
        for col in required_cols:
            if col not in df.columns:
                return False, f"Missing column: {col}"
        
        # 检查start < end
        if not (df['start'] < df['end']).all():
            return False, "start must be < end"
        
        # 检查has_cen
        if not df['has_cen'].isin([0, 1]).all():
            return False, "has_cen must be 0 or 1"
        
        # 检查NaN
        if df[required_cols].isnull().any().any():
            return False, "Contains NaN values"
        
        return True, "Valid"
    
    except Exception as e:
        return False, str(e)

# 使用示例
valid, msg = validate_csv("your_file.csv")
print(f"Validation: {valid}, {msg}")
```

## 数据预处理

### 归一化

模型内部会自动进行Z-score归一化：

```python
X_normalized = (X - mean) / std
```

其中mean和std从训练集计算得到。

### 注意事项

1. **不需要手动归一化**: 模型会自动处理
2. **保持原始尺度**: 提供原始的特征值
3. **一致性**: 确保所有文件使用相同的计算方法

## 数据获取

### 从BAM文件生成

如果您有原始的BAM比对文件，可以使用以下工具生成所需的CSV：

```bash
# 伪代码示例
# 实际工具需要根据您的数据处理流程调整

# 1. 计算k-mer覆盖度
for k in 64 128 256 512; do
    jellyfish count -m $k -s 100M -t 8 genome.fasta
    jellyfish dump mer_counts.jf > ${k}mer_counts.txt
done

# 2. 计算每个bin的统计值
python calculate_bin_statistics.py \
    --genome genome.fasta \
    --kmers 64 128 256 512 \
    --bin-size 10000 \
    --output chr1_multi_k_summary.csv
```

### 从已有数据转换

如果您有其他格式的数据，需要转换为本格式：

```python
import pandas as pd

# 示例：从多个文件合并
df_64 = pd.read_csv('64mer_stats.csv')
df_128 = pd.read_csv('128mer_stats.csv')
df_256 = pd.read_csv('256mer_stats.csv')
df_512 = pd.read_csv('512mer_stats.csv')

# 合并
result = pd.DataFrame({
    'start': df_64['start'],
    'end': df_64['end'],
    'has_cen': df_64['is_centromere'],
    '64_highlighted_percent': df_64['highlight_pct'],
    '64_coverage_depth_avg': df_64['coverage'],
    '128_highlighted_percent': df_128['highlight_pct'],
    '128_coverage_depth_avg': df_128['coverage'],
    '256_highlighted_percent': df_256['highlight_pct'],
    '256_coverage_depth_avg': df_256['coverage'],
    '512_highlighted_percent': df_512['highlight_pct'],
    '512_coverage_depth_avg': df_512['coverage'],
})

result.to_csv('chr1_multi_k_summary.csv', index=False)
```

## 常见问题

### Q: bin大小必须是10kb吗？

A: 不必须。可以是任意大小，但建议：
- 不要太小（<1kb）：特征不稳定
- 不要太大（>100kb）：分辨率不够
- 10kb是推荐值，平衡了性能和分辨率

### Q: 不同染色体的bin大小可以不同吗？

A: 可以，但建议保持一致以便模型学习统一的模式。

### Q: 特征值的量纲不统一怎么办？

A: 模型会自动归一化，不用担心。但要确保同一特征在所有文件中使用相同的计算方法。

### Q: 缺失部分k值的数据怎么办？

A: 所有8个特征列都是必需的。如果缺失某些k值，需要补充或使用其他k值替代。

### Q: 推理时必须提供has_cen列吗？

A: 不必须，但如果提供了，模型会计算评估指标。可以填-1表示未知。

### Q: 可以使用其他k值吗？

A: 当前模型训练时使用64/128/256/512。如果要使用其他k值，需要：
1. 修改数据格式
2. 修改模型输入维度
3. 重新训练模型

## 示例数据集

项目仓库中提供了示例数据（如果有的话）：

```
examples/data/
├── sample_chr1_multi_k_summary.csv
├── sample_chr2_multi_k_summary.csv
└── README.md
```

## 联系和支持

如果对数据格式有疑问：
- 查看项目文档
- 提交GitHub Issue
- 发送邮件至维护者

---

**最后更新**: 2024-12-19


