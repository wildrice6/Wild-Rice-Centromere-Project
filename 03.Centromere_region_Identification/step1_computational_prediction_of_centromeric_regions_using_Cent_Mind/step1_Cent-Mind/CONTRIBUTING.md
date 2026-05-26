# 贡献指南 / Contributing Guide

首先，感谢您考虑为本项目做出贡献！

## 如何贡献

### 报告Bug

如果您发现了Bug，请创建一个Issue并包含以下信息：
- Bug的详细描述
- 复现步骤
- 预期行为和实际行为
- 系统环境（操作系统、Python版本、PyTorch版本等）
- 如果可能，提供最小可复现示例

### 提出新功能

如果您有新功能的想法，请：
1. 先创建一个Issue讨论该功能
2. 说明功能的用途和预期效果
3. 等待维护者反馈后再开始实现

### 提交代码

1. **Fork仓库**
   ```bash
   # 在GitHub上点击Fork按钮
   # 然后克隆您的fork
   git clone https://github.com/YOUR_USERNAME/centromere_prediction.git
   cd centromere_prediction
   ```

2. **创建分支**
   ```bash
   git checkout -b feature/your-feature-name
   # 或
   git checkout -b fix/your-bug-fix
   ```

3. **进行更改**
   - 遵循现有的代码风格
   - 添加必要的注释和文档字符串
   - 如果可能，添加测试

4. **测试您的更改**
   ```bash
   # 运行基本测试
   python -m pytest tests/
   ```

5. **提交更改**
   ```bash
   git add .
   git commit -m "feat: add new feature X"
   # 或
   git commit -m "fix: resolve issue with Y"
   ```

   提交信息格式：
   - `feat:` 新功能
   - `fix:` Bug修复
   - `docs:` 文档更新
   - `style:` 代码格式（不影响功能）
   - `refactor:` 代码重构
   - `test:` 添加测试
   - `chore:` 构建/工具更新

6. **推送到GitHub**
   ```bash
   git push origin feature/your-feature-name
   ```

7. **创建Pull Request**
   - 访问您的GitHub仓库
   - 点击"New Pull Request"
   - 填写PR描述，说明您的更改
   - 等待代码审查

## 代码风格

### Python代码规范

- 遵循 PEP 8 规范
- 使用4个空格缩进
- 每行最多100个字符
- 使用有意义的变量名和函数名
- 为所有公共函数添加文档字符串

示例：
```python
def calculate_metrics(predictions: np.ndarray, labels: np.ndarray) -> dict:
    """
    计算评估指标
    
    Args:
        predictions: 预测值数组，形状为 (n_samples,)
        labels: 真实标签数组，形状为 (n_samples,)
    
    Returns:
        包含各项指标的字典
    """
    # 实现代码...
    return metrics
```

### 文档规范

- 为新功能添加文档
- 更新相关的README或文档文件
- 使用清晰的中英文说明
- 提供使用示例

## 开发环境设置

```bash
# 克隆仓库
git clone https://github.com/YOUR_USERNAME/centromere_prediction.git
cd centromere_prediction

# 创建虚拟环境
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 或
venv\Scripts\activate  # Windows

# 安装开发依赖
pip install -r requirements.txt
pip install -e .  # 以开发模式安装

# 安装开发工具（可选）
pip install pytest black flake8 mypy
```

## 测试

在提交PR之前，请确保：

```bash
# 运行代码格式检查
black src/ --check
flake8 src/

# 运行类型检查
mypy src/

# 运行测试
pytest tests/
```

## 代码审查流程

1. 提交PR后，维护者会进行代码审查
2. 可能会要求进行修改
3. 所有讨论解决后，PR会被合并
4. 您的贡献会被记录在项目历史中

## 行为准则

### 我们的承诺

为了营造开放和友好的环境，我们作为贡献者和维护者承诺：
- 尊重所有人
- 接受建设性的批评
- 专注于对社区最有利的事情
- 对其他社区成员表现出同理心

### 不可接受的行为

- 使用性化的语言或图像
- 人身攻击或政治攻击
- 公开或私下的骚扰
- 未经许可发布他人的私人信息
- 其他不专业或不受欢迎的行为

## 问题和讨论

- 对于Bug报告和功能请求，请使用GitHub Issues
- 对于一般性讨论，可以使用GitHub Discussions
- 对于紧急问题，可以发送邮件给维护者

## 许可

通过贡献代码，您同意您的贡献将在MIT许可证下发布。

## 联系方式

- GitHub Issues: [项目Issues页面]
- Email: your.email@example.com

## 致谢

感谢所有贡献者的支持和帮助！

您的名字会被添加到贡献者列表中。


