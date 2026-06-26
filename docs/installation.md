# 安装

## 环境要求

`cellscape` 要求 Python `>=3.12`。

主要依赖：

- `anndata`
- `scanpy`
- `squidpy`
- `matplotlib`
- `numpy`
- `pandas`
- `scipy`
- `seaborn`
- `pillow`
- `glasbey`
- `tqdm`

## 从 GitHub 安装

```bash
pip install git+https://github.com/hwr9912/cellscape.git
```

## 本地开发安装

在仓库根目录运行：

```bash
pip install -e .
```

开发环境：

```bash
pip install -e ".[dev]"
```

可选交互依赖：

```bash
pip install -e ".[interactive]"
```

## 验证安装

```python
import cellscape as cs

print(cs.__version__)
```
