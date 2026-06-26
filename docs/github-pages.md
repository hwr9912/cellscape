# GitHub Pages 部署

本仓库使用 MkDocs 组织文档。文档源文件位于 `docs/`，站点配置位于仓库根目录的
`mkdocs.yml`。

## 本地预览

```bash
pip install mkdocs mkdocs-material
mkdocs serve
```

启动后访问终端显示的本地地址，通常是：

```text
http://127.0.0.1:8000
```

## 构建静态站点

```bash
mkdocs build
```

构建结果会输出到 `site/` 目录。

## 手动部署

如果你有仓库推送权限，可以运行：

```bash
mkdocs gh-deploy
```

该命令会把构建后的站点推送到 `gh-pages` 分支。

## GitHub Actions 自动部署

也可以在 GitHub Pages 中选择 GitHub Actions 作为部署来源，并添加自动部署 workflow。

示例 `.github/workflows/docs.yml`：

```yaml
name: docs

on:
  push:
    branches: [main]

permissions:
  contents: read
  pages: write
  id-token: write

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - run: pip install mkdocs mkdocs-material
      - run: mkdocs build
      - uses: actions/upload-pages-artifact@v3
        with:
          path: site
      - uses: actions/deploy-pages@v4
```
