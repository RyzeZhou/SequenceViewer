# 序列查看器

一个 Python / PySide6 桌面 GUI，用于查看和检查 FASTA、SMILES 文本。

## 功能

- FASTA / SMILES 语法检查和高亮
- 原始视图和可编辑的格式化视图
- 按指定字符数折行显示序列
- 左右显示每行起止字符序号
- 上方显示 1、10、20、30 这样的坐标尺
- 在格式化视图修改后，切回原始视图或保存时会自动同步
- 参考百宝箱的 PySide6/Fusion/高 DPI 字体平滑设置，改善中文显示
- 序列编辑区支持切换等宽字体和字号
- 支持作为 Windows “打开方式”打开 FASTA 文件

## 依赖

```powershell
pip install PySide6
```

## 运行

```powershell
python .\sequence_viewer.pyw
```

直接打开文件：

```powershell
python .\sequence_viewer.pyw .\example.fasta
```

## 设置为 FASTA 打开方式

在 PowerShell 中运行：

```powershell
powershell -ExecutionPolicy Bypass -File .\register_open_with.ps1
```

然后右键 `.fasta` / `.fa` 文件，选择“打开方式”，选择“序列查看器”。
