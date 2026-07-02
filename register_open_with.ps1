$ErrorActionPreference = "Stop"

$appPath = Join-Path $PSScriptRoot "sequence_viewer.pyw"
if (-not (Test-Path $appPath)) {
    throw "找不到 sequence_viewer.pyw：$appPath"
}

$pythonw = Get-Command pythonw.exe -ErrorAction SilentlyContinue
if (-not $pythonw) {
    $python = Get-Command python.exe -ErrorAction SilentlyContinue
    if (-not $python) {
        throw "找不到 pythonw.exe 或 python.exe。请先安装 Python，并勾选 Add Python to PATH。"
    }
    $pythonwPath = $python.Source
} else {
    $pythonwPath = $pythonw.Source
}

$command = "`"$pythonwPath`" `"$appPath`" `"%1`""
$appKey = "HKCU:\Software\Classes\Applications\SequenceViewer.pyw"

New-Item -Path "$appKey\shell\open\command" -Force | Out-Null
Set-ItemProperty -Path "$appKey\shell\open\command" -Name "(default)" -Value $command
Set-ItemProperty -Path $appKey -Name "FriendlyAppName" -Value "序列查看器"

$extensions = ".fasta", ".fa", ".faa", ".fna", ".ffn", ".frn"
foreach ($ext in $extensions) {
    $openWithKey = "HKCU:\Software\Classes\$ext\OpenWithList\SequenceViewer.pyw"
    New-Item -Path $openWithKey -Force | Out-Null
}

Write-Host "已添加到 FASTA 文件的“打开方式”：序列查看器"
Write-Host "如果没有立刻显示，请右键 FASTA 文件 -> 打开方式 -> 选择其他应用 -> 更多应用。"
