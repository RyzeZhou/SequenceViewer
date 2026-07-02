import os
import platform
import re
import sys
from dataclasses import dataclass

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QColor, QFont, QFontDatabase, QSyntaxHighlighter, QTextCharFormat
from PySide6.QtWidgets import (
    QApplication, QComboBox, QFileDialog, QLabel, QListWidget,
    QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QSpinBox, QSplitter,
    QStatusBar, QStyleFactory, QTabWidget, QTextEdit, QToolBar, QVBoxLayout, QWidget,
)


DNA = set("ACGTRYSWKMBDHVN-")
RNA = set("ACGURYSWKMBDHVN-")
PROTEIN = set("ABCDEFGHIKLMNPQRSTVWXYZ*-")
SMILES_ATOMS = {
    "B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "H",
    "Si", "Na", "Li", "K", "Mg", "Ca", "Fe", "Zn", "Cu", "Mn",
    "Al", "Se", "As", "Ag", "Au", "Hg", "Pb", "Sn",
}


@dataclass
class Record:
    label: str
    sequence: str
    line: int


class SequenceHighlighter(QSyntaxHighlighter):
    def __init__(self, document, mode_getter):
        super().__init__(document)
        self.mode_getter = mode_getter
        self.formats = {
            "header": self._fmt("#075985", bold=True),
            "comment": self._fmt("#6b7280"),
            "a": self._fmt("#047857"),
            "c": self._fmt("#2563eb"),
            "g": self._fmt("#b45309"),
            "t": self._fmt("#be123c"),
            "u": self._fmt("#7c3aed"),
            "ambiguous": self._fmt("#6d28d9"),
            "gap": self._fmt("#64748b"),
            "atom": self._fmt("#0f766e"),
            "bond": self._fmt("#7c2d12"),
            "branch": self._fmt("#1d4ed8"),
            "ring": self._fmt("#9333ea"),
            "error": self._fmt("#991b1b", bg="#fee2e2", underline=True),
            "index": self._fmt("#475569", bg="#f8fafc"),
            "ruler": self._fmt("#6b7280", bg="#f3f4f6"),
        }

    def _fmt(self, fg, bg=None, bold=False, underline=False):
        fmt = QTextCharFormat()
        fmt.setForeground(QColor(fg))
        if bg:
            fmt.setBackground(QColor(bg))
        if bold:
            fmt.setFontWeight(QFont.Bold)
        if underline:
            fmt.setUnderlineStyle(QTextCharFormat.SingleUnderline)
        return fmt

    def highlightBlock(self, text):
        mode = self.mode_getter()
        stripped = text.strip()
        if not stripped:
            return
        if stripped.startswith(">"):
            self.setFormat(0, len(text), self.formats["header"])
            return
        if stripped.startswith("#") or stripped.startswith(";"):
            self.setFormat(0, len(text), self.formats["comment"])
            return
        if self._is_ruler_line(text):
            self.setFormat(0, len(text), self.formats["ruler"])
            return

        seq_start, seq_end = self._sequence_span(text)
        if seq_start > 0:
            self.setFormat(0, seq_start, self.formats["index"])
        if seq_end < len(text):
            self.setFormat(seq_end, len(text) - seq_end, self.formats["index"])

        for i in range(seq_start, seq_end):
            char = text[i]
            tag = self._fasta_tag(char.upper()) if mode == "fasta" else self._smiles_tag(char)
            if tag:
                self.setFormat(i, 1, self.formats[tag])

    def _sequence_span(self, text):
        match = re.match(r"^\s*\d+\s{2,}(.+?)\s{2,}\d+\s*$", text)
        if match:
            return match.start(1), match.end(1)
        return 0, len(text)

    def _is_ruler_line(self, line):
        stripped = line.strip()
        return bool(stripped) and all(char.isdigit() or char.isspace() for char in stripped)

    def _fasta_tag(self, char):
        if char == "A": return "a"
        if char == "C": return "c"
        if char == "G": return "g"
        if char == "T": return "t"
        if char == "U": return "u"
        if char == "-": return "gap"
        if char in DNA or char in RNA or char in PROTEIN: return "ambiguous"
        return "error" if not char.isspace() else None

    def _smiles_tag(self, char):
        if char.isalpha(): return "atom"
        if char.isdigit() or char == "%": return "ring"
        if char in "()[]": return "branch"
        if char in "-=#$:/\\.": return "bond"
        return None


class SequenceViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.current_path = None
        self._syncing = False
        self._formatted_dirty = False
        self._last_mode = "fasta"
        self.setWindowTitle("序列查看器")
        self.resize(1240, 800)
        self._build_ui()
        self._connect()
        self.load_fasta_sample()

    def _build_ui(self):
        toolbar = QToolBar("工具")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)

        open_action = QAction("打开", self)
        save_action = QAction("保存", self)
        save_as_action = QAction("另存为", self)
        open_action.triggered.connect(self.open_file)
        save_action.triggered.connect(self.save_file)
        save_as_action.triggered.connect(self.save_as)
        toolbar.addAction(open_action)
        toolbar.addAction(save_action)
        toolbar.addAction(save_as_action)
        toolbar.addSeparator()

        toolbar.addWidget(QLabel("类型"))
        self.mode_box = QComboBox()
        self.mode_box.addItems(["auto", "fasta", "smiles"])
        toolbar.addWidget(self.mode_box)
        toolbar.addWidget(QLabel("每行字符"))
        self.wrap_spin = QSpinBox()
        self.wrap_spin.setRange(10, 500)
        self.wrap_spin.setSingleStep(10)
        self.wrap_spin.setValue(60)
        toolbar.addWidget(self.wrap_spin)
        toolbar.addSeparator()

        toolbar.addWidget(QLabel("序列字体"))
        self.font_box = QComboBox()
        self.font_box.addItems(self._available_monospace_fonts())
        toolbar.addWidget(self.font_box)
        toolbar.addWidget(QLabel("字号"))
        self.font_size_spin = QSpinBox()
        self.font_size_spin.setRange(8, 24)
        self.font_size_spin.setValue(11)
        toolbar.addWidget(self.font_size_spin)
        toolbar.addSeparator()

        fasta_btn = QPushButton("示例 FASTA")
        smiles_btn = QPushButton("示例 SMILES")
        fasta_btn.clicked.connect(self.load_fasta_sample)
        smiles_btn.clicked.connect(self.load_smiles_sample)
        toolbar.addWidget(fasta_btn)
        toolbar.addWidget(smiles_btn)

        root = QWidget()
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(10, 8, 10, 8)
        self.path_label = QLabel("未打开文件")
        root_layout.addWidget(self.path_label)

        splitter = QSplitter(Qt.Horizontal)
        root_layout.addWidget(splitter, 1)

        self.tabs = QTabWidget()
        self.raw_edit = QPlainTextEdit()
        self.raw_edit.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.formatted_edit = QPlainTextEdit()
        self.formatted_edit.setLineWrapMode(QPlainTextEdit.NoWrap)
        self._apply_sequence_font()
        self.tabs.addTab(self.raw_edit, "原始")
        self.tabs.addTab(self.formatted_edit, "格式化")
        splitter.addWidget(self.tabs)

        side = QWidget()
        side_layout = QVBoxLayout(side)
        side_layout.setContentsMargins(10, 0, 0, 0)
        side_layout.addWidget(QLabel("诊断"))
        self.diagnostics = QListWidget()
        side_layout.addWidget(self.diagnostics, 2)
        side_layout.addWidget(QLabel("统计"))
        self.stats = QTextEdit()
        self.stats.setReadOnly(True)
        self.stats.setMaximumHeight(150)
        side_layout.addWidget(self.stats)
        side_layout.addWidget(QLabel("记录"))
        self.records = QListWidget()
        side_layout.addWidget(self.records, 2)
        splitter.addWidget(side)
        splitter.setSizes([850, 320])

        self.setCentralWidget(root)
        self.setStatusBar(QStatusBar())
        self.raw_highlighter = SequenceHighlighter(self.raw_edit.document(), self.current_mode)
        self.formatted_highlighter = SequenceHighlighter(self.formatted_edit.document(), self.current_mode)
        self._apply_style()

    def _connect(self):
        self.raw_edit.textChanged.connect(self._on_raw_changed)
        self.formatted_edit.textChanged.connect(self._on_formatted_changed)
        self.tabs.currentChanged.connect(self._on_tab_changed)
        self.mode_box.currentTextChanged.connect(lambda _text: self.refresh_from_raw())
        self.wrap_spin.valueChanged.connect(lambda _value: self.refresh_from_raw())
        self.font_box.currentTextChanged.connect(lambda _text: self._apply_sequence_font())
        self.font_size_spin.valueChanged.connect(lambda _value: self._apply_sequence_font())
        self.diagnostics.itemClicked.connect(self._jump_from_item)
        self.records.itemClicked.connect(self._jump_from_item)

        self.refresh_timer = QTimer(self)
        self.refresh_timer.setSingleShot(True)
        self.refresh_timer.setInterval(250)
        self.refresh_timer.timeout.connect(self.refresh_from_raw)

    def _available_monospace_fonts(self):
        installed = set(QFontDatabase.families())
        preferred = [
            "Cascadia Mono", "Cascadia Code", "Consolas", "JetBrains Mono",
            "Fira Code", "Source Code Pro", "DejaVu Sans Mono", "Courier New",
            "Lucida Console", "NSimSun",
        ]
        found = [name for name in preferred if name in installed]
        return found or ["Consolas", "Courier New"]

    def _apply_sequence_font(self):
        if not hasattr(self, "raw_edit"):
            return
        family = self.font_box.currentText() if hasattr(self, "font_box") else "Consolas"
        size = self.font_size_spin.value() if hasattr(self, "font_size_spin") else 11
        font = QFont(family, size)
        font.setStyleHint(QFont.Monospace)
        font.setFixedPitch(True)
        font.setStyleStrategy(QFont.PreferAntialias)
        font.setHintingPreference(QFont.PreferNoHinting)
        self.raw_edit.setFont(font)
        self.formatted_edit.setFont(font)
        self.raw_edit.document().setDefaultFont(font)
        self.formatted_edit.document().setDefaultFont(font)
        self.raw_edit.viewport().setFont(font)
        self.formatted_edit.viewport().setFont(font)

    def _apply_style(self):
        self.setStyleSheet("""
            QWidget { background-color: #f5f5f5; color: #1e1e1e; font-size: 10pt; }
            QMainWindow { background-color: #fafafa; }
            QToolBar { background-color: #e8e8e8; border: none; border-bottom: 1px solid #d0d0d0; spacing: 6px; padding: 4px; }
            QPushButton, QToolButton { background-color: #ffffff; border: 1px solid #c0c0c0; border-radius: 4px; padding: 6px 12px; }
            QPushButton:hover, QToolButton:hover { background-color: #e8f0fe; border-color: #0078d4; }
            QComboBox, QSpinBox { background-color: #ffffff; border: 1px solid #c0c0c0; border-radius: 4px; padding: 5px; }
            QPlainTextEdit, QTextEdit, QListWidget { background-color: #ffffff; border: 1px solid #d0d0d0; selection-background-color: #0078d4; selection-color: white; }
            QTabWidget::pane { border: 1px solid #d0d0d0; background: #ffffff; }
            QTabBar::tab { background: #ececec; padding: 8px 18px; border: 1px solid #d0d0d0; border-bottom: none; }
            QTabBar::tab:selected { background: #ffffff; color: #0078d4; }
        """)

    def current_mode(self):
        text = self.mode_box.currentText()
        if text != "auto":
            return text
        useful = [line.strip() for line in self.raw_edit.toPlainText().splitlines() if line.strip()]
        return "fasta" if useful and useful[0].startswith(">") else "smiles"

    def _on_raw_changed(self):
        if self._syncing:
            return
        self._formatted_dirty = False
        self.refresh_timer.start()

    def _on_formatted_changed(self):
        if self._syncing:
            return
        self._formatted_dirty = True
        self.statusBar().showMessage("格式化视图已修改，切回原始视图或保存时会同步")

    def _on_tab_changed(self, index):
        if index == 0:
            self.sync_formatted_to_raw()
        elif not self._formatted_dirty:
            self.refresh_from_raw()

    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "打开序列文件", "",
            "序列文件 (*.fasta *.fa *.faa *.fna *.ffn *.frn *.smi *.smiles *.txt);;所有文件 (*.*)"
        )
        if path:
            self.load_file(path)

    def load_file(self, path):
        try:
            try:
                with open(path, "r", encoding="utf-8-sig") as handle:
                    content = handle.read()
            except UnicodeDecodeError:
                with open(path, "r", encoding="gbk", errors="replace") as handle:
                    content = handle.read()
        except OSError as exc:
            QMessageBox.critical(self, "无法打开文件", str(exc))
            return
        self.current_path = path
        self.path_label.setText(path)
        self._guess_mode(path, content)
        self._set_raw_text(content)

    def save_file(self):
        self.sync_formatted_to_raw()
        if not self.current_path:
            return self.save_as()
        try:
            with open(self.current_path, "w", encoding="utf-8", newline="") as handle:
                handle.write(self.raw_edit.toPlainText())
            self.statusBar().showMessage(f"已保存：{self.current_path}")
        except OSError as exc:
            QMessageBox.critical(self, "无法保存文件", str(exc))

    def save_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "另存为", "", "FASTA (*.fasta);;SMILES (*.smi);;文本 (*.txt);;所有文件 (*.*)")
        if path:
            self.current_path = path
            self.path_label.setText(path)
            self.save_file()

    def _guess_mode(self, path, content):
        ext = os.path.splitext(path)[1].lower()
        if ext in {".smi", ".smiles"}:
            self.mode_box.setCurrentText("smiles")
        elif ext in {".fasta", ".fa", ".faa", ".fna", ".ffn", ".frn"}:
            self.mode_box.setCurrentText("fasta")
        else:
            self.mode_box.setCurrentText("fasta" if content.lstrip().startswith(">") else "auto")

    def _set_raw_text(self, text):
        self._syncing = True
        self.raw_edit.setPlainText(text)
        self._syncing = False
        self._formatted_dirty = False
        self.refresh_from_raw()

    def refresh_from_raw(self):
        if self._syncing:
            return
        text = self.raw_edit.toPlainText()
        mode = self.current_mode()
        self._last_mode = mode
        if mode == "fasta":
            diagnostics, stats, records = self._analyze_fasta(text)
        else:
            diagnostics, stats, records = self._analyze_smiles(text)
        if not self._formatted_dirty:
            self._render_formatted(records)
        self._show_sidebars(diagnostics, stats, records)
        self.raw_highlighter.rehighlight()
        self.formatted_highlighter.rehighlight()
        self.statusBar().showMessage(f"{mode.upper()} | {len(diagnostics)} 个问题 | 每行 {self.wrap_spin.value()} 字符")

    def sync_formatted_to_raw(self):
        if not self._formatted_dirty:
            return
        mode = self._last_mode or self.current_mode()
        raw = self._formatted_to_raw(self.formatted_edit.toPlainText(), mode)
        self._syncing = True
        self.raw_edit.setPlainText(raw)
        self._syncing = False
        self._formatted_dirty = False
        self.refresh_from_raw()

    def _analyze_fasta(self, text):
        diagnostics, records = [], []
        seq_count = total = gc_count = 0
        current_header = None
        current_line = 1
        current_parts = []
        seen_before_header = False

        def finish():
            if current_header is not None:
                records.append(Record(current_header or "(空标题)", "".join(current_parts), current_line))

        for line_no, raw in enumerate(text.splitlines(), 1):
            stripped = raw.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            if stripped.startswith(">"):
                if current_header is not None and not current_parts:
                    diagnostics.append(("警告", line_no - 1, 1, "上一条 FASTA 记录没有序列"))
                finish()
                current_header = stripped[1:].strip()
                current_line = line_no
                current_parts = []
                seq_count += 1
                if not current_header:
                    diagnostics.append(("错误", line_no, 1, "FASTA 标题不能为空"))
                continue
            if current_header is None:
                seen_before_header = True
                diagnostics.append(("错误", line_no, 1, "序列行出现在第一个 FASTA 标题之前"))
            sequence = re.sub(r"\s+", "", stripped).upper()
            current_parts.append(sequence)
            for col, char in enumerate(stripped, 1):
                if not char.isspace() and not self._valid_fasta_char(char.upper()):
                    diagnostics.append(("错误", line_no, col, f"非法 FASTA 字符 '{char}'"))
            total += len(sequence.replace("-", ""))
            gc_count += sequence.count("G") + sequence.count("C")
        finish()
        if not text.strip():
            diagnostics.append(("提示", 1, 1, "请输入或打开 FASTA/SMILES 内容"))
        elif not seq_count and not seen_before_header:
            diagnostics.append(("错误", 1, 1, "没有找到 FASTA 标题行，标题应以 > 开头"))
        elif current_header is not None and not current_parts:
            diagnostics.append(("警告", max(1, len(text.splitlines())), 1, "最后一条 FASTA 记录没有序列"))
        gc = (gc_count / total * 100) if total else 0
        stats = f"类型：FASTA\n记录数：{seq_count}\n序列字符数：{total}\nGC 含量：{gc:.2f}%\n格式化：每行 {self.wrap_spin.value()} 字符"
        return diagnostics, stats, records

    def _valid_fasta_char(self, char):
        return char in DNA or char in RNA or char in PROTEIN

    def _analyze_smiles(self, text):
        diagnostics, records = [], []
        atom_count, ring_counts = 0, {}
        atoms_upper = {atom.upper() for atom in SMILES_ATOMS}
        for line_no, raw in enumerate(text.splitlines(), 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smiles = parts[0]
            label = " ".join(parts[1:]) or smiles
            records.append(Record(label, smiles, line_no))
            branch_stack = []
            i = 0
            while i < len(smiles):
                char = smiles[i]
                if char == "[":
                    close = smiles.find("]", i + 1)
                    if close == -1:
                        diagnostics.append(("错误", line_no, i + 1, "方括号没有闭合"))
                    else:
                        if not smiles[i + 1:close]:
                            diagnostics.append(("错误", line_no, i + 1, "空的原子方括号"))
                        atom_count += 1
                        i = close
                elif smiles[i:i + 2] in SMILES_ATOMS:
                    atom_count += 1
                    i += 1
                elif char.isalpha():
                    if char.upper() in atoms_upper or char in "bcnops":
                        atom_count += 1
                    else:
                        diagnostics.append(("错误", line_no, i + 1, f"未知原子符号 '{char}'"))
                elif char == "(":
                    branch_stack.append((line_no, i + 1))
                elif char == ")":
                    if branch_stack:
                        branch_stack.pop()
                    else:
                        diagnostics.append(("错误", line_no, i + 1, "多余的右圆括号"))
                elif char.isdigit():
                    ring_counts[char] = ring_counts.get(char, 0) + 1
                elif char == "%":
                    token = smiles[i + 1:i + 3]
                    if len(token) == 2 and token.isdigit():
                        ring_counts[token] = ring_counts.get(token, 0) + 1
                        i += 2
                    else:
                        diagnostics.append(("错误", line_no, i + 1, "% 后面应跟两位环编号"))
                elif char not in "-=#$:/\\.[]+-":
                    diagnostics.append(("错误", line_no, i + 1, f"非法 SMILES 字符 '{char}'"))
                i += 1
            for branch_line, col in branch_stack:
                diagnostics.append(("错误", branch_line, col, "圆括号没有闭合"))
        for ring, count in sorted(ring_counts.items()):
            if count % 2:
                diagnostics.append(("错误", 1, 1, f"环编号 {ring} 出现奇数次，可能没有闭合"))
        if not text.strip():
            diagnostics.append(("提示", 1, 1, "请输入或打开 FASTA/SMILES 内容"))
        stats = f"类型：SMILES\n分子/记录数：{len(records)}\n原子 token 数：{atom_count}\n环编号种类：{len(ring_counts)}\n格式化：每行 {self.wrap_spin.value()} 字符"
        return diagnostics, stats, records

    def _render_formatted(self, records):
        width = self.wrap_spin.value()
        max_len = max((len(record.sequence) for record in records), default=0)
        index_width = max(8, len(str(max_len)) + 2)
        lines = []
        for record in records:
            lines.append(f">{record.label}")
            lines.append(" " * index_width + "  " + self._make_ruler(width))
            if not record.sequence:
                lines.append(" " * index_width + "  (空序列)")
            for offset in range(0, len(record.sequence), width):
                chunk = record.sequence[offset:offset + width]
                start = offset + 1
                end = offset + len(chunk)
                lines.append(f"{start:>{index_width}}  {chunk}  {end:<{index_width}}")
            lines.append("")
        self._syncing = True
        self.formatted_edit.setPlainText("\n".join(lines))
        self._syncing = False

    def _formatted_to_raw(self, text, mode):
        records = []
        current_header = None
        chunks = []
        expected_width = self.wrap_spin.value()
        index_width = 8

        def finish():
            nonlocal index_width
            if current_header is not None:
                records.append((current_header, "".join(chunks)))

        for raw_line in text.splitlines():
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                finish()
                current_header = stripped[1:].strip()
                chunks = []
                index_width = 8
                continue
            if stripped == "(空序列)" or self._is_ruler_line(raw_line):
                if self._is_ruler_line(raw_line):
                    index_width = self._detect_index_width(raw_line)
                continue
            chunk = self._extract_sequence_column(raw_line, mode, expected_width, index_width)
            if chunk:
                if current_header is None:
                    current_header = "sequence"
                chunks.append(chunk)
        finish()

        if mode == "smiles":
            return "\n".join(f"{seq} {label}".rstrip() for label, seq in records) + ("\n" if records else "")
        lines = []
        width = self.wrap_spin.value()
        for label, seq in records:
            lines.append(f">{label}")
            for offset in range(0, len(seq), width):
                lines.append(seq[offset:offset + width])
        return "\n".join(lines) + ("\n" if lines else "")

    def _is_ruler_line(self, line):
        stripped = line.strip()
        return bool(stripped) and all(char.isdigit() or char.isspace() for char in stripped)

    def _detect_index_width(self, ruler_line):
        match = re.match(r"^(\s*)\S", ruler_line)
        if not match:
            return 8
        return max(1, len(match.group(1)) - 2)

    def _extract_sequence_column(self, line, mode, expected_width, index_width):
        fixed_start = index_width + 2
        if re.match(r"^\s*\d+\s{2,}", line) and len(line) > fixed_start:
            tail = line[fixed_start:]
            text = re.sub(r"\s{2,}\d+\s*$", "", tail)[:expected_width]
        else:
            match = re.match(r"^\s*\d+\s{2,}(\S+)(?:\s{2,}\d+\s*)?$", line)
            if match:
                text = match.group(1)[:expected_width]
            else:
                text = line.strip()[:expected_width]
        text = re.sub(r"\s+", "", text)
        return text.upper() if mode == "fasta" else text

    def _make_ruler(self, width):
        chars = [" "] * width
        for number in [1] + list(range(10, width + 1, 10)):
            label = str(number)
            start = max(0, number - len(label))
            for offset, char in enumerate(label):
                if start + offset < width:
                    chars[start + offset] = char
        return "".join(chars)

    def _show_sidebars(self, diagnostics, stats, records):
        self.diagnostics.clear()
        if not diagnostics:
            self.diagnostics.addItem("没有发现语法问题")
        else:
            for level, line, col, msg in diagnostics:
                item_text = f"{level} L{line}:C{col} {msg}"
                self.diagnostics.addItem(item_text)
                self.diagnostics.item(self.diagnostics.count() - 1).setData(Qt.UserRole, (line, col))
        self.stats.setPlainText(stats)
        self.records.clear()
        for record in records:
            self.records.addItem(f"L{record.line}  {record.label}")
            self.records.item(self.records.count() - 1).setData(Qt.UserRole, (record.line, 1))

    def _jump_from_item(self, item):
        pos = item.data(Qt.UserRole)
        if not pos:
            return
        line, col = pos
        self.tabs.setCurrentIndex(0)
        cursor = self.raw_edit.textCursor()
        block = self.raw_edit.document().findBlockByNumber(max(0, line - 1))
        cursor.setPosition(block.position() + max(0, col - 1))
        self.raw_edit.setTextCursor(cursor)
        self.raw_edit.setFocus()

    def load_fasta_sample(self):
        self.current_path = None
        self.path_label.setText("示例 FASTA")
        self.mode_box.setCurrentText("fasta")
        self._set_raw_text(">seq1 human example\nATGCGTACGTAGCTAGCTAGNNN---ATGCGTACGTAGCTAGCTAGNNN---\n>seq2 protein example\nMEEPQSDPSVEPPLSQETFSDLWKLLPEN\n")

    def load_smiles_sample(self):
        self.current_path = None
        self.path_label.setText("示例 SMILES")
        self.mode_box.setCurrentText("smiles")
        self._set_raw_text("CCO ethanol\nc1ccccc1 benzene\nCC(=O)O acetic_acid\nC1CC invalid_unclosed_ring\n")


def main():
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )
    app = QApplication(sys.argv)
    app.setStyle(QStyleFactory.create("Fusion"))
    if platform.system() == "Windows":
        font = QFont("Microsoft YaHei", 10)
        font.setStyleStrategy(QFont.PreferAntialias)
        font.setHintingPreference(QFont.PreferNoHinting)
        app.setFont(font)
    window = SequenceViewer()
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        window.load_file(sys.argv[1])
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
