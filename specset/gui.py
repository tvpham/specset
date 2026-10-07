"""PySide6 graphical interface for specset."""

import contextlib
import html
import math
import re
import sys
import traceback
from pathlib import Path

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QObject, QPoint, QSortFilterProxyModel, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPalette, QPen, QTextCursor
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QLineEdit, QListWidget, QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton,
    QSpinBox, QSplitter, QStackedWidget, QTabWidget, QTableView, QVBoxLayout, QWidget,
)

# Field kinds
DIR = "dir"              # existing directory
FILE = "file"            # existing file
OUT_DIR = "out_dir"      # output directory (may not exist yet)
OUT_FILE = "out_file"    # output file
INT = "int"
OPT_INT = "opt_int"      # 0 means None ("All")
CHOICE = "choice"
BOOL = "bool"
KEY_VALUES = "key_values"  # one KEY=VALUE per line -> dict

_INT_MAX = 2**31 - 1

# Each command: utils function name, title, description, and fields as
# (keyword argument, label, kind, options).
COMMANDS = [
    {
        "func": "maxquant_to_specset",
        "title": "MaxQuant → specset",
        "help": "Convert a MaxQuant txt directory (msms.txt and evidence.txt) to a specset.",
        "fields": [
            ("txt_dir", "MaxQuant txt directory", DIR, {}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("additional_columns", "Additional columns", KEY_VALUES, {"placeholder": "msms.txt column=specset name, one per line\ne.g. Gene Names=genes"}),
            ("max_files_in_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("n", "Max rows (-1 = all)", INT, {"default": -1, "min": -1}),
            ("overwrite", "Overwrite output", BOOL, {}),
        ],
    },
    {
        "func": "easypqp_to_specset",
        "title": "EasyPQP → specset",
        "help": "Convert an EasyPQP library TSV to a specset.",
        "fields": [
            ("easypqp_file", "EasyPQP library", FILE, {"filter": "TSV files (*.tsv *.txt);;All files (*)"}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("max_no_files_per_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("overwrite", "Overwrite output", BOOL, {}),
        ],
    },
    {
        "func": "diannlib_to_specset",
        "title": "DIA-NN library → specset",
        "help": "Convert a DIA-NN parquet library to a specset.",
        "fields": [
            ("diannlib_file", "DIA-NN library", FILE, {"filter": "Parquet files (*.parquet);;All files (*)"}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("max_no_files_per_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("overwrite", "Overwrite output", BOOL, {}),
        ],
    },
    {
        "func": "diannreport_to_specset",
        "title": "DIA-NN report → specset",
        "help": "Convert a DIA-NN parquet report to a specset. With 'One specset per run', create one specset per run, "
                "named <Run>.specset in the output directory.",
        "fields": [
            ("diannreport", "DIA-NN report", FILE, {"filter": "Parquet files (*.parquet);;All files (*)"}),
            ("out_dir", "Output specset", OUT_DIR, {"parent_if": "by_Run",
                                                    "tooltip": "With 'One specset per run', the directory in which the <Run>.specset directories are created."}),
            ("by_Run", "One specset per run", BOOL, {}),
            ("max_no_files_per_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("overwrite", "Overwrite output", BOOL, {"tooltip": "Delete an output specset first if it exists and is not empty. If unchecked, a non-empty output specset is an error."}),
        ],
    },
    {
        "func": "psmlist_to_specset",
        "title": "PSM list → specset",
        "help": "Convert a list of files containing stripped MaxQuant PSMs to a specset.",
        "fields": [
            ("psmlist_file", "PSM list file", FILE, {"filter": "Text files (*.txt *.tsv);;All files (*)"}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("max_files_in_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("overwrite", "Overwrite output", BOOL, {}),
        ],
    },
    {
        "func": "specset_to_speclib",
        "title": "Specset → spectral library",
        "help": "Convert a specset to a spectral library TSV.",
        "fields": [
            ("specset_dir", "Specset", DIR, {}),
            ("lib_name", "Output library", OUT_FILE, {"filter": "TSV files (*.tsv);;All files (*)"}),
            ("additional_columns", "Additional columns", KEY_VALUES, {"placeholder": "specset field=library column, one per line\ne.g. proteins=Protein.Ids"}),
            ("consensus_method", "Consensus method", CHOICE, {"choices": ["", "mean", "median"], "editable": True,
                                                              "tooltip": "Empty: PSM with most fragments. 'mean', 'median', or a specset field to select the PSM with the max value."}),
            ("n", "Max peptides", OPT_INT, {}),
        ],
    },
    {
        "func": "merge_specsets",
        "title": "Merge specsets",
        "help": "Merge specset A and specset B into a new specset. A and B are left unchanged.",
        "fields": [
            ("specset_a", "Specset A", DIR, {}),
            ("specset_b", "Specset B", DIR, {}),
            ("specset_out", "Output specset", OUT_DIR, {}),
            ("mode", "Mode", CHOICE, {"choices": ["both", "ignore"],
                                      "tooltip": "For peptides in both: 'both' keeps the PSMs from A and B, 'ignore' skips the ones from B."}),
            ("max_files_in_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("flush_threshold", "Flush threshold", INT, {"default": 1000, "min": 1}),
            ("overwrite", "Overwrite output", BOOL, {}),
        ],
    },
    {
        "func": "summarize",
        "title": "Summarize",
        "help": "Print summary statistics for one specset, or two specsets with their peptide overlap.",
        "fields": [
            ("specset_a", "Specset A", DIR, {}),
            ("specset_b", "Specset B (optional)", DIR, {"optional": True}),
        ],
    },
]


class _PathEdit(QWidget):
    """Line edit with a Browse button."""

    def __init__(self, kind, file_filter="", parent=None):
        super().__init__(parent)
        self.kind = kind
        self.file_filter = file_filter or "All files (*)"
        self.edit = QLineEdit()
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.edit)
        layout.addWidget(browse)

    def _browse(self):
        start = self.edit.text() or str(Path.home())
        if self.kind in (DIR, OUT_DIR):
            path = QFileDialog.getExistingDirectory(self, "Select directory", start)
        elif self.kind == FILE:
            path, _ = QFileDialog.getOpenFileName(self, "Select file", start, self.file_filter)
        else:
            path, _ = QFileDialog.getSaveFileName(self, "Save as", start, self.file_filter)
        if path:
            self.edit.setText(path)

    def text(self):
        return self.edit.text().strip()


class _CommandPage(QWidget):
    """Form for one command."""

    run_requested = Signal(str, dict)

    def __init__(self, command, parent=None):
        super().__init__(parent)
        self.command = command
        self.widgets = {}

        layout = QVBoxLayout(self)
        title = QLabel(f"<h2>{command['title']}</h2><p>{command['help']}</p>")
        title.setWordWrap(True)
        layout.addWidget(title)

        form = QFormLayout()
        for name, label, kind, opts in command["fields"]:
            if kind in (DIR, FILE, OUT_DIR, OUT_FILE):
                w = _PathEdit(kind, opts.get("filter", ""))
            elif kind in (INT, OPT_INT):
                w = QSpinBox()
                w.setRange(opts.get("min", 0), _INT_MAX)
                if kind == OPT_INT:
                    w.setSpecialValueText("All")
                w.setValue(opts.get("default", 0))
            elif kind == CHOICE:
                w = QComboBox()
                w.addItems(opts["choices"])
                w.setEditable(opts.get("editable", False))
            elif kind == BOOL:
                w = QCheckBox()
                w.setChecked(opts.get("default", False))
                if name == "overwrite":
                    w.setToolTip("Delete the output first if it exists and is not empty. If unchecked, a non-empty output is an error.")
            elif kind == KEY_VALUES:
                w = QPlainTextEdit()
                w.setPlaceholderText(opts.get("placeholder", "KEY=VALUE, one per line"))
                w.setFixedHeight(70)
            else:
                raise ValueError(f"unknown field kind: {kind}")
            if "tooltip" in opts:
                w.setToolTip(opts["tooltip"])
            self.widgets[name] = (kind, opts, w)
            form.addRow(label, w)
        layout.addLayout(form)

        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._on_run)
        row = QHBoxLayout()
        row.addStretch()
        row.addWidget(self.run_button)
        layout.addLayout(row)
        layout.addStretch()

    def values(self):
        """Collect keyword arguments; raise ValueError with a message for invalid input."""
        kwargs = {}
        for name, (kind, opts, w) in self.widgets.items():
            if kind in (DIR, FILE, OUT_DIR, OUT_FILE):
                path = w.text()
                if not path:
                    if opts.get("optional"):
                        kwargs[name] = None
                        continue
                    raise ValueError(f"Please fill in '{self._label(name)}'.")
                if kind == DIR and not Path(path).is_dir():
                    raise ValueError(f"Directory does not exist: {path}")
                if kind == FILE and not Path(path).is_file():
                    raise ValueError(f"File does not exist: {path}")
                kwargs[name] = path
            elif kind == INT:
                kwargs[name] = w.value()
            elif kind == OPT_INT:
                kwargs[name] = w.value() or None
            elif kind == CHOICE:
                kwargs[name] = w.currentText().strip() or None
            elif kind == BOOL:
                kwargs[name] = w.isChecked()
            elif kind == KEY_VALUES:
                d = {}
                for line in w.toPlainText().splitlines():
                    line = line.strip()
                    if not line:
                        continue
                    if "=" not in line:
                        raise ValueError(f"'{self._label(name)}': expected KEY=VALUE, got '{line}'")
                    key, value = line.split("=", 1)
                    d[key.strip()] = value.strip()
                kwargs[name] = d
        return kwargs

    def output_dir(self, kwargs):
        """Return (output directory, whether it is a parent of the output specsets), or (None, False)."""
        for name, (kind, opts, _) in self.widgets.items():
            if kind == OUT_DIR:
                return kwargs.get(name), opts.get("parent", False) or bool(kwargs.get(opts.get("parent_if", "")))
        return None, False

    def _label(self, name):
        return next(label for n, label, _, _ in self.command["fields"] if n == name)

    def _on_run(self):
        try:
            kwargs = self.values()
        except ValueError as e:
            QMessageBox.warning(self, "Invalid input", str(e))
            return
        self.run_requested.emit(self.command["func"], kwargs)


class _Stream(QObject):
    """File-like object that forwards writes as a Qt signal (thread safe)."""

    text_written = Signal(str)

    def write(self, s):
        if s:
            self.text_written.emit(s)
        return len(s)

    def flush(self):
        pass


class _Worker(QObject):
    """Runs one specset.utils function in a background thread."""

    finished = Signal(bool, str)

    def __init__(self, func_name, kwargs, stream):
        super().__init__()
        self.func_name = func_name
        self.kwargs = kwargs
        self.stream = stream

    def run(self):
        try:
            with contextlib.redirect_stdout(self.stream), contextlib.redirect_stderr(self.stream):
                from . import utils
                getattr(utils, self.func_name)(**self.kwargs)
        except Exception:
            self.finished.emit(False, traceback.format_exc())
        else:
            self.finished.emit(True, "")


_FRAG_COLS = ["FragmentLossType", "FragmentNumber", "FragmentType", "FragmentCharge", "RelativeIntensity"]


def _read_psm_file(path):
    """Read a PSM file into a dict of field -> value (string)."""
    psm = {}
    with open(path) as f:
        for line in f:
            key, _, value = line.rstrip("\n").partition("\t")
            if key:
                psm[key] = value
    return psm


def _parse_ms2(ms2):
    """Parse the ms2 field into a list of fragment dicts."""
    ms2 = ms2.strip()
    if not ms2:
        return []
    rows = [item.split(",") for item in ms2.rstrip(";").split(";") if item]
    frags = []
    for row in rows:
        if len(row) != len(_FRAG_COLS):
            continue
        loss, number, ftype, charge, intensity = row
        try:
            frags.append({"type": str(ftype), "number": int(number), "charge": int(charge),
                          "loss": str(loss), "intensity": float(intensity)})
        except ValueError:
            continue
    return frags


_aa_mass = None


def _fragment_mz(pep, frags):
    """Add an 'mz' key to each fragment; None where it cannot be calculated."""
    global _aa_mass
    try:
        from .mz import generate_aa_mass, get_ion_mz
        from .utils import _to_mz_seq
        if _aa_mass is None:
            _aa_mass = generate_aa_mass()
        seq = _to_mz_seq(pep)
    except Exception:
        seq = None
    for frag in frags:
        frag["mz"] = None
        if seq is None:
            continue
        try:
            frag["mz"] = get_ion_mz(seq, frag["type"], frag["number"], frag["charge"], frag["loss"], _aa_mass)
        except Exception:
            pass


def _to_int(s):
    """Return s as int if possible (for numeric sorting), else unchanged."""
    try:
        return int(s)
    except ValueError:
        return s


def _fragment_label(frag):
    label = f"{frag['type']}{frag['number']}"
    if frag["loss"] and frag["loss"] != "noloss":
        label += f"-{frag['loss']}"
    if frag["charge"] > 1:
        label += "+" * frag["charge"]
    return label


def _read_index(specset_path):
    """Read _index.tsv of a specset into a list of (peptide, location)."""
    rows = []
    with open(Path(specset_path) / "_index.tsv") as f:
        for line in f:
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                rows.append((parts[0], parts[1].strip()))
    return rows


def _read_location_index(location_path):
    """Read __index.tsv of a location into a list of (file name, charge)."""
    rows = []
    with open(Path(location_path) / "__index.tsv") as f:
        for line in f:
            parts = line.split()
            if parts:
                rows.append((parts[0], _to_int(parts[1]) if len(parts) > 1 else None))
    return rows


def _load_psm(psm_path):
    """Read a PSM file; return (fields, fragments with m/z). Raises OSError or ValueError."""
    psm = _read_psm_file(psm_path)
    frags = _parse_ms2(psm.get("ms2", ""))
    QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
    try:
        _fragment_mz(psm.get("pep", ""), frags)
    finally:
        QApplication.restoreOverrideCursor()
    return psm, frags


class _TableModel(QAbstractTableModel):
    """Read-only table model over a list of row tuples."""

    def __init__(self, headers, parent=None):
        super().__init__(parent)
        self.headers = headers
        self.rows = []

    def set_rows(self, rows):
        self.beginResetModel()
        self.rows = rows
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.headers)

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        value = self.rows[index.row()][index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            if value is None:
                return ""
            if isinstance(value, float):
                return f"{value:.4f}"
            return str(value)
        if role == Qt.ItemDataRole.UserRole:  # sort key
            return value if value is not None else ""
        if role == Qt.ItemDataRole.TextAlignmentRole and isinstance(value, (int, float)):
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return self.headers[section]
        return None


def _make_table(model, sortable=True):
    """Create a table view (with a sort/filter proxy) for a _TableModel."""
    proxy = QSortFilterProxyModel()
    proxy.setSourceModel(model)
    proxy.setSortRole(Qt.ItemDataRole.UserRole)
    proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
    view = QTableView()
    view.setModel(proxy)
    view.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    view.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    view.setSortingEnabled(sortable)
    if sortable:
        view.sortByColumn(-1, Qt.SortOrder.AscendingOrder)  # show rows in model order until a header is clicked
    view.setAlternatingRowColors(True)
    view.setWordWrap(False)
    view.verticalHeader().hide()
    view.verticalHeader().setDefaultSectionSize(view.fontMetrics().height() + 6)
    view.horizontalHeader().setStretchLastSection(True)
    return view, proxy


class _SpectrumWidget(QWidget):
    """Stick plot of fragment intensities, by m/z when available.

    With a second list of fragments, draws a mirror plot: the first PSM pointing
    up and the second pointing down, each scaled to its own maximum intensity.
    """

    _COLORS = {"b": QColor(31, 119, 180), "y": QColor(214, 39, 40)}
    _MAX_LABELS = 25

    def __init__(self, parent=None):
        super().__init__(parent)
        self.frags = []
        self.mirror = None
        self.captions = ("A", "B")
        self.setMinimumHeight(220)

    def set_fragments(self, frags, mirror=None):
        self.frags = frags
        self.mirror = mirror
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        painter.fillRect(self.rect(), pal.color(QPalette.ColorRole.Base))
        text_color = pal.color(QPalette.ColorRole.Text)

        fm = painter.fontMetrics()
        left, right, top, bottom = fm.horizontalAdvance("100%") + 12, 16, fm.height() + 12, fm.height() + 14
        w = self.width() - left - right
        h = self.height() - top - bottom
        if w <= 10 or h <= 10:
            return

        mirrored = self.mirror is not None
        sides = [(1, self.frags)] + ([(-1, self.mirror)] if mirrored else [])
        baseline = top + h / 2 if mirrored else top + h
        half = h / 2 if mirrored else h

        all_frags = [f for _, frags in sides for f in frags]
        painter.setPen(QPen(text_color, 1))
        if not all_frags:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No fragments")
            return

        painter.drawLine(left, int(baseline), left + w, int(baseline))
        painter.drawLine(left, top, left, top + h)

        by_mz = all(f["mz"] is not None for f in all_frags)
        if by_mz:
            def xpos(f):
                return f["mz"]
            x_min, x_max = 0.0, max(f["mz"] for f in all_frags) * 1.05
            x_label = "m/z"
        else:
            # no m/z: spread fragments evenly by type and number; the same fragment gets the same position on both sides
            keys = sorted({(f["type"], f["number"], f["charge"], f["loss"]) for f in all_frags})
            positions = {key: i + 1 for i, key in enumerate(keys)}

            def xpos(f):
                return positions[(f["type"], f["number"], f["charge"], f["loss"])]
            x_min, x_max = 0.0, len(keys) + 1.0
            x_label = "fragment (m/z not available)"

        def px(x):
            return left + (x - x_min) / (x_max - x_min) * w

        # y axis ticks
        for sign, _ in sides:
            for frac in (0, 0.5, 1.0):
                if sign < 0 and frac == 0:
                    continue
                y = baseline - sign * frac * half
                painter.drawLine(left - 4, int(y), left, int(y))
                label = f"{int(frac * 100)}%"
                painter.drawText(left - 6 - fm.horizontalAdvance(label), int(y + fm.ascent() / 2 - 1), label)

        # x axis ticks
        if by_mz:
            step = _nice_step(x_max - x_min, max(w // 90, 2))
            tick = step
            while tick < x_max:
                x = px(tick)
                painter.drawLine(int(x), int(baseline), int(x), int(baseline) + 4)
                label = f"{tick:g}"
                painter.drawText(int(x - fm.horizontalAdvance(label) / 2), top + h + 4 + fm.ascent(), label)
                tick += step
        painter.drawText(left + w - fm.horizontalAdvance(x_label), self.height() - 2, x_label)

        # captions for the mirror plot
        if mirrored:
            painter.drawText(left + 6, top + fm.ascent(), self.captions[0])
            painter.drawText(left + 6, top + h - 4, self.captions[1])

        # peaks
        for sign, frags in sides:
            if not frags:
                continue
            y_max = max(f["intensity"] for f in frags) or 1.0
            label_set = set(sorted(range(len(frags)), key=lambda i: -frags[i]["intensity"])[:self._MAX_LABELS])
            for i, frag in enumerate(frags):
                color = self._COLORS.get(frag["type"], text_color)
                x = px(xpos(frag))
                y = baseline - sign * frag["intensity"] / y_max * half
                painter.setPen(QPen(color, 2))
                painter.drawLine(int(x), int(baseline), int(x), int(y))
                if i in label_set:
                    label_y = max(y - 3, fm.ascent()) if sign > 0 else min(y + fm.ascent() + 1, top + h)
                    painter.drawText(int(x + 2), int(label_y), _fragment_label(frag))


def _nice_step(span, n_ticks):
    raw = span / n_ticks
    magnitude = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    for m in (1, 2, 5, 10):
        if raw <= m * magnitude:
            return m * magnitude
    return 10 * magnitude


class _PsmPanel(QGroupBox):
    """PSMs of one peptide (__index.tsv) aligned with the Fields and Fragments tables of the selected PSM."""

    fragments_changed = Signal(list)

    def __init__(self, title, parent=None):
        super().__init__(title, parent)
        self.location_path = None
        self._count = ""

        self.label = QLabel()
        self.label.setTextFormat(Qt.TextFormat.RichText)
        self.label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.label.setWordWrap(True)
        self.psm_model = _TableModel(["File"])  # rows are (file name, charge); charge is used for filtering
        self.psm_view, self.psm_proxy = _make_table(self.psm_model)
        self.field_model = _TableModel(["Field", "Value"])
        self.field_view, _ = _make_table(self.field_model, sortable=False)
        self.frag_model = _TableModel(["Fragment", "Type", "Number", "Charge", "Loss", "m/z", "Intensity"])
        self.frag_view, _ = _make_table(self.frag_model)
        self.psm_view.selectionModel().currentRowChanged.connect(self._on_psm_selected)

        self.tables = tables = QSplitter(Qt.Orientation.Horizontal)
        tables.addWidget(self.psm_view)
        tables.addWidget(self.field_view)
        tables.addWidget(self.frag_view)
        # PSM list just wide enough for an 8-digit file name
        psm_width = (self.psm_view.fontMetrics().horizontalAdvance("00000000") + 16 + 2 * self.psm_view.frameWidth()
                     + self.psm_view.verticalScrollBar().sizeHint().width())
        tables.setSizes([psm_width, 330, 500])

        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(tables, 1)
        self.clear()

    def align_fragments_to(self, widget):
        """Resize the Fields table so that the Fragments table starts at the left edge of widget."""
        window = self.window()
        delta = widget.mapTo(window, QPoint(0, 0)).x() - self.frag_view.mapTo(window, QPoint(0, 0)).x()
        sizes = self.tables.sizes()
        if sizes[1] + delta > 0 and sizes[2] - delta > 0:
            sizes[1] += delta
            sizes[2] -= delta
            self.tables.setSizes(sizes)

    def set_location(self, location_path, charge=None):
        """Show the PSMs of a peptide location, optionally only those with the given charge, and select the first one."""
        self.clear()
        self.location_path = Path(location_path)
        try:
            rows = _read_location_index(self.location_path)
        except OSError as e:
            self._set_details(f"Cannot read __index.tsv: {html.escape(str(e))}")
            return
        if charge is not None:
            rows = [row for row in rows if row[1] == charge]
            self._count = f"{len(rows)} PSMs with charge {charge}"
        else:
            self._count = f"{len(rows)} PSMs"
        self.psm_model.set_rows(rows)
        if rows:
            self.psm_view.selectRow(0)
        else:
            self._set_details("")

    def clear(self):
        self.location_path = None
        self._count = ""
        self.psm_model.set_rows([])
        self._clear_psm()
        self._set_details("Select a peptide")

    def _clear_psm(self):
        self.field_model.set_rows([])
        self.frag_model.set_rows([])
        self.fragments_changed.emit([])

    def _set_details(self, text):
        self.label.setText(" &nbsp;·&nbsp; ".join(t for t in (self._count, text) if t))

    def _on_psm_selected(self, current, previous):
        if not current.isValid() or self.location_path is None:
            self._clear_psm()
            return
        file_name = self.psm_model.rows[self.psm_proxy.mapToSource(current).row()][0]
        psm_path = self.location_path / f"{file_name}.tsv"
        try:
            psm, frags = _load_psm(psm_path)
        except (OSError, ValueError) as e:
            self._clear_psm()
            self._set_details(f"Cannot read {html.escape(str(psm_path))}: {html.escape(str(e))}")
            return

        self._set_details(f"<b>{html.escape(psm.get('pep', '?'))}</b> &nbsp; charge {html.escape(psm.get('z', '?'))} &nbsp; "
                          f"<span style='color:gray'>{html.escape(str(psm_path))}</span>")
        self.field_model.set_rows([(k, v) for k, v in psm.items() if k != "ms2"]
                                  + [("ms2", f"{len(frags)} fragments")])
        self.frag_model.set_rows([(_fragment_label(f), f["type"], f["number"], f["charge"], f["loss"], f["mz"], f["intensity"])
                                  for f in frags])
        self.field_view.resizeColumnToContents(0)
        self.fragments_changed.emit(frags)


class _ChargeBox(QGroupBox):
    """Charges of the selected peptide, with a description on top; emits the selected charge (None when cleared)."""

    charge_selected = Signal(object)

    def __init__(self, title, parent=None):
        super().__init__(title, parent)
        self.label = QLabel()
        self.label.setWordWrap(True)
        self.model = _TableModel(["Charge"])
        self.view, self.proxy = _make_table(self.model)
        self.view.selectionModel().currentRowChanged.connect(self._on_selected)
        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.view)

    def set_charges(self, charges, description):
        """Show the charges and select the first one."""
        self.model.set_rows([(str(z), z) for z in charges])  # shown as text; the charge itself is emitted
        self.label.setText(description)
        if charges:
            self.view.selectRow(0)

    def clear(self, description=""):
        self.model.set_rows([])
        self.label.setText(description)

    def _on_selected(self, current, previous):
        if not current.isValid():
            self.charge_selected.emit(None)
            return
        self.charge_selected.emit(self.model.rows[self.proxy.mapToSource(current).row()][1])


class _ViewerWidget(QWidget):
    """Viewer tab: browse peptides, PSMs and spectra of a specset."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.specset_path = None
        self.location_path = None

        # specset selection
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("Specset directory")
        self.path_edit.returnPressed.connect(self._open_from_edit)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        open_button = QPushButton("Open")
        open_button.clicked.connect(self._open_from_edit)
        path_row = QHBoxLayout()
        path_row.addWidget(QLabel("Specset:"))
        path_row.addWidget(self.path_edit, 1)
        path_row.addWidget(browse)
        path_row.addWidget(open_button)

        # _index.tsv: peptides
        self.pep_model = _TableModel(["Peptide", "Location"])
        self.pep_view, self.pep_proxy = _make_table(self.pep_model)
        self.pep_filter = QLineEdit()
        self.pep_filter.setPlaceholderText("Filter peptides…")
        self.pep_filter.setClearButtonEnabled(True)
        self._filter_timer = QTimer(self)
        self._filter_timer.setSingleShot(True)
        self._filter_timer.setInterval(300)
        self._filter_timer.timeout.connect(self._apply_filter)
        self.pep_filter.textChanged.connect(self._filter_timer.start)
        self.pep_count = QLabel()
        pep_box = QGroupBox("Peptides (_index.tsv)")
        pep_layout = QVBoxLayout(pep_box)
        pep_layout.addWidget(self.pep_filter)
        pep_layout.addWidget(self.pep_view)
        pep_layout.addWidget(self.pep_count)
        self.pep_view.selectionModel().currentRowChanged.connect(self._on_peptide_selected)

        self.spectrum = _SpectrumWidget()

        # charges of the selected peptide (from __index.tsv)
        self.charge_box = _ChargeBox("Charges")
        self.charge_box.charge_selected.connect(self._on_charge_selected)

        # __index.tsv of the selected peptide and the selected PSM
        self.psm_panel = _PsmPanel("PSMs (__index.tsv)")
        self.psm_panel.fragments_changed.connect(self.spectrum.set_fragments)

        # top: peptides and plot; bottom: charges, PSMs and tables
        top = QSplitter(Qt.Orientation.Horizontal)
        top.addWidget(pep_box)
        top.addWidget(self.spectrum)
        top.setSizes([380, 800])

        bottom = QSplitter(Qt.Orientation.Horizontal)
        bottom.addWidget(self.charge_box)
        bottom.addWidget(self.psm_panel)
        bottom.setStretchFactor(1, 1)
        bottom.setSizes([110, 1070])

        main_split = QSplitter(Qt.Orientation.Vertical)
        main_split.addWidget(top)
        main_split.addWidget(bottom)
        main_split.setSizes([400, 350])

        layout = QVBoxLayout(self)
        layout.addLayout(path_row)
        layout.addWidget(main_split, 1)
        self._aligned = False

    def showEvent(self, event):
        super().showEvent(event)
        if not self._aligned:
            self._aligned = True
            QTimer.singleShot(0, lambda: self.psm_panel.align_fragments_to(self.spectrum))

    def _browse(self):
        start = self.path_edit.text() or str(Path.home())
        path = QFileDialog.getExistingDirectory(self, "Select specset directory", start)
        if path:
            self.path_edit.setText(path)
            self.open_specset(path)

    def _open_from_edit(self):
        path = self.path_edit.text().strip()
        if path:
            self.open_specset(path)

    def open_specset(self, path):
        specset_path = Path(path)
        index_file = specset_path / "_index.tsv"
        if not index_file.is_file():
            QMessageBox.warning(self, "Not a specset", f"No _index.tsv found in:\n{path}")
            return

        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            rows = _read_index(specset_path)
        finally:
            QApplication.restoreOverrideCursor()

        self.specset_path = specset_path
        self.path_edit.setText(str(specset_path))
        self.pep_view.sortByColumn(-1, Qt.SortOrder.AscendingOrder)  # keep file order
        self.pep_model.set_rows(rows)
        self._clear_charges()
        self._update_pep_count()
        self.window().statusBar().showMessage(f"Opened {specset_path}")
        if rows:
            self.pep_view.selectRow(0)

    def _apply_filter(self):
        self.pep_proxy.setFilterFixedString(self.pep_filter.text().strip())
        self._update_pep_count()

    def _update_pep_count(self):
        total = self.pep_model.rowCount()
        shown = self.pep_proxy.rowCount()
        self.pep_count.setText(f"{total} peptides" if shown == total else f"{shown} of {total} peptides")

    def _clear_charges(self, description=""):
        self.location_path = None
        self.charge_box.clear(description)
        self.psm_panel.clear()

    def _on_peptide_selected(self, current, previous):
        self._clear_charges()
        if not current.isValid() or self.specset_path is None:
            return
        location_path = self.specset_path / self.pep_model.rows[self.pep_proxy.mapToSource(current).row()][1]
        try:
            charges = sorted({z for _, z in _read_location_index(location_path)})
        except OSError as e:
            self._clear_charges(f"Cannot read __index.tsv: {e}")
            return
        self.location_path = location_path
        self.charge_box.set_charges(charges, f"{len(charges)} charges")

    def _on_charge_selected(self, charge):
        if charge is None or self.location_path is None:
            self.psm_panel.clear()
            return
        self.psm_panel.set_location(self.location_path, charge)


class _CompareViewerWidget(QWidget):
    """Compare tab: peptides shared by two specsets, their common charges, and a mirror plot of one PSM from each."""

    _SIDES = ("A", "B")

    def __init__(self, parent=None):
        super().__init__(parent)
        self.specset_paths = [None, None]
        self.locations = [{}, {}]            # peptide -> location, per specset
        self.current_pep = None
        self.frags = [[], []]
        self._totals = ""

        # specset selection
        self.path_edits = [_PathEdit(DIR), _PathEdit(DIR)]
        paths = QFormLayout()
        for side, edit in zip(self._SIDES, self.path_edits):
            edit.edit.setPlaceholderText(f"Specset {side} directory")
            edit.edit.returnPressed.connect(self._compare)
            paths.addRow(f"Specset {side}:", edit)
        compare = QPushButton("Compare")
        compare.clicked.connect(self._compare)
        path_row = QHBoxLayout()
        path_row.addLayout(paths, 1)
        path_row.addWidget(compare, 0, Qt.AlignmentFlag.AlignBottom)

        # overlapping peptides
        self.pep_model = _TableModel(["Peptide"])
        self.pep_view, self.pep_proxy = _make_table(self.pep_model)
        self.pep_filter = QLineEdit()
        self.pep_filter.setPlaceholderText("Filter peptides…")
        self.pep_filter.setClearButtonEnabled(True)
        self._filter_timer = QTimer(self)
        self._filter_timer.setSingleShot(True)
        self._filter_timer.setInterval(300)
        self._filter_timer.timeout.connect(self._apply_filter)
        self.pep_filter.textChanged.connect(self._filter_timer.start)
        self.pep_count = QLabel()
        pep_box = QGroupBox("Overlapping peptides (_index.tsv)")
        pep_layout = QVBoxLayout(pep_box)
        pep_layout.addWidget(self.pep_filter)
        pep_layout.addWidget(self.pep_view)
        pep_layout.addWidget(self.pep_count)
        self.pep_view.selectionModel().currentRowChanged.connect(self._on_peptide_selected)

        # charges of the selected peptide present in both specsets (from __index.tsv)
        self.charge_box = _ChargeBox("Common charges")
        self.charge_box.charge_selected.connect(self._on_charge_selected)

        # mirror plot
        self.spectrum = _SpectrumWidget()

        top = QSplitter(Qt.Orientation.Horizontal)
        top.addWidget(pep_box)
        top.addWidget(self.spectrum)
        top.setSizes([460, 720])

        # per specset: __index.tsv of the selected peptide and the selected PSM
        self.psm_panels = []
        panels = QSplitter(Qt.Orientation.Vertical)
        for i, side in enumerate(self._SIDES):
            panel = _PsmPanel(f"Specset {side}: PSMs (__index.tsv)")
            panel.fragments_changed.connect(lambda frags, i=i: self._on_fragments_changed(i, frags))
            self.psm_panels.append(panel)
            panels.addWidget(panel)

        # bottom: common charges on the left, spanning both PSM panels
        bottom = QSplitter(Qt.Orientation.Horizontal)
        bottom.addWidget(self.charge_box)
        bottom.addWidget(panels)
        bottom.setStretchFactor(1, 1)
        bottom.setSizes([110, 1070])

        main_split = QSplitter(Qt.Orientation.Vertical)
        main_split.addWidget(top)
        main_split.addWidget(bottom)
        main_split.setSizes([380, 460])

        layout = QVBoxLayout(self)
        layout.addLayout(path_row)
        layout.addWidget(main_split, 1)
        self._aligned = False

    def showEvent(self, event):
        super().showEvent(event)
        if not self._aligned:
            self._aligned = True
            QTimer.singleShot(0, self._align_panels)

    def _align_panels(self):
        for panel in self.psm_panels:
            panel.align_fragments_to(self.spectrum)

    def _compare(self):
        texts = [edit.text() for edit in self.path_edits]
        if not all(texts):
            QMessageBox.warning(self, "Missing specset", "Please select both specsets.")
            return
        paths = [Path(t) for t in texts]
        for side, path in zip(self._SIDES, paths):
            if not (path / "_index.tsv").is_file():
                QMessageBox.warning(self, "Not a specset", f"Specset {side}: no _index.tsv found in:\n{path}")
                return

        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            indexes = [_read_index(path) for path in paths]
        finally:
            QApplication.restoreOverrideCursor()

        self.specset_paths = paths
        self.locations = [dict(index) for index in indexes]
        overlap = [(pep,) for pep, _ in indexes[0] if pep in self.locations[1]]
        self._totals = f"A: {len(self.locations[0])}, B: {len(self.locations[1])}"

        self.pep_view.sortByColumn(-1, Qt.SortOrder.AscendingOrder)  # keep the order of specset A
        self.pep_model.set_rows(overlap)
        self._update_pep_count()
        self.window().statusBar().showMessage(f"Compared {paths[0]} and {paths[1]}: {len(overlap)} overlapping peptides")
        if overlap:
            self.pep_view.selectRow(0)
        else:
            self._clear_charges()

    def _apply_filter(self):
        self.pep_proxy.setFilterFixedString(self.pep_filter.text().strip())
        self._update_pep_count()

    def _update_pep_count(self):
        total = self.pep_model.rowCount()
        shown = self.pep_proxy.rowCount()
        count = f"{total}" if shown == total else f"{shown} of {total}"
        self.pep_count.setText(f"{count} overlapping peptides ({self._totals})")

    def _clear_psms(self):
        for panel in self.psm_panels:
            panel.clear()

    def _clear_charges(self, description=""):
        self.current_pep = None
        self.charge_box.clear(description)
        self._clear_psms()

    def _on_peptide_selected(self, current, previous):
        self._clear_charges()
        if not current.isValid() or None in self.specset_paths:
            return
        pep = self.pep_model.rows[self.pep_proxy.mapToSource(current).row()][0]
        try:
            charges = [{z for _, z in _read_location_index(self.specset_paths[i] / self.locations[i][pep])} for i in range(2)]
        except OSError as e:
            self._clear_charges(f"Cannot read __index.tsv: {e}")
            return
        common = sorted(charges[0] & charges[1])
        self.current_pep = pep
        self.charge_box.set_charges(common, f"{len(common)} common")

    def _on_charge_selected(self, charge):
        if charge is None or self.current_pep is None:
            self._clear_psms()
            return
        for i, panel in enumerate(self.psm_panels):
            panel.set_location(self.specset_paths[i] / self.locations[i][self.current_pep], charge)

    def _on_fragments_changed(self, i, frags):
        self.frags[i] = frags
        self.spectrum.set_fragments(self.frags[0], self.frags[1])


class _ToolsWidget(QWidget):
    """Tools tab: run specset commands."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.thread = None
        self.worker = None
        self._carriage_return = False

        self.command_list = QListWidget()
        self.pages = QStackedWidget()
        for command in COMMANDS:
            page = _CommandPage(command)
            page.run_requested.connect(self._run)
            self.command_list.addItem(command["title"])
            self.pages.addWidget(page)
        self.command_list.currentRowChanged.connect(self.pages.setCurrentIndex)
        self.command_list.setCurrentRow(0)
        self.command_list.setFixedWidth(220)

        top = QWidget()
        top_layout = QHBoxLayout(top)
        top_layout.addWidget(self.command_list)
        top_layout.addWidget(self.pages, 1)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(QFont("Monospace"))
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        clear = QPushButton("Clear log")
        clear.clicked.connect(self.log.clear)

        bottom = QWidget()
        bottom_layout = QVBoxLayout(bottom)
        bottom_layout.addWidget(QLabel("Log"))
        bottom_layout.addWidget(self.log)
        row = QHBoxLayout()
        row.addWidget(self.progress, 1)
        row.addWidget(clear)
        bottom_layout.addLayout(row)

        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(top)
        splitter.addWidget(bottom)
        splitter.setSizes([420, 280])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)

        self.stream = _Stream()
        self.stream.text_written.connect(self._append_log)

    def _run(self, func_name, kwargs):
        if self.thread is not None:
            return

        out_dir, is_parent = self.pages.currentWidget().output_dir(kwargs)
        if kwargs.get("overwrite") and out_dir and Path(out_dir).is_dir() and any(Path(out_dir).iterdir()):
            if is_parent:
                message = f"{out_dir}\n\nis not empty. Existing specsets in it with the same name as an output specset will be deleted. Continue?"
            else:
                message = f"{out_dir}\n\nexists and is not empty. All of its contents will be deleted. Continue?"
            answer = QMessageBox.question(self, "Overwrite output?", message)
            if answer != QMessageBox.StandardButton.Yes:
                return

        self._set_running(True)
        self._append_log(f"\n=== {func_name}({', '.join(f'{k}={v!r}' for k, v in kwargs.items())})\n")

        self.thread = QThread()
        self.worker = _Worker(func_name, kwargs, self.stream)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self._on_finished)
        self.worker.finished.connect(self.thread.quit)
        self.thread.finished.connect(self._on_thread_finished)
        self.thread.start()

    def _on_finished(self, ok, error):
        if ok:
            self.progress.setRange(0, 100)
            self.progress.setValue(100)
            self._append_log("\n=== Finished.\n")
            self._status("Finished")
        else:
            self.progress.setRange(0, 100)
            self.progress.setValue(0)
            self._append_log(f"\n=== Failed.\n{error}")
            self._status("Failed")
            QMessageBox.critical(self, "Error", error.strip().splitlines()[-1])

    def _on_thread_finished(self):
        self.thread.deleteLater()
        self.worker.deleteLater()
        self.thread = None
        self.worker = None
        self._set_running(False)

    def _set_running(self, running):
        for i in range(self.pages.count()):
            self.pages.widget(i).run_button.setEnabled(not running)
        if running:
            self.progress.setRange(0, 0)  # busy until a percentage is printed
            self._status("Running…")

    def _append_log(self, text):
        """Append text to the log; a carriage return overwrites the current line, like a terminal."""
        cursor = self.log.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        for part in re.split(r"(\r|\n)", text):
            if part == "\n":
                cursor.insertText("\n")
                self._carriage_return = False
            elif part == "\r":
                self._carriage_return = True
            elif part:
                if self._carriage_return:
                    cursor.movePosition(QTextCursor.MoveOperation.StartOfBlock, QTextCursor.MoveMode.KeepAnchor)
                    cursor.removeSelectedText()
                    self._carriage_return = False
                cursor.insertText(part)
        self.log.setTextCursor(cursor)
        self.log.ensureCursorVisible()

        if self.thread is not None:
            percents = re.findall(r"(\d+)%", text)
            if percents:
                self.progress.setRange(0, 100)
                self.progress.setValue(min(int(percents[-1]), 100))

    def is_running(self):
        return self.thread is not None

    def _status(self, message):
        self.window().statusBar().showMessage(message)


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("specset")
        self.resize(1200, 800)

        self.viewer = _ViewerWidget()
        self.compare_viewer = _CompareViewerWidget()
        self.tools = _ToolsWidget()
        tabs = QTabWidget()
        tabs.addTab(self.viewer, "Viewer")
        tabs.addTab(self.compare_viewer, "Compare")
        tabs.addTab(self.tools, "Tools")
        self.setCentralWidget(tabs)
        self.statusBar().showMessage("Ready")

    def closeEvent(self, event):
        if self.tools.is_running():
            QMessageBox.information(self, "Task running", "A task is still running. Please wait until it finishes.")
            event.ignore()
        else:
            event.accept()


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
