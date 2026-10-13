"""PySide6 graphical interface for specset."""

import contextlib
import html
import math
import re
import sys
import traceback
from pathlib import Path

from PySide6.QtCore import (
    QAbstractTableModel, QMarginsF, QModelIndex, QObject, QPoint, QRect, QRectF, QSize, QSortFilterProxyModel, Qt, QThread,
    QTimer, Signal,
)
from PySide6.QtGui import QColor, QFont, QIcon, QPageLayout, QPageSize, QPainter, QPalette, QPdfWriter, QPen, QTextCursor
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QLineEdit, QListWidget, QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar, QProgressDialog, QPushButton,
    QSizePolicy, QSpinBox, QSplitter, QStackedWidget, QTabWidget, QTableView, QVBoxLayout, QWidget,
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
KEY_VALUES = "key_values"  # one KEY=VALUE per line -> dict ("multi": values as lists, "number": values as floats)
TEXT = "text"              # single line of text -> str ("" if empty)
LINES = "lines"            # one value per line -> list of str

_INT_MAX = 2**31 - 1

# Each command: utils function name, title, description, and fields as
# (keyword argument, label, kind, options). Layout options: "panel": "left", "middle" or "right"
# puts a field in a column of the side-by-side Filters panel, label above the box; "lines": n
# makes a text box n lines high; "stretch": True makes a text box fill the height of its panel
# column. Optional "prepare": name of a
# function in this module that turns the form values into the keyword arguments of the utils
# function.
COMMANDS = [
    {
        "func": "maxquant_to_specset",
        "title": "MaxQuant → specset",
        "help": "Convert a MaxQuant txt directory (msms.txt and evidence.txt) to a specset.",
        "fields": [
            ("txt_dir", "MaxQuant txt directory", DIR, {}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("additional_columns", "Additional columns", KEY_VALUES, {"placeholder": "msms.txt column=specset name, one per line\ne.g. Gene Names=genes",
                                                                       "default": "Score=score\nProteins=proteins\nGene Names=genes"}),
            ("max_files_in_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
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
        ],
    },
    {
        "func": "diannreport_to_specset",
        "title": "DIA-NN report → specset",
        "help": "Convert a DIA-NN parquet report to a specset.",
        "fields": [
            ("diannreport", "DIA-NN report", FILE, {"filter": "Parquet files (*.parquet);;All files (*)"}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("column", "Select rows on column", TEXT, {"panel": "left", "placeholder": "Optional, e.g. Run",
                                                      "tooltip": "Report column to select rows on. Leave empty to import all rows."}),
            ("values", "Column values", LINES, {"panel": "left", "stretch": True, "placeholder": "One value per line, e.g.\nplasma1\nplasma2", "required_with": "column",
                                               "tooltip": "Import only the rows whose value in the column above is one of these."}),
            ("string_equal", "STRING keep if equal", KEY_VALUES, {
                "panel": "middle", "lines": 6,
                "multi": True,
                "placeholder": "column=value, one per line (optional)\ne.g. Run=plasma1\nRun=plasma2",
                "tooltip": "Keep only rows whose column equals one of the values, compared as text.",
            }),
            ("string_not_equal", "STRING drop if equal", KEY_VALUES, {
                "panel": "middle", "lines": 6,
                "multi": True,
                "placeholder": "column=value, one per line (optional)",
                "tooltip": "Drop rows whose column equals one of the values, compared as text.",
            }),
            ("double_less", "NUMERIC keep if less than or equal to", KEY_VALUES, {
                "panel": "right", "lines": 6,
                "number": True,
                "default": "Q.Value=0.01\nPG.Q.Value=0.05\nLib.Peptidoform.Q.Value=0.01\n"
                           "Global.Q.Value=0.01\nGlobal.PG.Q.Value=0.01",
                "placeholder": "column=number, one per line (optional)\ne.g. Q.Value=0.01",
                "tooltip": "Keep only rows whose column, read as a number, is less than or equal to the number.",
            }),
            ("double_greater", "NUMERIC keep if greater than or equal to", KEY_VALUES, {
                "panel": "right", "lines": 6,
                "number": True,
                "placeholder": "column=number, one per line (optional)",
                "tooltip": "Keep only rows whose column, read as a number, is greater than or equal to the number.",
            }),
            ("max_no_files_per_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
        ],
        "prepare": "_column_values_to_string_equal",
    },
    {
        "func": "psmlist_to_specset",
        "title": "PSM list → specset",
        "help": "Convert a list of files containing stripped MaxQuant PSMs to a specset.",
        "fields": [
            ("psmlist_file", "PSM list file", FILE, {"filter": "Text files (*.txt *.tsv);;All files (*)"}),
            ("out_dir", "Output specset", OUT_DIR, {}),
            ("max_files_in_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
        ],
    },
    {
        "func": "specset_to_speclib",
        "title": "Specset → spectral library",
        "help": "Convert a specset to a spectral library TSV.",
        "fields": [
            ("specset_dir", "Specset", DIR, {}),
            ("lib_name", "Output library", OUT_FILE, {"filter": "TSV files (*.tsv);;All files (*)"}),
            ("additional_columns", "Additional columns", KEY_VALUES, {"placeholder": "specset field=library column, one per line\ne.g. proteins=Protein.Ids",
                                                                       "default": "proteins=Protein.Ids\ngenes=Genes"}),
            ("consensus_method", "Consensus method", CHOICE, {"choices": ["", "mean", "median"], "default": "score", "editable": True,
                                                              "placeholder": "Specify a specset field for max value, e.g. score, or leave empty for max number of fragments",
                                                              "tooltip": "Empty: PSM with most fragments, mean, median, or a specset field for selection with the max value."}),
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
            ("mode", "Mode", CHOICE, {"choices": ["both", "ignore"], "default": "ignore",
                                      "tooltip": "For peptides in both: 'both' keeps the PSMs from A and B, 'ignore' skips the ones from B."}),
            ("max_files_in_directory", "Files per folder", INT, {"default": 5000, "min": 1}),
            ("flush_threshold", "Flush threshold", INT, {"default": 1000, "min": 1}),
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


# Last folder used in a file dialog, remembered for this session only
_last_dir = ""


def _start_dir():
    """Folder to start a file dialog in: the last folder used in this session, or the home folder."""
    return _last_dir if _last_dir and Path(_last_dir).is_dir() else str(Path.home())


def _remember_dir(path):
    """Remember the folder containing a picked file or directory, for the next file dialog."""
    global _last_dir
    _last_dir = str(Path(path).parent)


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
        start = self.edit.text() or _start_dir()
        if self.kind in (DIR, OUT_DIR):
            path = QFileDialog.getExistingDirectory(self, "Select directory", start)
        elif self.kind == FILE:
            path, _ = QFileDialog.getOpenFileName(self, "Select file", start, self.file_filter)
        else:
            path, _ = QFileDialog.getSaveFileName(self, "Save as", start, self.file_filter)
        if path:
            self.edit.setText(path)
            _remember_dir(path)

    def text(self):
        return self.edit.text().strip()


def _set_box_height(box, opts):
    """Height of a text box: fills its panel column ("stretch"), n lines ("lines"), or 70 pixels."""
    if opts.get("stretch"):
        # ignore the box's own preferred height, so the other columns set the panel height
        box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        box.setMinimumHeight(70)
    elif "lines" in opts:
        lines_height = opts["lines"] * box.fontMetrics().lineSpacing()
        box.setFixedHeight(math.ceil(lines_height + 2 * box.document().documentMargin() + 2 * box.frameWidth()))
    else:
        box.setFixedHeight(70)


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
        panel_names = [p for p in ("left", "middle", "right") if any(o.get("panel") == p for *_, o in command["fields"])]
        panels = None  # column layouts by panel name, created at the first field with a panel
        stretched = set()  # panels with a box that fills the column height
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
                if w.isEditable() and "placeholder" in opts:
                    w.lineEdit().setPlaceholderText(opts["placeholder"])  # grayed-out tip while the box is empty
                if "default" in opts:
                    w.setCurrentText(opts["default"])
            elif kind == BOOL:
                w = QCheckBox()
                w.setChecked(opts.get("default", False))
            elif kind == KEY_VALUES:
                w = QPlainTextEdit()
                w.setPlaceholderText(opts.get("placeholder", "KEY=VALUE, one per line"))
                w.setPlainText(opts.get("default", ""))
                _set_box_height(w, opts)
            elif kind == TEXT:
                w = QLineEdit()
                w.setPlaceholderText(opts.get("placeholder", ""))
            elif kind == LINES:
                w = QPlainTextEdit()
                w.setPlaceholderText(opts.get("placeholder", "One value per line"))
                _set_box_height(w, opts)
            else:
                raise ValueError(f"unknown field kind: {kind}")
            if "tooltip" in opts:
                w.setToolTip(opts["tooltip"])
            self.widgets[name] = (kind, opts, w)
            if "panel" in opts:
                if panels is None:
                    panels = {p: QVBoxLayout() for p in panel_names}
                    columns = QHBoxLayout()
                    for p in panel_names:
                        columns.addLayout(panels[p], 1)  # equal widths
                    form.addRow("Filters", columns)  # in the field column, aligned with the other boxes
                panels[opts["panel"]].addWidget(QLabel(label))
                panels[opts["panel"]].addWidget(w, 1 if opts.get("stretch") else 0)
                if opts.get("stretch"):
                    stretched.add(opts["panel"])
            else:
                form.addRow(label, w)
        if panels is not None:
            for p, panel in panels.items():
                if p not in stretched:
                    panel.addStretch()  # keep the columns top-aligned
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
                    key, value = (part.strip() for part in line.split("=", 1))
                    if opts.get("number"):
                        try:
                            value = float(value)
                        except ValueError:
                            raise ValueError(f"'{self._label(name)}': expected KEY=NUMBER, got '{line}'")
                    if opts.get("multi"):  # repeated keys collect their values in a list
                        d.setdefault(key, []).append(value)
                    else:
                        d[key] = value
                kwargs[name] = d
            elif kind == TEXT:
                kwargs[name] = w.text().strip()
            elif kind == LINES:
                kwargs[name] = [line.strip() for line in w.toPlainText().splitlines() if line.strip()]

        # fields that must be filled in when another field is, e.g. values with column
        for name, (kind, opts, w) in self.widgets.items():
            other = opts.get("required_with")
            if other and kwargs.get(other) and not kwargs.get(name):
                raise ValueError(f"Please fill in '{self._label(name)}' when '{self._label(other)}' is given.")
        return kwargs

    def _label(self, name):
        return next(label for n, label, _, _ in self.command["fields"] if n == name)

    def _on_run(self):
        try:
            kwargs = self.values()
        except ValueError as e:
            QMessageBox.warning(self, "Invalid input", str(e))
            return
        prepare = self.command.get("prepare")
        if prepare:
            kwargs = globals()[prepare](kwargs)
        self.run_requested.emit(self.command["func"], kwargs)


def _column_values_to_string_equal(kwargs):
    """Turn the GUI's column and values fields into a string_equal filter; both filters must hold if both are given."""
    column = kwargs.pop("column")
    values = kwargs.pop("values")
    if column:
        string_equal = kwargs["string_equal"]
        string_equal[column] = [v for v in string_equal[column] if v in values] if column in string_equal else values
    return kwargs


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


def _iter_index(specset_path):
    """Yield (peptide, location) from _index.tsv of a specset."""
    with open(Path(specset_path) / "_index.tsv") as f:
        for line in f:
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                yield parts[0], parts[1].strip()


def _read_index(specset_path):
    """Read _index.tsv of a specset into a list of (peptide, location)."""
    return list(_iter_index(specset_path))


def _read_location_index(location_path):
    """Read __index.tsv of a location into a list of (file name, charge)."""
    rows = []
    with open(Path(location_path) / "__index.tsv") as f:
        for line in f:
            parts = line.split()
            if parts:
                rows.append((parts[0], _to_int(parts[1]) if len(parts) > 1 else None))
    return rows


def _read_psm(psm_path):
    """Read a PSM file; return (fields, fragments with m/z). Raises OSError or ValueError."""
    psm = _read_psm_file(psm_path)
    frags = _parse_ms2(psm.get("ms2", ""))
    _fragment_mz(psm.get("pep", ""), frags)
    return psm, frags


def _load_psm(psm_path):
    """Like _read_psm, with a wait cursor (the first m/z calculation can take a while)."""
    QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
    try:
        return _read_psm(psm_path)
    finally:
        QApplication.restoreOverrideCursor()


_FIELD_HEADERS = ["Field", "Value"]
_FRAGMENT_HEADERS = ["Fragment", "Type", "Number", "Charge", "Loss", "m/z", "Intensity"]


def _field_rows(psm, frags):
    return [(k, v) for k, v in psm.items() if k != "ms2"] + [("ms2", f"{len(frags)} fragments")]


def _fragment_rows(frags):
    return [(_fragment_label(f), f["type"], f["number"], f["charge"], f["loss"], f["mz"], f["intensity"]) for f in frags]


def _format_cell(value):
    """Text of a table cell, as shown in the tables and exported to PDF."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


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
            return _format_cell(value)
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
        pal = self.palette()
        self.paint(painter, self.width(), self.height(), pal.color(QPalette.ColorRole.Base), pal.color(QPalette.ColorRole.Text))

    def paint(self, painter, width, height, background, text_color):
        """Draw the plot in the rectangle (0, 0, width, height) of painter; also used for PDF and SVG export."""
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(QRectF(0, 0, width, height), background)

        fm = painter.fontMetrics()
        left, right, top, bottom = fm.horizontalAdvance("100%") + 12, 16, fm.height() + 12, fm.height() + 14
        w = width - left - right
        h = height - top - bottom
        if w <= 10 or h <= 10:
            return

        mirrored = self.mirror is not None
        sides = [(1, self.frags)] + ([(-1, self.mirror)] if mirrored else [])
        baseline = top + h / 2 if mirrored else top + h
        half = h / 2 if mirrored else h

        all_frags = [f for _, frags in sides for f in frags]
        painter.setPen(QPen(text_color, 1))
        if not all_frags:
            painter.drawText(QRectF(0, 0, width, height), Qt.AlignmentFlag.AlignCenter, "No fragments")
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
        painter.drawText(left + w - fm.horizontalAdvance(x_label), height - 2, x_label)

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


def _ask_export_path(parent, caption, default_name, suffix, file_filter):
    """Ask for a file to export to; return its path with the suffix, or None if cancelled."""
    start = Path(_start_dir()) / default_name
    path, _ = QFileDialog.getSaveFileName(parent, caption, str(start), file_filter)
    if not path:
        return None
    if not path.lower().endswith(suffix):
        path += suffix
    _remember_dir(path)
    return path


def _export_name(psm):
    """Default export file name (without suffix) for a PSM, e.g. PEPTIDEK_z2."""
    name = re.sub(r'[\\/:*?"<>|\s()]+', "_", f"{psm.get('pep', 'spectrum')}_z{psm.get('z', '')}").strip("_")
    return name or "spectrum"


def _export_svg(spectrum, path, title):
    """Write the plot as an SVG, at the size it has on screen."""
    from PySide6.QtSvg import QSvgGenerator

    width, height = max(spectrum.width(), 600), max(spectrum.height(), 300)
    generator = QSvgGenerator()
    generator.setFileName(path)
    generator.setSize(QSize(width, height))
    generator.setViewBox(QRect(0, 0, width, height))
    generator.setTitle(title)
    painter = QPainter(generator)
    try:
        spectrum.paint(painter, width, height, QColor("white"), QColor("black"))
    finally:
        painter.end()


class _PdfReport:
    """Write text, plots and tables to an A4 PDF, starting new pages as needed."""

    _FLAGS = int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop) | int(Qt.TextFlag.TextWrapAnywhere)

    def __init__(self, path):
        self.writer = QPdfWriter(path)
        self.writer.setResolution(96)  # so that sizes match the screen
        self.writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
        self.writer.setPageMargins(QMarginsF(15, 15, 15, 15), QPageLayout.Unit.Millimeter)
        self.painter = QPainter(self.writer)
        if not self.painter.isActive():
            raise OSError(f"Cannot write {path}")
        self.width = self.writer.width()
        self.height = self.writer.height()
        self.y = 0

    def _font(self, size, bold=False):
        font = QFont()
        font.setPointSizeF(size)
        font.setBold(bold)
        return font

    def _text_height(self, text, font, width):
        self.painter.setFont(font)
        return self.painter.fontMetrics().boundingRect(QRect(0, 0, int(width), 100000), self._FLAGS, text).height()

    def _new_page_if_needed(self, height):
        if self.y > 0 and self.y + height > self.height:
            self.writer.newPage()
            self.y = 0

    def text(self, text, size=10, bold=False, color="black"):
        font = self._font(size, bold)
        height = self._text_height(text, font, self.width)
        self._new_page_if_needed(height)
        self.painter.setFont(font)
        self.painter.setPen(QColor(color))
        self.painter.drawText(QRectF(0, self.y, self.width, height), self._FLAGS, text)
        self.y += height + 4

    def plot(self, spectrum, aspect=0.5):
        """Draw a _SpectrumWidget plot across the page width, as vector graphics."""
        height = int(self.width * aspect)
        self._new_page_if_needed(height)
        self.painter.save()
        self.painter.translate(0, self.y)
        self.painter.setFont(self._font(8))
        spectrum.paint(self.painter, self.width, height, QColor("white"), QColor("black"))
        self.painter.restore()
        self.y += height + 12

    def table(self, title, headers, rows, weights):
        """Draw a table with a title; rows are never split, and the header is repeated on a new page."""
        title_font, font, bold = self._font(10, True), self._font(8), self._font(8, True)
        widths = [self.width * w / sum(weights) for w in weights]
        pad = 3
        header = [str(h) for h in headers]
        cells = [[_format_cell(v) for v in row[:len(headers)]] for row in rows]

        def row_height(row, row_font):
            return max(self._text_height(c, row_font, w - 2 * pad) for c, w in zip(row, widths)) + 2 * pad

        def draw_row(row, row_font, height, shade):
            if shade:
                self.painter.fillRect(QRectF(0, self.y, self.width, height), QColor(235, 235, 235))
            self.painter.setFont(row_font)
            self.painter.setPen(QColor("black"))
            x = 0
            for c, w in zip(row, widths):
                self.painter.drawText(QRectF(x + pad, self.y + pad, w - 2 * pad, height - 2 * pad), self._FLAGS, c)
                x += w
            self.painter.setPen(QPen(QColor(200, 200, 200), 0.5))
            self.painter.drawLine(QRectF(0, self.y, self.width, height).bottomLeft(), QRectF(0, self.y, self.width, height).bottomRight())
            self.y += height

        header_height = row_height(header, bold)
        title_height = self._text_height(title, title_font, self.width) + 4
        first_height = row_height(cells[0], font) if cells else 0
        self._new_page_if_needed(title_height + header_height + first_height)  # keep the title with the first row
        self.text(title, 10, bold=True)
        draw_row(header, bold, header_height, True)
        for row in cells:
            height = row_height(row, font)
            if self.y + height > self.height:
                self.writer.newPage()
                self.y = 0
                draw_row(header, bold, header_height, True)
            draw_row(row, font, height, False)
        if not cells:
            self.text("(none)", 8, color="gray")
        self.y += 12

    def new_page(self):
        """Continue on a new page, unless the current page is still empty."""
        if self.y > 0:
            self.writer.newPage()
            self.y = 0

    def fields_table(self, title, rows):
        self.table(title, _FIELD_HEADERS, rows, [1, 3])

    def fragments_table(self, title, rows):
        self.table(title, _FRAGMENT_HEADERS, rows, [1.3, 0.8, 1, 1, 1, 1.3, 1.3])

    def close(self):
        self.painter.end()


def _write_pages(parent, path, jobs, write_page, title):
    """
    Write each job to a new page of a PDF with write_page(report, job), showing a progress
    dialog with Cancel. Return True if completed; on cancel, the partial PDF is deleted.
    """
    progress = QProgressDialog("Writing the PDF…", "Cancel", 0, len(jobs), parent)
    progress.setWindowTitle(title)
    progress.setWindowModality(Qt.WindowModality.WindowModal)
    progress.setMinimumDuration(0)
    completed = False
    report = _PdfReport(path)
    try:
        for n, job in enumerate(jobs):
            progress.setValue(n)
            if progress.wasCanceled():
                break
            report.new_page()
            write_page(report, job)
        else:
            completed = True
    finally:
        report.close()
        progress.close()
    if not completed:
        Path(path).unlink(missing_ok=True)
    return completed


def _collect_jobs(parent, items, collect, title):
    """Call collect(item) for each item with a progress dialog; return the concatenated jobs, or None if cancelled."""
    progress = QProgressDialog("Collecting PSMs…", "Cancel", 0, len(items), parent)
    progress.setWindowTitle(title)
    progress.setWindowModality(Qt.WindowModality.WindowModal)
    progress.setMinimumDuration(500)
    jobs = []
    try:
        for n, item in enumerate(items):
            if n % 100 == 0:
                progress.setValue(n)
                if progress.wasCanceled():
                    return None
            jobs.extend(collect(item))
    finally:
        progress.close()
    return jobs


def _shown_rows(proxy, model):
    """Source rows of a table model in the order shown, filtered and sorted, by its proxy."""
    return [model.rows[proxy.mapToSource(proxy.index(r, 0)).row()] for r in range(proxy.rowCount())]


class _PsmPanel(QGroupBox):
    """PSMs of one peptide (__index.tsv) aligned with the Fields and Fragments tables of the selected PSM."""

    fragments_changed = Signal(list)

    def __init__(self, title, parent=None):
        super().__init__(title, parent)
        self.location_path = None
        self._count = ""
        self.psm = None          # fields of the selected PSM
        self.psm_path = None

        self.label = QLabel()
        self.label.setTextFormat(Qt.TextFormat.RichText)
        self.label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.label.setWordWrap(True)
        self.psm_model = _TableModel(["File"])  # rows are (file name, charge); charge is used for filtering
        self.psm_view, self.psm_proxy = _make_table(self.psm_model)
        self.field_model = _TableModel(_FIELD_HEADERS)
        self.field_view, _ = _make_table(self.field_model, sortable=False)
        self.frag_model = _TableModel(_FRAGMENT_HEADERS)
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
        self.psm = None
        self.psm_path = None
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
        self.field_model.set_rows(_field_rows(psm, frags))
        self.frag_model.set_rows(_fragment_rows(frags))
        self.field_view.resizeColumnToContents(0)
        self.psm = psm
        self.psm_path = psm_path
        self.fragments_changed.emit(frags)

    def add_tables_to(self, report, title_prefix=""):
        """Add the Fields and Fragments tables of the selected PSM to a _PdfReport."""
        report.fields_table(f"{title_prefix}Fields", self.field_model.rows)
        report.fragments_table(f"{title_prefix}Fragments", self.frag_model.rows)


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
        export_pdf = QPushButton("Export PDF…")
        export_pdf.setToolTip("Save the plot and the tables of the selected PSM as a PDF")
        export_pdf.clicked.connect(self._on_export_pdf)
        export_svg = QPushButton("Export SVG…")
        export_svg.setToolTip("Save the plot as an SVG")
        export_svg.clicked.connect(self._on_export_svg)
        export_all = QPushButton("Export all PDF…")
        export_all.setToolTip("Save every PSM of the peptides currently listed as a PDF, one PSM per page")
        export_all.clicked.connect(self._on_export_all_pdf)
        path_row.addWidget(QLabel("Specset:"))
        path_row.addWidget(self.path_edit, 1)
        path_row.addWidget(browse)
        path_row.addWidget(open_button)
        path_row.addSpacing(16)
        path_row.addWidget(export_pdf)
        path_row.addWidget(export_svg)
        path_row.addWidget(export_all)

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
        start = self.path_edit.text() or _start_dir()
        path = QFileDialog.getExistingDirectory(self, "Select specset directory", start)
        if path:
            _remember_dir(path)
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

    def _on_export_pdf(self):
        psm = self.psm_panel.psm
        if psm is None:
            QMessageBox.information(self, "Export PDF", "Select a PSM first.")
            return
        path = _ask_export_path(self, "Export PDF", _export_name(psm) + ".pdf", ".pdf", "PDF files (*.pdf)")
        if not path:
            return
        try:
            report = _PdfReport(path)
            try:
                report.text(f"{psm.get('pep', '?')}    charge {psm.get('z', '?')}", 13, bold=True)
                report.plot(self.spectrum)
                self.psm_panel.add_tables_to(report)
            finally:
                report.close()
        except Exception as e:
            QMessageBox.critical(self, "Export PDF", f"Export failed: {e}")
            return
        self.window().statusBar().showMessage(f"Exported {path}")

    def _on_export_svg(self):
        psm = self.psm_panel.psm
        if psm is None:
            QMessageBox.information(self, "Export SVG", "Select a PSM first.")
            return
        path = _ask_export_path(self, "Export SVG", _export_name(psm) + ".svg", ".svg", "SVG files (*.svg)")
        if not path:
            return
        try:
            _export_svg(self.spectrum, path, f"{psm.get('pep', '?')} charge {psm.get('z', '?')}")
        except Exception as e:
            QMessageBox.critical(self, "Export SVG", f"Export failed: {e}")
            return
        self.window().statusBar().showMessage(f"Exported {path}")

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

    def _on_export_all_pdf(self):
        title = "Export all PDF"
        peptides = _shown_rows(self.pep_proxy, self.pep_model)
        if self.specset_path is None or not peptides:
            QMessageBox.information(self, title, "Open a specset first; the peptides currently listed are exported.")
            return

        def collect(row):
            pep, location = row
            location_path = self.specset_path / location
            try:
                return [(pep, location_path / f"{file_name}.tsv", file_name) for file_name, _ in _read_location_index(location_path)]
            except OSError:
                return []

        jobs = _collect_jobs(self, peptides, collect, title)
        if not jobs:
            if jobs is not None:
                QMessageBox.information(self, title, "No PSMs found.")
            return
        answer = QMessageBox.question(self, title, f"Export {len(jobs)} PSMs of {len(peptides)} peptides, one PSM per page?")
        if answer != QMessageBox.StandardButton.Yes:
            return
        path = _ask_export_path(self, title, f"{self.specset_path.stem}_all.pdf", ".pdf", "PDF files (*.pdf)")
        if not path:
            return

        spectrum = _SpectrumWidget()  # off screen, used to draw the plots

        def write_page(report, job):
            pep, psm_path, file_name = job
            try:
                psm, frags = _read_psm(psm_path)
            except (OSError, ValueError) as e:
                report.text(f"{pep}: cannot read PSM {file_name}: {e}", 10, color="red")
                return
            report.text(f"{psm.get('pep', pep)}    charge {psm.get('z', '?')}", 13, bold=True)
            report.text(f"PSM {file_name}", 8, color="gray")
            spectrum.set_fragments(frags)
            report.plot(spectrum)
            report.fields_table("Fields", _field_rows(psm, frags))
            report.fragments_table("Fragments", _fragment_rows(frags))

        try:
            completed = _write_pages(self, path, jobs, write_page, title)
        except Exception as e:
            QMessageBox.critical(self, title, f"Export failed: {e}")
            return
        self.window().statusBar().showMessage(f"Exported {len(jobs)} PSMs to {path}" if completed else "Export cancelled")


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
        export_pdf = QPushButton("Export PDF…")
        export_pdf.setToolTip("Save the mirror plot and the tables of both PSMs as a PDF")
        export_pdf.clicked.connect(self._on_export_pdf)
        export_svg = QPushButton("Export SVG…")
        export_svg.setToolTip("Save the mirror plot as an SVG")
        export_svg.clicked.connect(self._on_export_svg)
        export_all = QPushButton("Export all PDF…")
        export_all.setToolTip("Save every pair of PSMs with a common charge, of the overlapping peptides currently listed, "
                              "as a PDF, one pair per page")
        export_all.clicked.connect(self._on_export_all_pdf)
        path_row = QHBoxLayout()
        path_row.addLayout(paths, 1)
        path_row.addWidget(compare, 0, Qt.AlignmentFlag.AlignBottom)
        path_row.addSpacing(16)
        path_row.addWidget(export_pdf, 0, Qt.AlignmentFlag.AlignBottom)
        path_row.addWidget(export_svg, 0, Qt.AlignmentFlag.AlignBottom)
        path_row.addWidget(export_all, 0, Qt.AlignmentFlag.AlignBottom)

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

        # read all of A; of B, keep only the peptides in A, so that a large B is not loaded into memory
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            index_a = _read_index(paths[0])
            locations_a = dict(index_a)
            locations_b = {}
            n_b = 0
            for pep, location in _iter_index(paths[1]):
                n_b += 1
                if pep in locations_a:
                    locations_b[pep] = location
        finally:
            QApplication.restoreOverrideCursor()

        self.specset_paths = paths
        self.locations = [locations_a, locations_b]
        overlap = [(pep,) for pep, _ in index_a if pep in locations_b]
        self._totals = f"A: {len(locations_a)}, B: {n_b}"

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

    def _selected_psm(self):
        """Fields of the first selected PSM (A, else B), or None."""
        return next((panel.psm for panel in self.psm_panels if panel.psm is not None), None)

    def _on_export_pdf(self):
        psm = self._selected_psm()
        if psm is None:
            QMessageBox.information(self, "Export PDF", "Select a PSM first.")
            return
        path = _ask_export_path(self, "Export PDF", _export_name(psm) + "_compare.pdf", ".pdf", "PDF files (*.pdf)")
        if not path:
            return
        try:
            report = _PdfReport(path)
            try:
                report.text(f"{psm.get('pep', '?')}    charge {psm.get('z', '?')}", 13, bold=True)
                report.plot(self.spectrum)
                for side, panel in zip(self._SIDES, self.psm_panels):
                    if panel.psm is not None:
                        panel.add_tables_to(report, f"Specset {side}: ")
            finally:
                report.close()
        except Exception as e:
            QMessageBox.critical(self, "Export PDF", f"Export failed: {e}")
            return
        self.window().statusBar().showMessage(f"Exported {path}")

    def _on_export_svg(self):
        psm = self._selected_psm()
        if psm is None:
            QMessageBox.information(self, "Export SVG", "Select a PSM first.")
            return
        path = _ask_export_path(self, "Export SVG", _export_name(psm) + "_compare.svg", ".svg", "SVG files (*.svg)")
        if not path:
            return
        try:
            _export_svg(self.spectrum, path, f"{psm.get('pep', '?')} charge {psm.get('z', '?')}: A up, B down")
        except Exception as e:
            QMessageBox.critical(self, "Export SVG", f"Export failed: {e}")
            return
        self.window().statusBar().showMessage(f"Exported {path}")

    def _on_export_all_pdf(self):
        title = "Export all PDF"
        peptides = [row[0] for row in _shown_rows(self.pep_proxy, self.pep_model)]
        if None in self.specset_paths or not peptides:
            QMessageBox.information(self, title, "Compare two specsets first; the overlapping peptides currently listed are exported.")
            return

        def collect(pep):
            """One job per pair of PSMs of A and B with the same charge."""
            location_paths = [self.specset_paths[i] / self.locations[i][pep] for i in range(2)]
            try:
                psms = [_read_location_index(path) for path in location_paths]
            except OSError:
                return []
            jobs = []
            for z in sorted({c for _, c in psms[0]} & {c for _, c in psms[1]}):
                files = [[f for f, c in psms[i] if c == z] for i in range(2)]
                for file_a in files[0]:
                    for file_b in files[1]:
                        jobs.append((pep, z, location_paths[0] / f"{file_a}.tsv", file_a, location_paths[1] / f"{file_b}.tsv", file_b))
            return jobs

        jobs = _collect_jobs(self, peptides, collect, title)
        if not jobs:
            if jobs is not None:
                QMessageBox.information(self, title, "No pairs of PSMs with a common charge found.")
            return
        answer = QMessageBox.question(self, title, f"Export {len(jobs)} pairs of PSMs of {len(peptides)} peptides, one pair per page?")
        if answer != QMessageBox.StandardButton.Yes:
            return
        name = f"{self.specset_paths[0].stem}_vs_{self.specset_paths[1].stem}.pdf"
        path = _ask_export_path(self, title, name, ".pdf", "PDF files (*.pdf)")
        if not path:
            return

        spectrum = _SpectrumWidget()  # off screen, used to draw the mirror plots

        def write_page(report, job):
            pep, z, path_a, file_a, path_b, file_b = job
            report.text(f"{pep}    charge {z}", 13, bold=True)
            report.text(f"A (up): PSM {file_a}    B (down): PSM {file_b}", 8, color="gray")
            frags = []
            for side, psm_path, file_name in (("A", path_a, file_a), ("B", path_b, file_b)):
                try:
                    frags.append(_read_psm(psm_path)[1])
                except (OSError, ValueError) as e:
                    report.text(f"Specset {side}: cannot read PSM {file_name}: {e}", 9, color="red")
                    frags.append([])
            spectrum.set_fragments(frags[0], frags[1])
            report.plot(spectrum)
            for side, side_frags in zip(self._SIDES, frags):
                report.fragments_table(f"Specset {side}: Fragments", _fragment_rows(side_frags))

        try:
            completed = _write_pages(self, path, jobs, write_page, title)
        except Exception as e:
            QMessageBox.critical(self, title, f"Export failed: {e}")
            return
        self.window().statusBar().showMessage(f"Exported {len(jobs)} pairs of PSMs to {path}" if completed else "Export cancelled")


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


SPECSET_VERSION = "0.9.3"

# Third-party software used by specset: (name, version, licence, URL)
ACKNOWLEDGEMENTS = [
    ("lxml", "6.1.3", "BSD-3-Clause", "https://lxml.de/"),    
    ("pandas", "3.0.6", "BSD-3-Clause", "https://pandas.pydata.org/"),
    ("psm_utils", "1.5.5", "Apache-2.0", "https://github.com/CompOmics/psm_utils"),
    ("pyarrow (Apache Arrow)", "26.0.0", "Apache-2.0", "https://arrow.apache.org/"),
    ("PySide6 (Qt for Python)", "6.12.0", "LGPL-3.0", "https://doc.qt.io/qtforpython-6/"),
    ("pyteomics", "4.7.5", "Apache-2.0", "https://pyteomics.readthedocs.io/"),
    ("Unimod database (unimod.xml)", "Downloaded 9 October 2026", "Design Science License", "https://www.unimod.org/"),
]


class _AboutWidget(QWidget):

    def __init__(self):
        super().__init__()
        rows = "".join(
            f"<tr><td><a href='{url}'>{html.escape(name)}</a></td><td>{html.escape(version)}</td>"
            f"<td>{html.escape(licence)}</td></tr>"
            for name, version, licence, url in ACKNOWLEDGEMENTS
        )
        label = QLabel(
            f"<h2>specset {html.escape(SPECSET_VERSION)}</h2>"
            "<p>Toolkit for manipulation of sets of mass spectra.</p>"
            "<h3>Acknowledgements</h3>"
            "<p>specset was built and tested on the following open-source software. "
            "We thank their authors and contributors.</p>"            
            "<table cellspacing='0' cellpadding='8'>"
            "<tr><th align='left'>Package</th><th align='left'>Version</th><th align='left'>Licence</th></tr>"
            f"{rows}</table>"            
        )
        label.setTextFormat(Qt.RichText)
        label.setOpenExternalLinks(True)
        label.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        label.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.addWidget(label)


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
        tabs.addTab(_AboutWidget(), "About")
        self.setCentralWidget(tabs)
        self.statusBar().showMessage("Ready")

    def closeEvent(self, event):
        if self.tools.is_running():
            QMessageBox.information(self, "Task running", "A task is still running. Please wait until it finishes.")
            event.ignore()
        else:
            event.accept()


# Application icon, used for the window and the taskbar if present
ICON_FILE = Path(__file__).with_name("favicon.ico")


def main():
    if sys.platform == "win32":
        # own taskbar entry with our icon, instead of being grouped under python.exe
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("specset.gui")
    app = QApplication(sys.argv)
    if ICON_FILE.is_file():
        app.setWindowIcon(QIcon(str(ICON_FILE)))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
