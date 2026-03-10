"""
NetworkX-based molecular visualizer with Fortran ULTRA-FAST parser
"""

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import networkx as nx
from matplotlib.collections import LineCollection

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QLabel, QFileDialog, QGroupBox, QCheckBox, QApplication,
    QComboBox, QSpinBox, QTextEdit, QSplitter, QMessageBox,
    QRadioButton, QToolButton, QMenu, QProgressBar
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont, QAction

import time
from pathlib import Path

from .data_parser import ADOBMDData, Atom, Bond
from .exporter import ADOBMDExporter

# Try to import Fortran parser (ULTRA FAST)
try:
    from ..fortran_parser import get_fortran_parser
    FORTRAN_AVAILABLE = True
except ImportError:
    FORTRAN_AVAILABLE = False
    print("[Warning] Fortran parser not available")

# Try to import LuaJIT parser (FAST)
try:
    from ..fast_parser import LuaJITFastParser
    LUAJIT_AVAILABLE = True
except ImportError:
    LUAJIT_AVAILABLE = False
    print("[Warning] LuaJIT parser not available")

class GraphMolecularCanvas(FigureCanvas):
    """Molecular visualizer using NetworkX graph"""
    
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        
        # Graph structure
        self.graph = None
        self.pos = None
        self.data = None
        self.selected_atom = None
        self.highlighted_nodes = []
        
        # Current view
        self.current_view = 0  # 0=xy, 1=xz, 2=yz
        self.views = ['xy', 'xz', 'yz']
        
        # Style
        self.fig.patch.set_facecolor('white')
        self.ax.set_facecolor('white')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_aspect('equal')
        
        # Connect pick event
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
    
    def build_graph(self, data: ADOBMDData):
        """Build NetworkX graph from molecular data"""
        self.data = data
        self.graph = nx.Graph()
        
        # Add nodes (atoms) with attributes
        for atom in data.atoms:
            self.graph.add_node(
                atom.id,
                element=atom.element,
                charge=atom.charge,
                is_qm=atom.is_qm,
                x=atom.x, y=atom.y, z=atom.z,
                radius=data.ELEMENT_RADII.get(atom.element, 0.5),
                color=self._get_atom_color(atom, data)
            )
        
        # Add edges (bonds) with attributes
        for bond in data.bonds:
            if bond.atom1 in self.graph and bond.atom2 in self.graph:
                self.graph.add_edge(
                    bond.atom1, bond.atom2,
                    order=bond.order,
                    type=bond.type_id
                )
        
        # Create position dictionary for current view
        self.update_positions()
        
        return self.graph
    
    def _get_atom_color(self, atom, data):
        """Get color for atom based on attributes"""
        if atom.is_qm:
            return 'red'
        elif data.has_qm_region:
            return 'blue'
        else:
            return data.ELEMENT_COLORS.get(atom.element, 'gray')
    
    def update_positions(self):
        """Update 2D positions based on current view"""
        if not self.graph:
            return
        
        self.pos = {}
        view = self.views[self.current_view]
        
        for node, attr in self.graph.nodes(data=True):
            if view == 'xy':
                self.pos[node] = (attr['x'], attr['y'])
            elif view == 'xz':
                self.pos[node] = (attr['x'], attr['z'])
            else:  # yz
                self.pos[node] = (attr['y'], attr['z'])
    
    def plot_system(self, data: ADOBMDData, show_bonds=True, show_box=True,
                   atom_size=300, color_by='region', highlight_atoms=None):
        """
        Plot molecular system using NetworkX (ULTRA FAST)
        
        This uses nx.draw_networkx which renders everything in ONE draw call
        """
        self.ax.clear()
        self.ax.grid(True, alpha=0.3)
        
        # Build graph
        self.build_graph(data)
        
        # Prepare node colors
        node_colors = []
        highlight_colors = []
        
        for node, attr in self.graph.nodes(data=True):
            # Base color
            if color_by == 'region':
                color = attr['color']
            elif color_by == 'element':
                color = data.ELEMENT_COLORS.get(attr['element'], 'gray')
            elif color_by == 'charge':
                if attr['charge'] < -0.1:
                    color = 'darkblue'
                elif attr['charge'] > 0.1:
                    color = 'darkred'
                else:
                    color = 'gray'
            else:
                color = 'gray'
            
            # Highlight if selected
            if highlight_atoms and node in highlight_atoms:
                highlight_colors.append('yellow')
                node_colors.append(color)  # Keep original for reference
            else:
                node_colors.append(color)
                highlight_colors.append(color)
        
        # Calculate node sizes based on atomic radius
        node_sizes = [attr['radius'] * atom_size for attr in self.graph.nodes(data=True)]
        
        # Draw regular nodes
        nx.draw_networkx_nodes(
            self.graph, self.pos,
            node_color=node_colors,
            node_size=node_sizes,
            alpha=0.7,
            linewidths=0.5,
            edgecolors='black',
            ax=self.ax
        )
        
        # Draw highlighted nodes on top (if any)
        if highlight_atoms:
            highlight_nodes = [n for n in self.graph.nodes() if n in highlight_atoms]
            if highlight_nodes:
                highlight_sizes = [node_sizes[i] * 1.2 for i, n in enumerate(self.graph.nodes()) 
                                 if n in highlight_atoms]
                nx.draw_networkx_nodes(
                    self.graph, self.pos,
                    nodelist=highlight_nodes,
                    node_color='yellow',
                    node_size=highlight_sizes,
                    alpha=1.0,
                    linewidths=2,
                    edgecolors='red',
                    ax=self.ax
                )
        
        # Draw edges if requested (SINGLE DRAW CALL!)
        if show_bonds and data.bonds and self.graph.edges():
            # Color edges by bond order
            edge_colors = []
            for u, v, attr in self.graph.edges(data=True):
                if attr['order'] == 1:
                    edge_colors.append('gray')
                elif attr['order'] == 2:
                    edge_colors.append('orange')
                elif attr['order'] == 3:
                    edge_colors.append('red')
                else:
                    edge_colors.append('purple')
            
            nx.draw_networkx_edges(
                self.graph, self.pos,
                edge_color=edge_colors,
                width=1.0,
                alpha=0.5,
                ax=self.ax
            )
        
        # Draw box if requested
        if show_box and data.has_box:
            self._draw_box(data)
        
        # Labels
        xlabel, ylabel = self._get_axis_labels()
        self.ax.set_xlabel(f'{xlabel} (Å)')
        self.ax.set_ylabel(f'{ylabel} (Å)')
        
        # Title
        qm_info = f" - QM: {len(data.qm_indices)} atoms" if data.has_qm_region else ""
        parser_info = " [Fortran]" if hasattr(self, 'using_fortran') and self.using_fortran else ""
        self.ax.set_title(f"{data.filename}{parser_info} - {len(data.atoms)} atoms, "
                         f"{len(data.bonds)} bonds{qm_info}")
        
        # Auto-scale
        self._auto_scale()
        
        self.fig.tight_layout()
        self.draw()
    
    def _draw_box(self, data):
        """Draw 2D projection of simulation box"""
        lo = data.box_lo
        hi = data.box_hi
        
        # Project corners to 2D
        corners_3d = [
            [lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]],
            [hi[0], hi[1], lo[2]], [lo[0], hi[1], lo[2]],
            [lo[0], lo[1], hi[2]], [hi[0], lo[1], hi[2]],
            [hi[0], hi[1], hi[2]], [lo[0], hi[1], hi[2]]
        ]
        
        view = self.views[self.current_view]
        corners = []
        for c in corners_3d:
            if view == 'xy':
                corners.append([c[0], c[1]])
            elif view == 'xz':
                corners.append([c[0], c[2]])
            else:
                corners.append([c[1], c[2]])
        
        # Define edges
        edges = [(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4),
                 (0,4), (1,5), (2,6), (3,7)]
        
        # Draw edges as line collection (ONE DRAW CALL)
        segments = [[corners[i], corners[j]] for i, j in edges]
        lc = LineCollection(segments, colors='cyan', linewidths=0.5, alpha=0.3)
        self.ax.add_collection(lc)
    
    def _auto_scale(self):
        """Auto-scale with padding"""
        if not self.pos:
            return
        
        x_coords = [p[0] for p in self.pos.values()]
        y_coords = [p[1] for p in self.pos.values()]
        
        if not x_coords or not y_coords:
            return
        
        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)
        
        # Handle single atom case
        if x_min == x_max:
            x_min -= 5
            x_max += 5
        if y_min == y_max:
            y_min -= 5
            y_max += 5
        
        x_pad = (x_max - x_min) * 0.1
        y_pad = (y_max - y_min) * 0.1
        
        self.ax.set_xlim(x_min - x_pad, x_max + x_pad)
        self.ax.set_ylim(y_min - y_pad, y_max + y_pad)
    
    def _get_axis_labels(self):
        """Get axis labels for current view"""
        view = self.views[self.current_view]
        if view == 'xy':
            return 'X', 'Y'
        elif view == 'xz':
            return 'X', 'Z'
        else:
            return 'Y', 'Z'
    
    def cycle_view(self):
        """Switch to next view"""
        self.current_view = (self.current_view + 1) % 3
        if self.graph:
            self.update_positions()
            self.plot_system(self.data)
        return self.views[self.current_view]
    
    def on_pick(self, event):
        """Handle atom selection"""
        if not hasattr(event, 'ind') or not self.graph:
            return
        
        # Get selected node
        node = list(self.graph.nodes())[event.ind[0]]
        self.selected_atom = node
        
        # Highlight the selected atom
        self.plot_system(self.data, highlight_atoms=[node])
    
    def find_shortest_path(self, atom1, atom2):
        """Find shortest path between atoms (graph algorithm)"""
        if self.graph and nx.has_path(self.graph, atom1, atom2):
            return nx.shortest_path(self.graph, atom1, atom2)
        return None
    
    def get_degree_distribution(self):
        """Get degree distribution of the molecular graph"""
        if self.graph:
            return dict(self.graph.degree())
        return {}
    
    def find_cycles(self):
        """Find cycles in the molecular graph"""
        if self.graph:
            return list(nx.cycle_basis(self.graph))
        return []

class VisualizerWidget(QWidget):
    """Main widget with NetworkX-powered visualization and Fortran parser"""
    
    def __init__(self, plugin_instance=None):
        super().__init__()
        self.plugin = plugin_instance
        self.data = None
        self.exporter = None
        self.current_file = ""
        
        # Initialize parsers (Fortran FIRST - it's fastest!)
        self.fortran_parser = None
        self.luajit_parser = None
        
        # Try Fortran (ULTRA FAST)
        if FORTRAN_AVAILABLE:
            try:
                self.fortran_parser = get_fortran_parser()
                if self.fortran_parser:
                    print("[Ok] Fortran ultra-fast parser available!")
                else:
                    print("[Warning] Fortran parser not available")
            except Exception as e:
                print(f"[Warning] Fortran parser error: {e}")
        
        # Try LuaJIT (FAST) as fallback
        if not self.fortran_parser and LUAJIT_AVAILABLE:
            try:
                from ..fast_parser import LuaJITFastParser
                self.luajit_parser = LuaJITFastParser()
                if self.luajit_parser.initialize():
                    print("[Ok] LuaJIT fast parser available!")
                else:
                    self.luajit_parser = None
            except Exception as e:
                print(f"[Warning] LuaJIT parser error: {e}")
        
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Status bar for parser info
        status_bar = QHBoxLayout()
        parser_status = "Fortran" if self.fortran_parser else "⚡ LuaJIT" if self.luajit_parser else "Python"
        status_label = QLabel(f"Parser: {parser_status}")
        status_label.setStyleSheet("background-color: #e0e0e0; padding: 2px 10px; border-radius: 3px;")
        status_bar.addWidget(status_label)
        status_bar.addStretch()
        layout.addLayout(status_bar)
        
        # Toolbar
        toolbar = QHBoxLayout()
        
        self.load_btn = QPushButton("Load File")
        self.load_btn.clicked.connect(self.load_file)
        self.load_btn.setFixedWidth(120)
        toolbar.addWidget(self.load_btn)
        
        self.view_btn = QPushButton("View: XY")
        self.view_btn.clicked.connect(self.cycle_view)
        self.view_btn.setFixedWidth(100)
        toolbar.addWidget(self.view_btn)
        
        self.analyze_btn = QPushButton("Analyze Graph")
        self.analyze_btn.clicked.connect(self.show_graph_analysis)
        self.analyze_btn.setEnabled(False)
        self.analyze_btn.setFixedWidth(120)
        toolbar.addWidget(self.analyze_btn)
        
        toolbar.addStretch()
        
        self.file_label = QLabel("No file loaded")
        self.file_label.setStyleSheet("color: gray; padding: 5px;")
        toolbar.addWidget(self.file_label)
        
        layout.addLayout(toolbar)
        
        # Export toolbar
        export_bar = QHBoxLayout()
        export_bar.addWidget(QLabel("Export:"))
        
        formats = ['XYZ', 'SDF', 'PDB', 'MOL2', 'CIF', 'GRO']
        self.export_btns = {}
        for fmt in formats:
            btn = QPushButton(fmt)
            btn.clicked.connect(lambda checked, f=fmt.lower(): self.export_file(f))
            btn.setEnabled(False)
            btn.setFixedWidth(50)
            export_bar.addWidget(btn)
            self.export_btns[fmt.lower()] = btn
        
        export_bar.addStretch()
        
        self.viewer_menu = QToolButton()
        self.viewer_menu.setText("Launch Viewer ▼")
        self.viewer_menu.setPopupMode(QToolButton.InstantPopup)
        self.viewer_menu.setEnabled(False)
        
        menu = QMenu()
        for viewer in ['VMD', 'PyMOL', 'Avogadro', 'Jmol']:
            action = QAction(viewer, self)
            action.triggered.connect(lambda checked, v=viewer.lower(): self.launch_viewer(v))
            menu.addAction(action)
        self.viewer_menu.setMenu(menu)
        export_bar.addWidget(self.viewer_menu)
        
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        self.progress.setFixedWidth(150)
        export_bar.addWidget(self.progress)
        
        layout.addLayout(export_bar)
        
        # Main splitter
        splitter = QSplitter(Qt.Horizontal)
        
        # Left panel - Controls
        left = QWidget()
        left_layout = QVBoxLayout(left)
        
        # Display controls
        display_group = QGroupBox("Display")
        display = QVBoxLayout()
        
        # Color by
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color:"))
        self.color_combo = QComboBox()
        self.color_combo.addItems(["region", "element", "charge"])
        self.color_combo.currentTextChanged.connect(self.update_plot)
        color_row.addWidget(self.color_combo)
        display.addLayout(color_row)
        
        # Size
        size_row = QHBoxLayout()
        size_row.addWidget(QLabel("Size:"))
        self.size_spin = QSpinBox()
        self.size_spin.setRange(50, 1000)
        self.size_spin.setValue(300)
        self.size_spin.valueChanged.connect(self.update_plot)
        size_row.addWidget(self.size_spin)
        display.addLayout(size_row)
        
        # Options
        self.show_bonds = QCheckBox("Show Bonds")
        self.show_bonds.setChecked(True)
        self.show_bonds.toggled.connect(self.update_plot)
        display.addWidget(self.show_bonds)
        
        self.show_box = QCheckBox("Show Box")
        self.show_box.setChecked(True)
        self.show_box.toggled.connect(self.update_plot)
        display.addWidget(self.show_box)
        
        display_group.setLayout(display)
        left_layout.addWidget(display_group)
        
        # Filter controls
        filter_group = QGroupBox("Filter")
        filter = QVBoxLayout()
        
        self.show_all = QRadioButton("Show All")
        self.show_all.setChecked(True)
        self.show_all.toggled.connect(self.filter_atoms)
        filter.addWidget(self.show_all)
        
        self.show_qm = QRadioButton("QM Only")
        self.show_qm.toggled.connect(self.filter_atoms)
        filter.addWidget(self.show_qm)
        
        self.show_mm = QRadioButton("MM Only")
        self.show_mm.toggled.connect(self.filter_atoms)
        filter.addWidget(self.show_mm)
        
        filter_group.setLayout(filter)
        left_layout.addWidget(filter_group)
        
        # Statistics
        stats_group = QGroupBox("Statistics")
        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setMaximumHeight(200)
        self.stats_text.setFont(QFont("Courier New", 9))
        stats_group.setLayout(QVBoxLayout())
        stats_group.layout().addWidget(self.stats_text)
        left_layout.addWidget(stats_group)
        
        left_layout.addStretch()
        
        # Right panel - Canvas
        right = QWidget()
        right_layout = QVBoxLayout(right)
        self.canvas = GraphMolecularCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        right_layout.addWidget(self.toolbar)
        right_layout.addWidget(self.canvas)
        
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([350, 650])
        
        layout.addWidget(splitter)
    
    def load_file_fortran(self, file_path):
        """Load file using Fortran (ULTRA FAST - microseconds!)"""
        try:
            result = self.fortran_parser.parse_file(file_path)
            if not result:
                return False
            
            # Convert to ADOBMDData
            self.data = ADOBMDData()
            self.data.filename = Path(file_path).name
            
            # Convert atoms
            self.data.atoms = []
            for a in result['atoms']:
                atom = Atom(
                    id=int(a['id']),
                    type_id=int(a['type_id']),
                    molecule=int(a['molecule']),
                    element=str(a['element']),
                    x=float(a['x']),
                    y=float(a['y']),
                    z=float(a['z']),
                    charge=float(a['charge']),
                    is_qm=bool(a['is_qm'])
                )
                self.data.atoms.append(atom)
            
            # Convert bonds
            self.data.bonds = []
            if 'bonds' in result:
                for b in result['bonds']:
                    bond = Bond(
                        id=int(b['id']),
                        type_id=int(b['type_id']),
                        atom1=int(b['atom1']),
                        atom2=int(b['atom2']),
                        order=int(b.get('order', b['type_id']))
                    )
                    self.data.bonds.append(bond)
            
            # Set box
            if result['header']['has_box']:
                self.data.has_box = True
                self.data.box_lo = result['header']['box_lo']
                self.data.box_hi = result['header']['box_hi']
            
            # Set QM indices
            self.data.qm_indices = result.get('qm_indices', [])
            self.data.has_qm_region = result.get('has_qm', False)
            
            # Set counts
            self.data.natoms = len(self.data.atoms)
            self.data.nbonds = len(self.data.bonds)
            
            # Mark that we used Fortran
            self.canvas.using_fortran = True
            return True
            
        except Exception as e:
            print(f"[Error] Fortran parser error: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def load_file_luajit(self, file_path):
        """Load file using LuaJIT (FAST)"""
        try:
            if not self.luajit_parser.parser_loaded:
                if not self.luajit_parser.initialize():
                    return False
            
            result = self.luajit_parser.parse_file(file_path)
            if not result:
                return False
            
            # Convert to ADOBMDData
            self.data = ADOBMDData()
            self.data.filename = Path(file_path).name
            
            # Convert atoms
            self.data.atoms = []
            for a in result['atoms']:
                atom = Atom(
                    id=int(a['id']),
                    type_id=int(a['type']),
                    molecule=int(a['mol']),
                    element=str(a['element']),
                    x=float(a['x']),
                    y=float(a['y']),
                    z=float(a['z']),
                    charge=float(a['charge']),
                    is_qm=bool(a['qm'])
                )
                self.data.atoms.append(atom)
            
            # Convert bonds
            self.data.bonds = []
            if 'bonds' in result:
                for b in result['bonds']:
                    bond = Bond(
                        id=int(b['id']),
                        type_id=int(b['type']),
                        atom1=int(b['a1']),
                        atom2=int(b['a2']),
                        order=int(b.get('order', b['type']))
                    )
                    self.data.bonds.append(bond)
            
            # Set box
            if result['header']['has_box']:
                self.data.has_box = True
                self.data.box_lo = result['header']['box_lo']
                self.data.box_hi = result['header']['box_hi']
            
            # Set QM indices
            self.data.qm_indices = result.get('qm_indices', [])
            self.data.has_qm_region = result.get('has_qm', False)
            
            self.data.natoms = len(self.data.atoms)
            self.data.nbonds = len(self.data.bonds)
            
            self.canvas.using_fortran = False
            return True
            
        except Exception as e:
            print(f"[Error] LuaJIT parser error: {e}")
            return False
    
    def load_file_python(self, file_path):
        """Load file using Python (slow fallback)"""
        self.data = ADOBMDData()
        if not self.data.load(file_path):
            return False
        self.canvas.using_fortran = False
        return True
    
    def load_file(self):
        """Load ADOBMD data file with automatic parser selection"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select ADOBMD Data File", "",
            "ADOBMD Files (*.data);;All Files (*)"
        )
        
        if not file_path:
            return
        
        self.current_file = file_path
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        QApplication.processEvents()
        
        start = time.time()
        success = False
        
        # Try parsers in order: Fortran -> LuaJIT -> Python
        if self.fortran_parser:
            success = self.load_file_fortran(file_path)
            parser_used = "Fortran"
        elif self.luajit_parser:
            success = self.load_file_luajit(file_path)
            parser_used = "LuaJIT"
        else:
            success = self.load_file_python(file_path)
            parser_used = "Python"
        
        if not success:
            QMessageBox.critical(self, "Error", "Failed to load file")
            self.progress.setVisible(False)
            return
        
        load_time = time.time() - start
        
        self.exporter = ADOBMDExporter(self.data)
        
        # Update file label with timing
        time_str = f"{load_time*1000:.1f}ms" if load_time < 1.0 else f"{load_time:.2f}s"
        self.file_label.setText(f"{self.data.filename} ({parser_used}: {time_str})")
        
        # Enable controls
        for btn in self.export_btns.values():
            btn.setEnabled(True)
        self.viewer_menu.setEnabled(True)
        self.analyze_btn.setEnabled(True)
        
        self.update_plot()
        self.update_statistics()
        
        self.progress.setVisible(False)
        
        # Show success message with timing
        speed = "⚡⚡ ULTRA FAST" if parser_used == "Fortran" else "⚡ FAST" if parser_used == "LuaJIT" else "🐍"
        QMessageBox.information(self, "Success", 
                               f"{speed} Loaded {len(self.data.atoms)} atoms, "
                               f"{len(self.data.bonds)} bonds in {time_str}")
    
    def update_plot(self):
        """Update visualization"""
        if not self.data:
            return
        
        # Filter atoms
        atoms_to_show = []
        if self.show_all.isChecked():
            atoms_to_show = self.data.atoms
        elif self.show_qm.isChecked():
            atoms_to_show = [a for a in self.data.atoms if a.is_qm]
        elif self.show_mm.isChecked():
            atoms_to_show = [a for a in self.data.atoms if not a.is_qm]
        
        # Create filtered data
        if len(atoms_to_show) != len(self.data.atoms):
            filtered = ADOBMDData()
            filtered.atoms = atoms_to_show
            filtered.bonds = self.data.bonds
            filtered.has_box = self.data.has_box
            filtered.box_lo = self.data.box_lo
            filtered.box_hi = self.data.box_hi
            filtered.filename = self.data.filename
            filtered.qm_indices = self.data.qm_indices
            filtered.has_qm_region = self.data.has_qm_region
            plot_data = filtered
        else:
            plot_data = self.data
        
        self.canvas.plot_system(
            plot_data,
            show_bonds=self.show_bonds.isChecked(),
            show_box=self.show_box.isChecked(),
            atom_size=self.size_spin.value(),
            color_by=self.color_combo.currentText()
        )
    
    def cycle_view(self):
        """Cycle through views"""
        view = self.canvas.cycle_view()
        self.view_btn.setText(f"View: {view.upper()}")
    
    def filter_atoms(self):
        """Filter atoms"""
        self.update_plot()
        self.update_statistics()
    
    def update_statistics(self):
        """Update statistics display"""
        if not self.data:
            return
        
        stats = self.data.get_statistics()
        
        # Graph analysis if available
        graph_stats = ""
        if self.canvas.graph:
            degrees = self.canvas.get_degree_distribution()
            avg_degree = sum(degrees.values()) / len(degrees) if degrees else 0
            cycles = self.canvas.find_cycles()
            
            graph_stats = f"\nGraph Analysis:\n"
            graph_stats += f"  Average degree: {avg_degree:.2f}\n"
            graph_stats += f"  Cycles found: {len(cycles)}\n"
            if cycles:
                graph_stats += f"  Largest cycle: {len(max(cycles, key=len))} atoms\n"
        
        # Parser info
        parser_info = ""
        if self.fortran_parser:
            parser_info = "\n Fortran acceleration"
        elif self.luajit_parser:
            parser_info = "\n⚡ LuaJIT acceleration"
        
        text = f"""File: {self.data.filename}
Atoms: {stats['atoms']}
Bonds: {stats['bonds']}
QM atoms: {stats['qm_atoms']}
MM atoms: {stats['mm_atoms']}

Elements:
{chr(10).join(f'  {k}: {v}' for k, v in sorted(stats['elements'].items()))}
{graph_stats}"""
        
        if stats['box']:
            box = stats['box']
            text += f"\nBox: {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f} Å³"
        
        text += parser_info
        
        self.stats_text.setText(text)
    
    def show_graph_analysis(self):
        """Show detailed graph analysis"""
        if not self.canvas.graph:
            return
        
        degrees = self.canvas.get_degree_distribution()
        cycles = self.canvas.find_cycles()
        
        # Calculate additional metrics
        avg_path_length = 0
        if nx.is_connected(self.canvas.graph):
            avg_path_length = nx.average_shortest_path_length(self.canvas.graph)
        
        msg = f"""📊 Molecular Graph Analysis

Nodes (atoms): {self.canvas.graph.number_of_nodes()}
Edges (bonds): {self.canvas.graph.number_of_edges()}
Connected: {nx.is_connected(self.canvas.graph)}
Average path length: {avg_path_length:.2f}

Degree Distribution:
{chr(10).join(f'  Atom {k}: {v} bonds' for k, v in sorted(degrees.items())[:20])}

Cycles found: {len(cycles)}
"""
        if cycles:
            msg += "\nLargest cycles:\n"
            sorted_cycles = sorted(cycles, key=len, reverse=True)
            for i, cycle in enumerate(sorted_cycles[:5]):
                msg += f"  Cycle {i+1}: {len(cycle)} atoms → {'-'.join(map(str, cycle[:5]))}...\n"
        
        QMessageBox.information(self, "Graph Analysis", msg)
    
    def export_file(self, format_type):
        """Export to specified format"""
        if not self.data or not self.exporter:
            return
        
        base = Path(self.current_file).stem
        ext_map = {
            'xyz': '.xyz', 'sdf': '.sdf', 'pdb': '.pdb',
            'mol2': '.mol2', 'cif': '.cif', 'gro': '.gro'
        }
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, f"Export as {format_type.upper()}",
            f"{base}{ext_map[format_type]}",
            f"{format_type.upper()} Files (*.{format_type})"
        )
        
        if file_path:
            try:
                getattr(self.exporter, f'to_{format_type}')(file_path)
                QMessageBox.information(self, "Success", f"Exported to: {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Export failed: {str(e)}")
    
    def launch_viewer(self, viewer):
        """Launch external viewer"""
        if self.exporter:
            self.exporter.launch_viewer(viewer)