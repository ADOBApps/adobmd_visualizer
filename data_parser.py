"""
Parser for ADOBMD unified data files

This file is part of Quantum Analysis Helper.
Quantum Analysis Helper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Copyright (C) [2026] Acxel David Orozco Baldomero
"""

import re
import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path

@dataclass
class Atom:
    """Atom data structure"""
    id: int
    type_id: int
    molecule: int
    element: str
    x: float
    y: float
    z: float
    charge: float
    is_qm: bool = False
    color: Tuple[float, float, float] = (0.5, 0.5, 0.5)
    radius: float = 0.5

@dataclass
class Bond:
    """Bond data structure"""
    id: int
    type_id: int
    atom1: int
    atom2: int
    order: int = 1  # 1=single, 2=double, 3=triple

class ADOBMDData:
    """ADOBMD unified data file parser"""
    
    # Element colors (RGB) for visualization
    ELEMENT_COLORS = {
        'H': (1.0, 1.0, 1.0),    # White
        'C': (0.2, 0.2, 0.2),    # Dark gray
        'N': (0.0, 0.0, 1.0),    # Blue
        'O': (1.0, 0.0, 0.0),    # Red
        'F': (0.0, 1.0, 0.0),    # Green
        'CL': (0.0, 1.0, 0.0),   # Green
        'BR': (0.5, 0.2, 0.0),   # Brown
        'I': (0.6, 0.0, 0.8),    # Purple
        'S': (1.0, 1.0, 0.0),    # Yellow
        'P': (1.0, 0.5, 0.0),    # Orange
        'SI': (0.5, 0.5, 0.5),   # Gray
        'FE': (1.0, 0.5, 0.0),   # Orange-brown
        'CU': (0.8, 0.5, 0.2),   # Copper
        'ZN': (0.6, 0.6, 0.8),   # Light blue-gray
        'AU': (1.0, 0.8, 0.0),   # Gold
        'AG': (0.8, 0.8, 0.8),   # Silver
        'CD': (0.0, 0.5, 0.5),   # Teal
        'TE': (0.8, 0.5, 0.2),   # Brown
        'SE': (0.8, 0.4, 0.0),   # Orange-brown
    }
    
    # Element radii (Angstroms) for visualization
    ELEMENT_RADII = {
        'H': 0.3, 'C': 0.7, 'N': 0.65, 'O': 0.6, 'F': 0.5,
        'CL': 1.0, 'BR': 1.15, 'I': 1.4, 'S': 1.0, 'P': 1.0,
        'SI': 1.1, 'FE': 1.4, 'CU': 1.4, 'ZN': 1.4, 'AU': 1.5,
        'AG': 1.5, 'CD': 1.5, 'TE': 1.4, 'SE': 1.4,
    }
    
    # QM/MM region colors
    QM_COLOR = (1.0, 0.2, 0.2)    # Red for QM
    MM_COLOR = (0.2, 0.6, 1.0)    # Blue for MM
    
    def __init__(self):
        self.filename = ""
        self.title = ""
        self.comments: List[str] = []
        
        # Counts
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.nimpropers = 0
        self.natom_types = 0
        self.nbond_types = 0
        
        # Box
        self.box_lo = np.array([0.0, 0.0, 0.0])
        self.box_hi = np.array([0.0, 0.0, 0.0])
        self.has_box = False
        
        # Masses
        self.masses: Dict[int, float] = {}  # type_id -> mass
        self.type_elements: Dict[int, str] = {}  # type_id -> element
        
        # Atoms and bonds
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        
        # QM region
        self.qm_indices: List[int] = []
        self.has_qm_region = False
    
    def load(self, filename: str) -> bool:
        """Load ADOBMD data file"""
        try:
            self.filename = str(Path(filename).name)
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            # Parse header comments
            i = 0
            while i < len(lines) and lines[i].startswith('#'):
                if 'Title:' in lines[i]:
                    self.title = lines[i].split('Title:')[-1].strip()
                elif 'Comment:' in lines[i]:
                    self.comments.append(lines[i].split('Comment:')[-1].strip())
                else:
                    self.comments.append(lines[i][1:].strip())
                i += 1
            
            # Parse counts
            while i < len(lines):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                
                if 'atoms' in line:
                    self.natoms = int(line.split()[0])
                elif 'bonds' in line:
                    self.nbonds = int(line.split()[0])
                elif 'angles' in line:
                    self.nangles = int(line.split()[0])
                elif 'dihedrals' in line:
                    self.ndihedrals = int(line.split()[0])
                elif 'impropers' in line:
                    self.nimpropers = int(line.split()[0])
                elif 'atom types' in line:
                    self.natom_types = int(line.split()[0])
                elif 'bond types' in line:
                    self.nbond_types = int(line.split()[0])
                elif 'xlo xhi' in line:
                    vals = line.split()
                    self.box_lo[0] = float(vals[0])
                    self.box_hi[0] = float(vals[1])
                    self.has_box = True
                elif 'ylo yhi' in line:
                    vals = line.split()
                    self.box_lo[1] = float(vals[0])
                    self.box_hi[1] = float(vals[1])
                elif 'zlo zhi' in line:
                    vals = line.split()
                    self.box_lo[2] = float(vals[0])
                    self.box_hi[2] = float(vals[1])
                elif line == 'Masses':
                    i = self._parse_masses(lines, i + 1)
                elif line == 'Atoms':
                    i = self._parse_atoms(lines, i + 1)
                elif line == 'Bonds':
                    i = self._parse_bonds(lines, i + 1)
                else:
                    i += 1
            
            return True
            
        except Exception as e:
            print(f"Error loading file: {e}")
            return False
    
    def _parse_masses(self, lines: List[str], start: int) -> int:
        """Parse Masses section"""
        i = start
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('#'):
                i += 1
                continue
            
            parts = line.split('#')[0].strip().split()
            if len(parts) >= 2:
                type_id = int(parts[0])
                mass = float(parts[1])
                self.masses[type_id] = mass
                
                # Try to extract element from comment
                if '#' in line:
                    comment = line.split('#')[1].strip()
                    # Look for element symbol (1-2 letters)
                    import re
                    match = re.search(r'([A-Z][a-z]?)', comment)
                    if match:
                        self.type_elements[type_id] = match.group(1)
                    else:
                        self.type_elements[type_id] = 'X'
                else:
                    self.type_elements[type_id] = 'X'
            
            i += 1
            if i < len(lines) and lines[i].strip() and not lines[i].strip()[0].isdigit():
                break
        
        return i
    
    def _parse_atoms(self, lines: List[str], start: int) -> int:
        """Parse Atoms section"""
        i = start
        atom_count = 0
        
        while i < len(lines) and atom_count < self.natoms:
            line = lines[i].strip()
            if not line or line.startswith('#'):
                i += 1
                continue
            
            # Parse atom line
            parts = line.split()
            is_qm = False
            
            if len(parts) >= 7:
                # Format: id mol type q x y z [region]
                atom_id = int(parts[0])
                molecule = int(parts[1])
                type_id = int(parts[2])
                charge = float(parts[3])
                x = float(parts[4])
                y = float(parts[5])
                z = float(parts[6])
                
                # Check for region tag
                if len(parts) >= 8:
                    region = parts[7].upper()
                    is_qm = (region == 'QM' or region == 'Q')
                
                # Get element from type
                element = self.type_elements.get(type_id, 'X')
                
                # Set color based on region or element
                if is_qm:
                    color = self.QM_COLOR
                elif self.has_qm_region:
                    color = self.MM_COLOR
                else:
                    color = self.ELEMENT_COLORS.get(element, (0.7, 0.7, 0.7))
                
                # Set radius based on element
                radius = self.ELEMENT_RADII.get(element, 0.5)
                
                atom = Atom(
                    id=atom_id,
                    type_id=type_id,
                    molecule=molecule,
                    element=element,
                    x=x, y=y, z=z,
                    charge=charge,
                    is_qm=is_qm,
                    color=color,
                    radius=radius
                )
                
                self.atoms.append(atom)
                atom_count += 1
                
                if is_qm:
                    self.qm_indices.append(atom_id)
            
            i += 1
        
        self.has_qm_region = len(self.qm_indices) > 0
        return i
    
    def _parse_bonds(self, lines: List[str], start: int) -> int:
        """Parse Bonds section"""
        i = start
        bond_count = 0
        
        while i < len(lines) and bond_count < self.nbonds:
            line = lines[i].strip()
            if not line or line.startswith('#'):
                i += 1
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                bond_id = int(parts[0])
                bond_type = int(parts[1])
                atom1 = int(parts[2])
                atom2 = int(parts[3])
                
                bond = Bond(
                    id=bond_id,
                    type_id=bond_type,
                    atom1=atom1,
                    atom2=atom2,
                    order=bond_type  # For now, use type as order
                )
                
                self.bonds.append(bond)
                bond_count += 1
            
            i += 1
        
        return i
    
    def get_center(self) -> np.ndarray:
        """Get center of atomic coordinates"""
        if not self.atoms:
            return np.zeros(3)
        
        coords = np.array([[a.x, a.y, a.z] for a in self.atoms])
        return np.mean(coords, axis=0)
    
    def get_extent(self) -> float:
        """Get maximum extent from center"""
        if not self.atoms:
            return 10.0
        
        center = self.get_center()
        coords = np.array([[a.x, a.y, a.z] for a in self.atoms])
        dists = np.linalg.norm(coords - center, axis=1)
        return float(np.max(dists)) + 2.0
    
    def get_bond_coordinates(self, bond: Bond) -> Tuple[np.ndarray, np.ndarray]:
        """Get coordinates of bond endpoints"""
        atom1 = next((a for a in self.atoms if a.id == bond.atom1), None)
        atom2 = next((a for a in self.atoms if a.id == bond.atom2), None)
        
        if atom1 and atom2:
            return np.array([atom1.x, atom1.y, atom1.z]), np.array([atom2.x, atom2.y, atom2.z])
        return None, None
    
    def get_statistics(self) -> Dict:
        """Get system statistics"""
        elements = {}
        for atom in self.atoms:
            elements[atom.element] = elements.get(atom.element, 0) + 1
        
        return {
            'atoms': len(self.atoms),
            'bonds': len(self.bonds),
            'qm_atoms': len(self.qm_indices),
            'mm_atoms': len(self.atoms) - len(self.qm_indices),
            'elements': elements,
            'box': [self.box_hi[i] - self.box_lo[i] for i in range(3)] if self.has_box else None
        }