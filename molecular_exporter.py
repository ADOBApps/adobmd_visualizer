"""
molecular_exporter.py - Python wrapper for Fortran molecular export library
Supports XYZ, SDF, PDB, MOL2, CIF, GRO + Schlegel diagram projections
"""

import ctypes
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import platform
import time

# Platform-specific library extension
if platform.system() == 'Windows':
    LIB_EXT = '.dll'
elif platform.system() == 'Darwin':
    LIB_EXT = '.dylib'
else:  # Linux
    LIB_EXT = '.so'

# Load Fortran library
lib_path = Path(__file__).parent / f"libmolexport{LIB_EXT}"
if lib_path.exists():
    lib = ctypes.CDLL(str(lib_path))
    print(f"[Ok] Loaded molecular export library: {lib_path}")
else:
    lib = None
    print(f"[Warning] Library not found: {lib_path}")

class ExportAtom(ctypes.Structure):
    """C-compatible atom structure for export"""
    _fields_ = [
        ('id', ctypes.c_int),
        ('type_id', ctypes.c_int),
        ('molecule', ctypes.c_int),
        ('element', ctypes.c_char * 3),
        ('x', ctypes.c_double),
        ('y', ctypes.c_double),
        ('z', ctypes.c_double),
        ('charge', ctypes.c_double),
        ('is_qm', ctypes.c_int),
        # Schlegel projections
        ('xy_x', ctypes.c_double),
        ('xy_y', ctypes.c_double),
        ('xz_x', ctypes.c_double),
        ('xz_z', ctypes.c_double),
        ('yz_y', ctypes.c_double),
        ('yz_z', ctypes.c_double),
        ('schlegel_x', ctypes.c_double),
        ('schlegel_y', ctypes.c_double),
    ]
    
    def to_dict(self):
        """Convert to dictionary"""
        return {
            'id': self.id,
            'type_id': self.type_id,
            'molecule': self.molecule,
            'element': self.element.decode('utf-8').strip('\x00'),
            'x': self.x,
            'y': self.y,
            'z': self.z,
            'charge': self.charge,
            'is_qm': bool(self.is_qm),
            'xy': (self.xy_x, self.xy_y),
            'xz': (self.xz_x, self.xz_z),
            'yz': (self.yz_y, self.yz_z),
            'schlegel': (self.schlegel_x, self.schlegel_y),
        }

class ExportBond(ctypes.Structure):
    """C-compatible bond structure"""
    _fields_ = [
        ('id', ctypes.c_int),
        ('type_id', ctypes.c_int),
        ('atom1', ctypes.c_int),
        ('atom2', ctypes.c_int),
        ('order', ctypes.c_int),
    ]
    
    def to_dict(self):
        return {
            'id': self.id,
            'type_id': self.type_id,
            'atom1': self.atom1,
            'atom2': self.atom2,
            'order': self.order,
        }

class MolecularExporter:
    """
    Molecular structure exporter with Schlegel diagram projections
    
    Exports to: XYZ, PDB, SDF, MOL2, CIF, GRO
    Schlegel projections: XY, XZ, YZ planes + perspective view
    """
    
    def __init__(self):
        self.lib = lib
        self.version = "Unknown"
        
        if self.lib:
            # Configure export_all_formats function
            self.lib.export_all_formats.argtypes = [
                ctypes.c_char_p,              # filename
                ctypes.c_int,                  # natoms
                ctypes.c_int,                  # nbonds
                ctypes.POINTER(ExportAtom),    # atoms
                ctypes.POINTER(ExportBond),    # bonds
                ctypes.c_int,                  # has_box
                ctypes.c_double * 3,           # box_lo
                ctypes.c_double * 3,           # box_hi
                ctypes.c_double * 3,           # projection_point
                ctypes.c_double,               # view_angle
            ]
            self.lib.export_all_formats.restype = None
            
            # Configure export_schlegel_only function
            self.lib.export_schlegel_only.argtypes = [
                ctypes.c_char_p,               # filename
                ctypes.c_int,                   # natoms
                ctypes.POINTER(ExportAtom),     # atoms
                ctypes.c_double * 3,            # projection_point
                ctypes.c_double,                # view_angle
            ]
            self.lib.export_schlegel_only.restype = None
            
            # Get version
            self.lib.get_version.argtypes = [ctypes.c_char_p]
            self.lib.get_version.restype = None
            version_buf = ctypes.create_string_buffer(256)
            self.lib.get_version(version_buf)
            self.version = version_buf.value.decode('utf-8')
            print(f"[Ok] Molecular exporter version: {self.version}")
    
    def export_all(self, base_filename: str, atoms: List[Dict], bonds: List[Dict],
                   box: Optional[Tuple[float, float, float]] = None,
                   projection_point: Optional[Tuple[float, float, float]] = None,
                   view_angle: float = 45.0) -> bool:
        """
        Export all formats including Schlegel projections
        
        Args:
            base_filename: Base name for output files (without extension)
            atoms: List of atom dictionaries
            bonds: List of bond dictionaries
            box: Optional (x, y, z) box dimensions
            projection_point: Point from which to project (default: above structure)
            view_angle: View angle in degrees for cone projection
        """
        if not self.lib:
            print("[Error] Library not loaded")
            return False
        
        start_time = time.time()
        
        # Convert to C-compatible arrays
        n_atoms = len(atoms)
        n_bonds = len(bonds)
        
        atom_array = (ExportAtom * n_atoms)()
        bond_array = (ExportBond * n_bonds)()
        
        for i, a in enumerate(atoms):
            atom_array[i].id = a['id']
            atom_array[i].type_id = a.get('type_id', 1)
            atom_array[i].molecule = a.get('molecule', 1)
            atom_array[i].element = a['element'].ljust(2)[:2].encode('utf-8') + b'\x00'
            atom_array[i].x = a['x']
            atom_array[i].y = a['y']
            atom_array[i].z = a['z']
            atom_array[i].charge = a.get('charge', 0.0)
            atom_array[i].is_qm = 1 if a.get('is_qm', False) else 0
        
        for i, b in enumerate(bonds):
            bond_array[i].id = b['id']
            bond_array[i].type_id = b.get('type_id', 1)
            bond_array[i].atom1 = b['atom1']
            bond_array[i].atom2 = b['atom2']
            bond_array[i].order = b.get('order', 1)
        
        # Box
        has_box = 1 if box else 0
        if box:
            box_lo = (ctypes.c_double * 3)(0.0, 0.0, 0.0)
            box_hi = (ctypes.c_double * 3)(box[0], box[1], box[2])
        else:
            # Calculate box from atoms
            xs = [a['x'] for a in atoms]
            ys = [a['y'] for a in atoms]
            zs = [a['z'] for a in atoms]
            min_coords = [min(xs), min(ys), min(zs)]
            max_coords = [max(xs), max(ys), max(zs)]
            padding = 5.0
            box_lo = (ctypes.c_double * 3)(min_coords[0] - padding,
                                            min_coords[1] - padding,
                                            min_coords[2] - padding)
            box_hi = (ctypes.c_double * 3)(max_coords[0] + padding,
                                            max_coords[1] + padding,
                                            max_coords[2] + padding)
        
        # Projection point (default: above the structure)
        if projection_point:
            proj = (ctypes.c_double * 3)(projection_point[0],
                                          projection_point[1],
                                          projection_point[2])
        else:
            # Place projection point above the highest atom
            max_z = max(a['z'] for a in atoms)
            proj = (ctypes.c_double * 3)(0.0, 0.0, max_z + 10.0)
        
        # Call Fortran
        filename_bytes = str(base_filename).encode('utf-8')
        self.lib.export_all_formats(
            filename_bytes,
            n_atoms,
            n_bonds,
            atom_array,
            bond_array,
            has_box,
            box_lo,
            box_hi,
            proj,
            view_angle
        )
        
        elapsed = time.time() - start_time
        print(f"✅ Exported all formats in {elapsed*1000:.1f}ms")
        print(f"   Files created: {base_filename}.* + projections")
        
        return True
    
    def export_schlegel_only(self, base_filename: str, atoms: List[Dict],
                             projection_point: Optional[Tuple[float, float, float]] = None,
                             view_angle: float = 45.0) -> bool:
        """
        Export only Schlegel diagram projections
        
        Creates: {base_filename}_xy.txt, _xz.txt, _yz.txt, _schlegel.txt, _qm_mm.txt
        """
        if not self.lib:
            return False
        
        n_atoms = len(atoms)
        atom_array = (ExportAtom * n_atoms)()
        
        for i, a in enumerate(atoms):
            atom_array[i].id = a['id']
            atom_array[i].type_id = a.get('type_id', 1)
            atom_array[i].molecule = a.get('molecule', 1)
            atom_array[i].element = a['element'].ljust(2)[:2].encode('utf-8') + b'\x00'
            atom_array[i].x = a['x']
            atom_array[i].y = a['y']
            atom_array[i].z = a['z']
            atom_array[i].charge = a.get('charge', 0.0)
            atom_array[i].is_qm = 1 if a.get('is_qm', False) else 0
        
        if projection_point:
            proj = (ctypes.c_double * 3)(projection_point[0],
                                          projection_point[1],
                                          projection_point[2])
        else:
            max_z = max(a['z'] for a in atoms)
            proj = (ctypes.c_double * 3)(0.0, 0.0, max_z + 10.0)
        
        filename_bytes = str(base_filename).encode('utf-8')
        self.lib.export_schlegel_only(filename_bytes, n_atoms, atom_array, proj, view_angle)
        
        print(f"✅ Exported Schlegel projections to {base_filename}_*.txt")
        return True
    
    def get_schlegel_coordinates(self, atoms: List[Dict],
                                 projection_point: Optional[Tuple[float, float, float]] = None,
                                 view_angle: float = 45.0) -> List[Dict]:
        """
        Calculate Schlegel projection coordinates without exporting files
        
        Returns atoms with added 'schlegel' coordinate field
        """
        n_atoms = len(atoms)
        atom_array = (ExportAtom * n_atoms)()
        
        for i, a in enumerate(atoms):
            atom_array[i].id = a['id']
            atom_array[i].type_id = a.get('type_id', 1)
            atom_array[i].molecule = a.get('molecule', 1)
            atom_array[i].element = a['element'].ljust(2)[:2].encode('utf-8') + b'\x00'
            atom_array[i].x = a['x']
            atom_array[i].y = a['y']
            atom_array[i].z = a['z']
            atom_array[i].charge = a.get('charge', 0.0)
            atom_array[i].is_qm = 1 if a.get('is_qm', False) else 0
        
        if projection_point:
            proj = (ctypes.c_double * 3)(projection_point[0],
                                          projection_point[1],
                                          projection_point[2])
        else:
            max_z = max(a['z'] for a in atoms)
            proj = (ctypes.c_double * 3)(0.0, 0.0, max_z + 10.0)
        
        # Call internal calculation (would need additional Fortran function)
        # For now, return atoms with placeholder
        result = atoms.copy()
        for i, a in enumerate(result):
            a['schlegel_xy'] = (a['x'], a['y'])  # Placeholder
            a['schlegel_xz'] = (a['x'], a['z'])
            a['schlegel_yz'] = (a['y'], a['z'])
        
        return result


# Singleton instance
_exporter = None

def get_molecular_exporter() -> Optional[MolecularExporter]:
    """Get or create MolecularExporter singleton"""
    global _exporter
    if _exporter is None and lib is not None:
        _exporter = MolecularExporter()
    return _exporter

# Test function
def test_exporter():
    """Test the molecular exporter with sample data"""
    exporter = get_molecular_exporter()
    if not exporter:
        print("[Error] Exporter not available")
        return
    
    # Sample methane molecule
    atoms = [
        {'id': 1, 'element': 'C', 'x': 0.000, 'y': 0.000, 'z': 0.000, 'is_qm': True},
        {'id': 2, 'element': 'H', 'x': 0.628, 'y': 0.628, 'z': 0.628, 'is_qm': True},
        {'id': 3, 'element': 'H', 'x': -0.628, 'y': -0.628, 'z': 0.628, 'is_qm': True},
        {'id': 4, 'element': 'H', 'x': -0.628, 'y': 0.628, 'z': -0.628, 'is_qm': True},
        {'id': 5, 'element': 'H', 'x': 0.628, 'y': -0.628, 'z': -0.628, 'is_qm': True},
    ]
    
    bonds = [
        {'id': 1, 'atom1': 1, 'atom2': 2, 'order': 1},
        {'id': 2, 'atom1': 1, 'atom2': 3, 'order': 1},
        {'id': 3, 'atom1': 1, 'atom2': 4, 'order': 1},
        {'id': 4, 'atom1': 1, 'atom2': 5, 'order': 1},
    ]
    
    # Export all formats
    exporter.export_all("test_molecule", atoms, bonds, view_angle=45.0)
    
    # Export just Schlegel projections
    exporter.export_schlegel_only("test_schlegel", atoms, view_angle=45.0)
    
    print("\n✅ Test complete - check generated files:")
    print("   test_molecule.xyz, .pdb, .sdf, .mol2, .cif, .gro")
    print("   test_schlegel_xy.txt, _xz.txt, _yz.txt, _schlegel.txt, _qm_mm.txt")

if __name__ == "__main__":
    test_exporter()