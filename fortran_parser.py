"""
fortran_parser.py - Ultra-fast ADOBMD file parser using Fortran library
Fixed memory management issues
"""

import ctypes
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import platform
import time
import os

# Platform-specific library extension
if platform.system() == 'Windows':
    LIB_EXT = '.dll'
elif platform.system() == 'Darwin':
    LIB_EXT = '.dylib'
else:  # Linux
    LIB_EXT = '.so'

# Load Fortran library
lib_path = Path(__file__).parent / f"libadobmd_parser{LIB_EXT}"
if lib_path.exists():
    try:
        lib = ctypes.CDLL(str(lib_path))
        print(f"[Ok] Loaded Fortran library: {lib_path}")
    except Exception as e:
        lib = None
        print(f"[Error] Failed to load Fortran library: {e}")
else:
    lib = None
    print(f"[Warning] Fortran library not found: {lib_path}")

class Atom(ctypes.Structure):
    """C-compatible atom structure"""
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
    ]

class Bond(ctypes.Structure):
    """C-compatible bond structure"""
    _fields_ = [
        ('id', ctypes.c_int),
        ('type_id', ctypes.c_int),
        ('atom1', ctypes.c_int),
        ('atom2', ctypes.c_int),
        ('order', ctypes.c_int),
    ]

class ParseResult(ctypes.Structure):
    """C-compatible parse result structure"""
    _fields_ = [
        ('natoms', ctypes.c_int),
        ('nbonds', ctypes.c_int),
        ('nqm', ctypes.c_int),
        ('box_lo', ctypes.c_double * 3),
        ('box_hi', ctypes.c_double * 3),
        ('has_box', ctypes.c_int),
        ('atoms', ctypes.POINTER(Atom)),
        ('bonds', ctypes.POINTER(Bond)),
        ('filename', ctypes.c_char * 256)
    ]

class FortranParser:
    """
    Ultra-fast parser using Fortran library
    Fixed memory management - proper cleanup
    """
    
    def __init__(self):
        self.lib = lib
        self.version = "Unknown"
        self._result_cache = None  # Store result to prevent premature GC
        
        if self.lib:
            # Configure parse_file function
            self.lib.parse_file.argtypes = [ctypes.c_char_p, ctypes.POINTER(ParseResult)]
            self.lib.parse_file.restype = None
            
            # Configure free_result function
            self.lib.free_result.argtypes = [ctypes.POINTER(ParseResult)]
            self.lib.free_result.restype = None
            
            # Get version
            try:
                self.lib.get_version.argtypes = [ctypes.c_char_p]
                self.lib.get_version.restype = None
                version_buf = ctypes.create_string_buffer(256)
                self.lib.get_version(version_buf)
                self.version = version_buf.value.decode('utf-8')
                print(f"[Ok] Fortran parser version: {self.version}")
            except:
                print("[Ok] Fortran parser loaded (version unknown)")
    
    def parse_file(self, filepath: str) -> Optional[Dict]:
        """
        Parse file using Fortran library
        
        Args:
            filepath: Path to .data file
            
        Returns:
            Dictionary with atoms, bonds, and metadata
        """
        if not self.lib:
            print("[Error] Fortran library not loaded")
            return None
        
        if not os.path.exists(filepath):
            print(f"[Error] File not found: {filepath}")
            return None
        
        start_time = time.time()
        
        # Create result structure
        result = ParseResult()
        filename_bytes = str(filepath).encode('utf-8')
        
        try:
            # Call Fortran parser
            self.lib.parse_file(filename_bytes, ctypes.byref(result))
            
            # Check for error
            if result.natoms < 0:
                print(f"[Error] Fortran parser failed for: {filepath}")
                return None
            
            # Extract atoms
            atoms = []
            for i in range(result.natoms):
                atom = result.atoms[i]
                atoms.append({
                    'id': atom.id,
                    'type_id': atom.type_id,
                    'molecule': atom.molecule,
                    'element': atom.element.decode('utf-8').strip('\x00'),
                    'x': atom.x,
                    'y': atom.y,
                    'z': atom.z,
                    'charge': atom.charge,
                    'is_qm': bool(atom.is_qm)
                })
            
            # Extract bonds
            bonds = []
            for i in range(result.nbonds):
                bond = result.bonds[i]
                bonds.append({
                    'id': bond.id,
                    'type_id': bond.type_id,
                    'atom1': bond.atom1,
                    'atom2': bond.atom2,
                    'order': bond.order
                })
            
            # Build result
            data = {
                'atoms': atoms,
                'bonds': bonds,
                'qm_indices': [a['id'] for a in atoms if a['is_qm']],
                'header': {
                    'has_box': bool(result.has_box),
                    'box_lo': [result.box_lo[i] for i in range(3)],
                    'box_hi': [result.box_hi[i] for i in range(3)],
                    'natoms': result.natoms,
                    'nbonds': result.nbonds
                },
                'has_qm': result.nqm > 0,
                'filename': result.filename.decode('utf-8').strip('\x00')
            }
            
            elapsed = time.time() - start_time
            print(f"⚡ Fortran parsed {result.natoms} atoms in {elapsed*1000:.1f}ms")
            
            # IMPORTANT: Store result to prevent premature garbage collection
            self._result_cache = result
            
            return data
            
        except Exception as e:
            print(f"[Error] Fortran parser exception: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def cleanup(self):
        """Explicitly clean up cached results"""
        self._result_cache = None
    
    def __del__(self):
        """Destructor - ensure cleanup"""
        self.cleanup()


# Singleton instance
_parser = None

def get_fortran_parser() -> Optional[FortranParser]:
    """Get or create Fortran parser singleton"""
    global _parser
    if _parser is None and lib is not None:
        _parser = FortranParser()
    return _parser