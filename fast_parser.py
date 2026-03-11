"""
fast_parser.py - Using C library directly
"""

import ctypes
import os
from pathlib import Path
from typing import Optional, Tuple
from .data_parser import ADOBMDData, Atom, Bond

# Load the C library
lib_path = Path(__file__).parent / "c" / "libadobmd_full.so"
lib = ctypes.CDLL(str(lib_path))

# Configure function signatures
lib.parse_adobmd_file.argtypes = [ctypes.c_char_p]
lib.parse_adobmd_file.restype = ctypes.c_void_p

lib.free_molecular_data.argtypes = [ctypes.c_void_p]
lib.free_molecular_data.restype = None

lib.export_all_formats.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
lib.export_all_formats.restype = None

lib.set_schlegel_params.argtypes = [ctypes.c_void_p, 
                                     ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                     ctypes.c_double, ctypes.c_double]
lib.set_schlegel_params.restype = None

lib.get_atom_data.argtypes = [ctypes.c_void_p, ctypes.c_int,
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.c_char_p,
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_double)]
lib.get_atom_data.restype = None

class FastParser:
    """Fast parser using C library"""
    
    def __init__(self):
        self._handle = None
    
    def parse_file(self, filepath: str) -> Optional[ADOBMDData]:
        """Parse ADOBMD file"""
        if not os.path.exists(filepath):
            print(f"[Error] File not found: {filepath}")
            return None
        
        # Parse with C library
        self._handle = lib.parse_adobmd_file(filepath.encode('utf-8'))
        if not self._handle:
            print("[Error] Parse failed")
            return None
        
        # We need to get the data from C
        # For now, we'll create a simple wrapper
        data = ADOBMDData()
        data.filename = os.path.basename(filepath)
        
        # In a real implementation, you'd add functions to get counts etc.
        # For now, let's create a simple approach
        
        return data
    
    def export_schlegel(self, base_filename: str, 
                        proj_x=0.0, proj_y=0.0, proj_z=10.0,
                        view_angle=45.0, scale=1.0):
        """Export Schlegel projections"""
        if not self._handle:
            print("[Error] No data loaded")
            return
        
        # Set parameters
        lib.set_schlegel_params(self._handle, proj_x, proj_y, proj_z, view_angle, scale)
        
        # Export
        lib.export_all_formats(self._handle, base_filename.encode('utf-8'))
        print(f"[Ok] Exported Schlegel projections to {base_filename}_*.txt")
    
    def cleanup(self):
        """Free C memory"""
        if self._handle:
            lib.free_molecular_data(self._handle)
            self._handle = None
    
    def __del__(self):
        self.cleanup()


# Simple test
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        parser = FastParser()
        data = parser.parse_file(sys.argv[1])
        if data:
            print("[Ok] Parse successful")
            # Export with Schlegel projections
            base = sys.argv[1].replace('.data', '')
            parser.export_schlegel(base)
        parser.cleanup()
    else:
        print("Usage: python fast_parser.py <file.data>")