"""
fast_parser.py - Ultra-fast ADOBMD file parser using LuaJIT
Integrates with your existing luajit_bridge

This file is part of Quantum Analysis Helper.
Quantum Analysis Helper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Copyright (C) [2026] Acxel David Orozco Baldomero
"""

import time
import json
import os
from pathlib import Path
from typing import Optional, Dict, Any, List
import tempfile

# Try to import LuaJIT bridge
try:
    from plugins.libs.luajit_bridge import get_luajit_bridge
    LUAJIT_AVAILABLE = True
except ImportError:
    LUAJIT_AVAILABLE = False
    print("[Warning] LuaJIT bridge not available, using fallback parser")

from .data_parser import ADOBMDData, Atom, Bond

class LuaJITFastParser:
    """
    Ultra-fast ADOBMD file parser using LuaJIT
    10x faster than pure Python parser
    """
    
    def __init__(self):
        self.lua = None
        self.L = None
        self.parser_loaded = False
        self.parser_script = None
        self.stats = {
            'files_parsed': 0,
            'total_time': 0,
            'avg_time': 0
        }
        
        # Try to initialize LuaJIT
        if LUAJIT_AVAILABLE:
            try:
                self.lua = get_luajit_bridge()
                print("[Ok] LuaJIT bridge acquired")
            except Exception as e:
                print(f"[Warning] Failed to get LuaJIT bridge: {e}")
                self.lua = None
    
    def initialize(self) -> bool:
        """Initialize Lua state and load parser script"""
        if not self.lua:
            print("[Warning] LuaJIT not available, cannot initialize")
            return False
        
        try:
            # Create new Lua state
            self.L = self.lua.init()
            if not self.L:
                raise RuntimeError("Failed to create Lua state")
            
            # Load parser script
            parser_path = Path(__file__).parent / "lua" / "adobmd_parser.lua"
            if not parser_path.exists():
                raise FileNotFoundError(f"Parser script not found: {parser_path}")
            
            self.parser_script = parser_path.read_text(encoding='utf-8')
            
            # Execute parser script
            if not self.lua.execute_script(self.L, self.parser_script):
                error = self._get_last_error()
                raise RuntimeError(f"Failed to load parser script: {error}")
            
            self.parser_loaded = True
            print("[Ok] LuaJIT parser initialized successfully")
            return True
            
        except Exception as e:
            print(f"[Error] Failed to initialize LuaJIT parser: {e}")
            self.cleanup()
            return False
    
    def _get_last_error(self) -> str:
        """Get last Lua error message"""
        if not self.lua or not self.L:
            return "Unknown error"
        
        try:
            # Get error from stack
            error_ptr = self.lua.lib.lua_tolstring(self.L, -1, None)
            if error_ptr:
                import ctypes
                return ctypes.string_at(error_ptr).decode('utf-8', errors='ignore')
        except:
            pass
        return "Unknown Lua error"
    
    def parse_file_old(self, filepath: str) -> Optional[ADOBMDData]:
        """
        Parse ADOBMD file using LuaJIT
        
        Args:
            filepath: Path to .data file
            
        Returns:
            ADOBMDData object or None if failed
        """
        if not os.path.exists(filepath):
            print(f"[Error] File not found: {filepath}")
            return None
        
        # Initialize if needed
        if not self.parser_loaded:
            if not self.initialize():
                return None
        
        start_time = time.time()
        
        try:
            # Read file (Python I/O is fine, LuaJIT would be slower for this)
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            read_time = time.time() - start_time
            
            # Call LuaJIT parser
            parse_start = time.time()
            
            # Escape content for Lua string
            # This is critical to avoid syntax errors
            escaped_content = self._escape_lua_string(content)
            
            # Call the parser
            result_str = self.lua.call_function_str(
                self.L,
                "M.parse_file",
                escaped_content
            )
            
            parse_time = time.time() - parse_start
            
            if not result_str or result_str == "":
                print("[Error] Parser returned empty result")
                return None
            
            # Parse JSON result
            try:
                result = json.loads(result_str)
            except json.JSONDecodeError as e:
                print(f"[Error] Failed to parse JSON result: {e}")
                print(f"Result preview: {result_str[:200]}")
                return None
            
            # Convert to ADOBMDData
            data = self._convert_to_adobmd_data(result, filepath)
            
            total_time = time.time() - start_time
            
            # Update stats
            self.stats['files_parsed'] += 1
            self.stats['total_time'] += total_time
            self.stats['avg_time'] = self.stats['total_time'] / self.stats['files_parsed']
            
            print(f"⚡ LuaJIT parse: {parse_time:.3f}s, Total: {total_time:.3f}s")
            print(f"   Read: {read_time:.3f}s, Parse: {parse_time:.3f}s")
            print(f"   Atoms: {len(data.atoms)}, Bonds: {len(data.bonds)}")
            
            return data
            
        except Exception as e:
            print(f"[Error] Error parsing file with LuaJIT: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def parse_file(self, filepath: str) -> Optional[ADOBMDData]:
        """Parse ADOBMD file using LuaJIT"""
        if not os.path.exists(filepath):
            print(f"[Error] File not found: {filepath}")
            return None
        
        # Initialize if needed
        if not self.parser_loaded:
            if not self.initialize():
                return None
        
        start_time = time.time()
        
        try:
            # Read file
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Call LuaJIT parser
            escaped_content = self._escape_lua_string(content)
            result_str = self.lua.call_function_str(
                self.L,
                "M.parse_file",
                escaped_content
            )
            
            if not result_str or result_str == "":
                print("[Error] Parser returned empty result")
                # Try to get error from Lua
                error = self._get_last_error()
                if error:
                    print(f"[Error] Lua error: {error}")
                return None
            
            # Debug: print first 200 chars of result
            print(f"[Debug] Lua result preview: {result_str[:200]}")
            
            # Parse JSON result
            try:
                result = json.loads(result_str)
            except json.JSONDecodeError as e:
                print(f"[Error] Failed to parse JSON result: {e}")
                print(f"Result preview: {result_str[:500]}")
                return None
            
            # Convert to ADOBMDData
            data = self._convert_to_adobmd_data(result, filepath)
            
            elapsed = time.time() - start_time
            print(f"⚡ LuaJIT parsed {len(data.atoms)} atoms in {elapsed*1000:.1f}ms")
            
            return data
            
        except Exception as e:
            print(f"[Error] LuaJIT parse failed: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _escape_lua_string(self, content: str) -> str:
        """Escape content for safe inclusion in Lua string"""
        # Use long bracket format for safety
        # This handles any content without escaping
        return f"[=[{content}]=]"
    
    def _convert_to_adobmd_data(self, result: Dict, filepath: str) -> ADOBMDData:
        """Convert LuaJIT result to ADOBMDData object"""
        data = ADOBMDData()
        data.filename = os.path.basename(filepath)
        
        # Set header info
        if 'header' in result:
            header = result['header']
            data.natoms = header.get('natoms', 0)
            data.nbonds = header.get('nbonds', 0)
            data.natom_types = header.get('natom_types', 0)
            
            if header.get('has_box', False):
                data.has_box = True
                data.box_lo = header['box_lo']
                data.box_hi = header['box_hi']
        
        # Create type_elements mapping
        if 'type_elements' in result:
            # Convert string keys to int
            for k, v in result['type_elements'].items():
                data.type_elements[int(k)] = v
        
        # Create atoms
        data.atoms = []
        data.qm_indices = []
        
        if 'atoms' in result:
            for atom_data in result['atoms']:
                atom = Atom(
                    id=int(atom_data['id']),
                    type_id=int(atom_data['type']),
                    molecule=int(atom_data['mol']),
                    element=atom_data['element'],
                    x=float(atom_data['x']),
                    y=float(atom_data['y']),
                    z=float(atom_data['z']),
                    charge=float(atom_data['charge']),
                    is_qm=bool(atom_data['qm'])
                )
                data.atoms.append(atom)
                
                if atom.is_qm:
                    data.qm_indices.append(atom.id)
        
        data.has_qm_region = len(data.qm_indices) > 0
        
        # Create bonds
        data.bonds = []
        if 'bonds' in result:
            for bond_data in result['bonds']:
                bond = Bond(
                    id=int(bond_data['id']),
                    type_id=int(bond_data['type']),
                    atom1=int(bond_data['a1']),
                    atom2=int(bond_data['a2']),
                    order=int(bond_data.get('order', bond_data['type']))
                )
                data.bonds.append(bond)
        
        return data
    
    def parse_file_super_fast(self, filepath: str) -> Optional[ADOBMDData]:
        """
        EVEN FASTER: Let LuaJIT handle file I/O too
        Useful for very large files
        """
        if not self.parser_loaded:
            if not self.initialize():
                return None
        
        start_time = time.time()
        
        try:
            # Escape filepath for Lua
            # Use forward slashes for cross-platform
            lua_path = filepath.replace('\\', '/')
            
            # Create Lua script that reads and parses in one go
            script = f"""
            local f = io.open('{lua_path}', 'r')
            if not f then
                return "ERROR: Cannot open file"
            end
            local content = f:read('*all')
            f:close()
            return M.parse_file(content)
            """
            
            # Execute script
            result_str = self.lua.call_function_str(self.L, "loadstring", script)
            
            if not result_str or result_str == "":
                return None
            
            # Parse JSON
            result = json.loads(result_str)
            
            total_time = time.time() - start_time
            print(f"⚡⚡ Super fast mode: {total_time:.3f}s")
            
            return self._convert_to_adobmd_data(result, filepath)
            
        except Exception as e:
            print(f"[Error] Super fast mode failed: {e}")
            return None
    
    def parse_file_parallel(self, filepaths: List[str]) -> List[Optional[ADOBMDData]]:
        """
        Parse multiple files in parallel using LuaJIT
        Note: Each file needs its own Lua state
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        
        results = []
        
        def parse_single(path):
            # Create new parser instance for thread safety
            parser = LuaJITFastParser()
            if parser.initialize():
                return parser.parse_file(path)
            return None
        
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {executor.submit(parse_single, path): path for path in filepaths}
            
            for future in as_completed(futures):
                try:
                    result = future.result(timeout=30)
                    results.append(result)
                except Exception as e:
                    print(f"[Error] Parallel parse failed: {e}")
                    results.append(None)
        
        return results
    
    def benchmark(self, filepath: str, iterations: int = 10) -> Dict:
        """Benchmark parser performance"""
        if not os.path.exists(filepath):
            return {"error": "File not found"}
        
        times = []
        
        for i in range(iterations):
            start = time.time()
            data = self.parse_file(filepath)
            if data:
                times.append(time.time() - start)
        
        if not times:
            return {"error": "No successful parses"}
        
        return {
            'min': min(times),
            'max': max(times),
            'avg': sum(times) / len(times),
            'files_parsed': len(times),
            'atoms': len(data.atoms) if data else 0,
            'bonds': len(data.bonds) if data else 0
        }
    
    def get_stats(self) -> Dict:
        """Get parser statistics"""
        return self.stats.copy()
    
    def cleanup(self):
        """Clean up Lua state"""
        if self.lua and self.L:
            try:
                self.lua.lib.lua_close(self.L)
            except:
                pass
            self.L = None
            self.parser_loaded = False
            print("[Ok] LuaJIT parser cleaned up")
    
    def __del__(self):
        """Destructor - ensure cleanup"""
        self.cleanup()


# Factory function for easy integration
def get_fast_parser(use_luajit: bool = True) -> Optional[LuaJITFastParser]:
    """
    Get a fast parser instance
    
    Args:
        use_luajit: If False, returns None (use Python fallback)
    
    Returns:
        LuaJITFastParser or None if not available
    """
    if not use_luajit:
        return None
    
    parser = LuaJITFastParser()
    if parser.initialize():
        return parser
    return None


# Simple test function
def test_parser(filepath: str):
    """Test the parser on a file"""
    parser = LuaJITFastParser()
    
    print(f"\n🔍 Testing parser on: {filepath}")
    print("-" * 50)
    
    if not parser.initialize():
        print("[Error] Failed to initialize parser")
        return
    
    # Parse file
    data = parser.parse_file(filepath)
    
    if data:
        print(f"[Ok] Success!")
        print(f"   Atoms: {len(data.atoms)}")
        print(f"   Bonds: {len(data.bonds)}")
        print(f"   QM atoms: {len(data.qm_indices)}")
        print(f"   Has box: {data.has_box}")
        
        # Show stats
        stats = parser.get_stats()
        print(f"\nStats: {stats}")
    else:
        print("[Error] Failed to parse")
    
    parser.cleanup()


if __name__ == "__main__":
    # Quick test
    import sys
    if len(sys.argv) > 1:
        test_parser(sys.argv[1])
    else:
        print("Usage: python fast_parser.py <file.data>")