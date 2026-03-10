# Linux
gfortran -shared -fPIC -O3 -o libadobmd_parser.so libadobmd_parser.f90
gfortran -shared -fPIC -O3 -o libmolexport.so molecular_export_full.f90

# macOS
gfortran -shared -O3 -o libadobmd_parser.dylib libadobmd_parser.f90

# Windows (with MinGW)
gfortran -shared -O3 -o libadobmd_parser.dll libadobmd_parser.f90

