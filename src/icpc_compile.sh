python2 ../translator/ops.py full_matrix.cpp
python2 ../translator/ops.py compressed.cpp
icpc -O3 -xHost -Ifull_matrix -Icompressed_cell_centric compressed_mm.cpp full_matrix_mm.cpp MPI/*.cpp multimat_abs.cpp -o multimat_abs -std=c++17 -I. -qopenmp -qopt-report=5
