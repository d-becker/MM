python2 ../translator/ops.py full_matrix.cpp
python2 ../translator/ops.py compressed.cpp
icpc -O3 -xHost -Ifull_matrix -Icompressed_cell_centric compact_alg_1.cpp compact_alg_2.cpp compact_alg_3.cpp compressed_mm.cpp full_matrix_mm.cpp MPI/*.cpp multimat_abs.cpp -o multimat_abs -std=c++17 -I. -qopenmp -qopt-report=5
