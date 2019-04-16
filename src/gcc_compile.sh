python2 ../translator/ops.py full_matrix.cpp
python2 ../translator/ops.py compressed.cpp
g++ -g -O3 -Ifull_matrix -Icompressed_cell_centric compact_alg_1.cpp compact_alg_2.cpp compact_alg_3.cpp compressed_mm.cpp full_matrix_mm.cpp MPI/*.cpp multimat_abs.cpp -o multimat_abs -std=c++11 -I. -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -fopt-info-vec -fopt-info-vec-missed -march=native -DOMP -DNDEBUG
