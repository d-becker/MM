#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <iostream>

#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Computation.hpp"
#include "full_matrix/Data.hpp"

using namespace std;
using namespace MM;
using namespace MM::compressed_cell_centric;

extern void compact_matrix_cell_centric(unsigned int sizex, unsigned int sizey, int Nmats, MM::full_matrix::Data<2> &d,
	MM::full_matrix::CellMatData<2>& f_rho, MM::full_matrix::CellMatData<2>& f_rho_mat_ave,
  MM::full_matrix::CellMatData<2>& f_p, MM::full_matrix::CellMatData<2>& f_Vf, MM::full_matrix::CellMatData<2>& f_t,
	MM::full_matrix::CellData<2>& f_V, MM::full_matrix::CellData<2>& f_x, MM::full_matrix::CellData<2>& f_y,
	MM::full_matrix::MatData<2>& f_n, MM::full_matrix::CellData<2>& f_rho_ave, vector<vector<size_t>> &mats)
{
  Data<2> data({sizex, sizey}, mats);

  CellMatData<2> rho = data.new_cell_mat_data();
  CellMatData<2> rho_mat_ave = data.new_cell_mat_data();
  CellMatData<2> p = data.new_cell_mat_data();
  CellMatData<2> Vf = data.new_cell_mat_data();
  CellMatData<2> t = data.new_cell_mat_data();

  CellData<2> V = data.new_cell_data();
  CellData<2> x = data.new_cell_data();
  CellData<2> y = data.new_cell_data();
  CellData<2> rho_ave = data.new_cell_data();

  MatData<2> n = data.new_mat_data();

  size_t compact_idx = 0;
  for (size_t j = 0; j <  sizey; j++) {
    for (size_t i = 0; i <  sizex; i++) {
      size_t idx = j * sizex + i;
      if (mats[idx].size() > 1) {
        for (size_t m = 0; m < mats[idx].size(); m++) {
          rho.mixed_storage_value_at(compact_idx) = f_rho.at(Coords<2>(i,j),mats[idx][m]);
          Vf.mixed_storage_value_at(compact_idx) = f_Vf.at(Coords<2>(i,j),mats[idx][m]);
          p.mixed_storage_value_at(compact_idx) = f_p.at(Coords<2>(i,j),mats[idx][m]);
          t.mixed_storage_value_at(compact_idx) = f_t.at(Coords<2>(i,j),mats[idx][m]);
          compact_idx++;
        }
      } else {
        rho.cell_value_at(idx) = f_rho.at(Coords<2>(i,j),mats[idx][0]);
        Vf.cell_value_at(idx) = f_Vf.at(Coords<2>(i,j),mats[idx][0]);
        p.cell_value_at(idx) = f_p.at(Coords<2>(i,j),mats[idx][0]);
        t.cell_value_at(idx) = f_t.at(Coords<2>(i,j),mats[idx][0]);
      }
        V.at(Coords<2>(i,j)) = f_V.at(Coords<2>(i,j));
        x.at(Coords<2>(i,j)) = f_x.at(Coords<2>(i,j));
        y.at(Coords<2>(i,j)) = f_y.at(Coords<2>(i,j));
    }
  }
  for (size_t mat = 0; mat < Nmats; mat++) {
    n.at(mat) = f_n.at(mat);
  }

	// Cell-centric algorithms
	// Computational loop 1 - average density in cell
  double t1 = omp_get_wtime();
  auto INC = [] (double left, double right) {
    return left + right;
  };
  IndexGenerator<2> index_generator({0, 0}, {sizex, sizey});
  Computation<2> computation(data, index_generator);
  computation.compute([] (double rho, double Vf, double V, ReduceProxy<double> rho_ave) {
      rho_ave << rho*Vf/V;
      },
      IN<CellMatData<2>>(rho),
      IN<CellMatData<2>>(Vf),
      IN<CellData<2>>(V),
      REDUCE<CellData<2>>(INC, rho_ave));

  printf("Compact matrix, cell centric, alg 1: %g sec\n", omp_get_wtime()-t1);

	// Computational loop 2 - Pressure for each cell and each material
  t1 = omp_get_wtime();
  Computation<2> computation3(data, index_generator);
  computation3.compute([](double &p, double n, double rho, double t, double Vf) {
      if (Vf > 0) p = (n * rho * t) / Vf;
      else p = 0;
      },
      OUT<CellMatData<2>>(p),
      IN<MatData<2>>(n),
      IN<CellMatData<2>>(rho),
      IN<CellMatData<2>>(t),
      IN<CellMatData<2>>(Vf));
  printf("Compact matrix, cell centric, alg 2: %g sec\n", omp_get_wtime()-t1);

	// Computational loop 3 - Average density of each material over neighborhood of each cell
  t1 = omp_get_wtime();
  Stencil<2> s9pt({{1,1},  {1,0},  {1,-1},
       {0,1},  {0,0},  {0,-1},
       {-1,1}, {-1,0}, {-1,-1}});
  IndexGenerator<2> index_generator2({1, 1}, {sizex-1, sizey-1});
  Computation<2> computation2(data, index_generator2);
  computation2.compute([](NeighProxy<CellData<2>> x, NeighProxy<CellData<2>> y,
                          NeighProxy<CellMatData<2>> Vf, NeighProxy<CellMatData<2>> rho,
                          double &rho_out) {
			double xo = x[{0,0}];
			double yo = y[{0,0}];

			// There are at most 9 neighbours in 2D case.
			double dsqr[9];

			for (int nj = -1; nj <= 1; nj++) {
				for (int ni = -1; ni <= 1; ni++) {

					dsqr[(nj+1)*3 + (ni+1)] = 0.0;

					// i: inner
					double xi = x[{ni,nj}];
					double yi = y[{ni,nj}];

					dsqr[(nj+1)*3 + (ni+1)] += (xo - xi) * (xo - xi);
					dsqr[(nj+1)*3 + (ni+1)] += (yo - yi) * (yo - yi);
				}
			}
      if (Vf[{0,0}] > 0.0) {
        double rho_sum = 0.0;
        int Nn = 0;

        for (int nj = -1; nj <= 1; nj++) {
          for (int ni = -1; ni <= 1; ni++) {

            if (Vf.has_neigh({ni,nj})) {
              rho_sum += rho[{ni,nj}] / dsqr[(nj+1)*3 + (ni+1)];
              Nn += 1;
            }
          }
        }
        rho_out = rho_sum / Nn;
      } else
        rho_out = 0.0;
    },
    NEIGH<CellData<2>>(x, s9pt),
    NEIGH<CellData<2>>(y, s9pt),
    NEIGH<CellMatData<2>>(Vf, s9pt),
    NEIGH<CellMatData<2>>(rho, s9pt),
    OUT<CellMatData<2>>(rho_mat_ave));
  printf("Compact matrix, cell centric, alg 3: %g sec\n", omp_get_wtime()-t1);
}
