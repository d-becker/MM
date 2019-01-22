#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <cmath>
#include <iostream>

#include "Arguments.hpp"
#include "Computation.hpp"
#include "Coords.hpp"
#include "Data.hpp"
#include "Datasets.hpp"
#include "IndexGenerator.hpp"

using namespace std;
using namespace MM;
using namespace MM::full_matrix;

void full_matrix_cell_centric(unsigned int sizex, unsigned int sizey, int Nmats, Data<2> &data,
	CellMatData<2>& rho, CellMatData<2>& rho_mat_ave, CellMatData<2>& p, CellMatData<2>& Vf, CellMatData<2>& t,
	CellData<2>& V, CellData<2>& x, CellData<2>& y,
	MatData<2>& n, CellData<2> &rho_ave)
{
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

  printf("Full matrix, material centric, alg 1: %g sec\n", omp_get_wtime()-t1);

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
  printf("Full matrix, material centric, alg 2: %g sec\n", omp_get_wtime()-t1);

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

            if (Vf[{ni,nj}] > 0.0) {
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
  printf("Full matrix, material centric, alg 3: %g sec\n", omp_get_wtime()-t1);
}
