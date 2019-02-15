//
// auto-generated by mm.py
//
#include <iostream>
#include "full_matrix/Computation.hpp"
#include "Arguments.hpp"
#include "Data.hpp"
#include "IndexGenerator.hpp"

using namespace MM;
using namespace MM::full_matrix;

#define NeighProxy NeighProxyDirect

//user function
void anonymusAt1885(NeighProxy<CellData<2>> x, NeighProxy<CellData<2>> y,
                          NeighProxy<CellMatData<2>> Vf, NeighProxy<CellMatData<2>> rho,
                          double &rho_out) {
			double xo = x[{0,0}];
			double yo = y[{0,0}];

			double dsqr[9];

			for (int nj = -1; nj <= 1; nj++) {
				for (int ni = -1; ni <= 1; ni++) {

					dsqr[(nj+1)*3 + (ni+1)] = 0.0;

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
    }

// host stub function
void mm_par_loop_anonymusAt1885(std::string name, Computation<2> &computation,
 NEIGH<CellData<2>> arg0, NEIGH<CellData<2>> arg1, NEIGH<CellMatData<2>> arg2, NEIGH<CellMatData<2>> arg3,
 OUT<CellMatData<2>> arg4) {
  const std::array<std::size_t, 2> &begin = computation.index_generator.get_begin();
  const std::array<std::size_t, 2> &end = computation.index_generator.get_end();
  std::array<std::size_t,2> shape0;
  const double * __restrict data0 = arg0.get_raw(shape0);
  std::array<std::size_t,2> shape1;
  const double * __restrict data1 = arg1.get_raw(shape1);
  std::array<std::size_t,2> shape2;
  const double * __restrict data2 = arg2.get_raw(shape2);
  std::array<std::size_t,2> shape3;
  const double * __restrict data3 = arg3.get_raw(shape3);
  std::array<std::size_t,2> shape4;
  double * __restrict data4 = arg4.get_raw(shape4);
  for (std::size_t mat_index = 0; mat_index < computation.data.get_mat_number(); mat_index++) {
    #pragma omp parallel for //private(rho_sum) //collapse(1)
    for (std::size_t j = begin[1]; j < end[1]; j++) {
      #pragma omp simd //private(rho_sum)
      for (std::size_t i = begin[0]; i < end[0]; i++) {
        const Coords<2> coords(i, j);
/*        anonymusAt1885(
          NeighProxy<CellData<2>>(shape0, &data0[i+j*shape0[0]]),
          NeighProxy<CellData<2>>(shape1, &data1[i+j*shape1[0]]),
          NeighProxy<CellMatData<2>>(shape2, &data2[i+j*shape2[0]+mat_index*shape2[0]*shape2[1]]),
          NeighProxy<CellMatData<2>>(shape3, &data3[i+j*shape3[0]+mat_index*shape3[0]*shape3[1]]),
          data4[i+j*shape4[0]+mat_index*shape4[0]*shape4[1]]);*/

          NeighProxy<CellData<2>> x(shape0, &data0[i+j*shape0[0]]);
          NeighProxy<CellData<2>> y(shape1, &data1[i+j*shape1[0]]);
          NeighProxy<CellMatData<2>> Vf(shape2, &data2[i+j*shape2[0]+mat_index*shape2[0]*shape2[1]]);
          NeighProxy<CellMatData<2>> rho(shape3, &data3[i+j*shape3[0]+mat_index*shape3[0]*shape3[1]]);
          //double &rho_out = data4[i+j*shape4[0]+mat_index*shape4[0]*shape4[1]];
       if (data2[i+j*shape2[0]] > 0.0) {
         // o: outer
         double xo = data0[i+j*shape0[0]];
         double yo = data1[i+j*shape1[0]];

         double rho_sum = 0.0;
         int Nn = 0;

         for (int nj = -1; nj <= 1; nj++) {
           for (int ni = -1; ni <= 1; ni++) {
             if (data2[(i+ni)+(j+nj)*shape2[0]] > 0.0) {
               double dsqr = 0.0;

               // i: inner
               double xi = data0[(i+ni)+(j+nj)*shape0[0]];
               double yi = data1[(i+ni)+(j+nj)*shape1[0]];

               dsqr += (xo - xi) * (xo - xi);
               dsqr += (yo - yi) * (yo - yi);

               rho_sum += data3[(i+ni)+(j+nj)*shape3[0]] / dsqr;
               Nn += 1;
             }
           }
         }

         data4[i+j*shape4[0]+mat_index*shape4[0]*shape4[1]] = rho_sum / Nn;
       }
       else {
         data4[i+j*shape4[0]+mat_index*shape4[0]*shape4[1]] = 0.0;
       }
    /* 
			double xo = x[{0,0}];
			double yo = y[{0,0}];

			double dsqr[9];

			for (int nj = -1; nj <= 1; nj++) {
				for (int ni = -1; ni <= 1; ni++) {

					dsqr[(nj+1)*3 + (ni+1)] = 0.0;

					double xi = x[{ni,nj}];
					double yi = y[{ni,nj}];

					dsqr[(nj+1)*3 + (ni+1)] += (xo - xi) * (xo - xi);
					dsqr[(nj+1)*3 + (ni+1)] += (yo - yi) * (yo - yi);
				}
			}
      if (Vf[{0,0}] > 0.0) {
        int Nn = 0;
        double rho_sum = 0.0;

        for (int nj = -1; nj <= 1; nj++) {
          for (int ni = -1; ni <= 1; ni++) {

            if (Vf[{ni,nj}] > 0.0) {
              rho_sum += rho[{ni,nj}] / dsqr[(nj+1)*3 + (ni+1)];
              Nn += 1;
            }
          }
        }
        data4[i+j*shape4[0]+mat_index*shape4[0]*shape4[1]] = rho_sum / Nn;
      } else
        data4[i+j*shape4[0]+mat_index*shape4[0]*shape4[1]] = 0.0;
    */
      }
    }
  }
}