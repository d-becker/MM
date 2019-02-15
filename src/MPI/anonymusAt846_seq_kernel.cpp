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
inline void anonymusAt846(double rho, double Vf, double V, double &rho_ave) {
      rho_ave += rho*Vf/V;
      }

// host stub function
void mm_par_loop_anonymusAt846(std::string name, Computation<2> &computation,
 IN<CellMatData<2>> arg0, IN<CellMatData<2>> arg1, IN<CellData<2>> arg2, REDUCE<CellData<2>> arg3) {
  const std::array<std::size_t, 2> &begin = computation.index_generator.get_begin();
  const std::array<std::size_t, 2> &end = computation.index_generator.get_end();
  std::array<std::size_t,2> shape0;
  const double * __restrict data0 = arg0.get_raw(shape0);
  std::array<std::size_t,2> shape1;
  const double * __restrict data1 = arg1.get_raw(shape1);
  std::array<std::size_t,2> shape2;
  double * __restrict data2 = arg2.get_raw(shape2);
  double rho_ave = 0;
  for (std::size_t mat_index = 0; mat_index < computation.data.get_mat_number(); mat_index++) {
    #pragma omp parallel for collapse(1) reduction(+:rho_ave)
    for (std::size_t j = begin[1]; j < end[1]; j++) {
      #pragma omp simd reduction(+:rho_ave)
      for (std::size_t i = begin[0]; i < end[0]; i++) {
        anonymusAt846(
          data0[i+j*shape0[0]+mat_index*shape0[0]*shape0[1]],
          data1[i+j*shape1[0]+mat_index*shape1[0]*shape1[1]],
          data2[i+j*shape2[0]],
          rho_ave);
      }
    }
  }
  const Coords<2> coords(0ul, 0ul);
  arg3.get(coords, 0) << rho_ave;
}