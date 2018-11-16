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
inline void anonymusAt1275(double &p, double n, double rho, double t, double Vf) {
      if (Vf > 0) p = (n * rho * t) / Vf;
      else p = 0;
      }

// host stub function
void mm_par_loop_anonymusAt1275(std::string name, Computation<2> &computation,
 OUT<CellMatData<2>> arg0, IN<MatData<>> arg1, IN<CellMatData<2>> arg2, IN<CellMatData<2>> arg3,
 IN<CellMatData<2>> arg4) {
  const std::array<std::size_t, 2> &begin = computation.index_generator.get_begin();
  const std::array<std::size_t, 2> &end = computation.index_generator.get_end();
  std::array<std::size_t,2> shape0;
  std::vector<double *> data0 = arg0.get_raw(shape0);
  std::array<std::size_t,1> shape1;
  std::vector<double *> data1 = arg1.get_raw(shape1);
  std::array<std::size_t,2> shape2;
  std::vector<double *> data2 = arg2.get_raw(shape2);
  std::array<std::size_t,2> shape3;
  std::vector<double *> data3 = arg3.get_raw(shape3);
  std::array<std::size_t,2> shape4;
  std::vector<double *> data4 = arg4.get_raw(shape4);
  for (std::size_t mat_index = 0; mat_index < computation.data.get_mat_number(); mat_index++) {
#pragma omp parallel for collapse(1)
    for (std::size_t j = begin[1]; j < end[1]; j++) {
      for (std::size_t i = begin[0]; i < end[0]; i++) {
        const Coords<2> coords(i, j);
        anonymusAt1275(
          data0[mat_index][i+j*shape0[0]],
          data1[0][mat_index],
          data2[mat_index][i+j*shape2[0]],
          data3[mat_index][i+j*shape3[0]],
          data4[mat_index][i+j*shape4[0]]);
      }
    }
  }
}
