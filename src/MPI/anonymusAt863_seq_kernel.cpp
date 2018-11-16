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

template <typename dtype = double>
class ReduceProxyDirectAdd {
public:
	ReduceProxyDirectAdd(dtype& p_reduced_value)
		: reduced_value(p_reduced_value)
		{
		}
		
	void operator<<(dtype value) {
		reduced_value += value;
	}

private:		
	dtype& reduced_value;
};



#define NeighProxy NeighProxyDirect
#define ReduceProxy ReduceProxyDirectAdd

//user function
inline void anonymusAt863(double rho, double Vf, double V, ReduceProxy<double> rho_ave) {
      rho_ave << rho*Vf/V;
      }

// host stub function
void mm_par_loop_anonymusAt863(std::string name, Computation<2> &computation,
 IN<CellMatData<2>> arg0, IN<CellMatData<2>> arg1, IN<CellData<2>> arg2, REDUCE<CellData<2>> arg3) {
  const std::array<std::size_t, 2> &begin = computation.index_generator.get_begin();
  const std::array<std::size_t, 2> &end = computation.index_generator.get_end();
  std::array<std::size_t,2> shape0;
  std::vector<double *> data0 = arg0.get_raw(shape0);
  std::array<std::size_t,2> shape1;
  std::vector<double *> data1 = arg1.get_raw(shape1);
  std::array<std::size_t,2> shape2;
  std::vector<double *> data2 = arg2.get_raw(shape2);
  std::array<std::size_t,2> shape3;
  std::vector<double *> data3 = arg3.get_raw(shape3);
  for (std::size_t mat_index = 0; mat_index < computation.data.get_mat_number(); mat_index++) {
#pragma omp parallel for
    for (std::size_t j = begin[1]; j < end[1]; j++) {
#pragma omp simd
      for (std::size_t i = begin[0]; i < end[0]; i++) {
        const Coords<2> coords(i, j);
        anonymusAt863(
          data0[mat_index][i+j*shape0[0]],
          data1[mat_index][i+j*shape1[0]],
          data2[0][i+j*shape2[0]],
          ReduceProxyDirectAdd<double>(data3[0][i+j*shape3[0]]));
      }
    }
  }
}
