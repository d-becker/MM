//
// auto-generated by mm.py
//
#include <iostream>
#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Computation.hpp"

using namespace MM;
using namespace MM::compressed_cell_centric;


//user function
inline void anonymusAt3831(NeighProxy<CellData<2>> x, NeighProxy<CellData<2>> y,
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

            if (Vf.has_neigh({ni,nj})) {
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
void mm_par_loop_anonymusAt3831(std::string name, Computation<2> &computation,
 NEIGH<CellData<2>> arg0, NEIGH<CellData<2>> arg1, NEIGH<CellMatData<2>> arg2, NEIGH<CellMatData<2>> arg3,
 OUT<CellMatData<2>> arg4) {
  const std::array<std::size_t, 2> &begin = computation.index_generator.get_begin();
  const std::array<std::size_t, 2> &end = computation.index_generator.get_end();

  CompressedDataStructure &structure = computation.data.structure;
  std::array<std::size_t,2> shape = computation.data.get_size();
  double * __restrict data0 = arg0.get_raw();
  double * __restrict data1 = arg1.get_raw();
  double * __restrict data2 = arg2.get_raw();
  double *__restrict data2_list = arg2.get_raw_list();
  double * __restrict data3 = arg3.get_raw();
  double *__restrict data3_list = arg3.get_raw_list();
  double * __restrict data4 = arg4.get_raw();
  double *__restrict data4_list = arg4.get_raw_list();
  #pragma omp parallel for collapse(1)
  for (std::size_t j = begin[1]; j < end[1]; j++) {
    #pragma omp simd
    for (std::size_t i = begin[0]; i < end[0]; i++) {
      if (structure.structure[i+j*shape[0]].nmats == 1) {
        std::size_t mat_index = structure.structure[i+j*shape[0]].imat;
        anonymusAt3831(
          NeighProxy<CellData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::SINGLE_MAT, i+j*shape[0]), arg0.dataset, arg0.stencil),
          NeighProxy<CellData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::SINGLE_MAT, i+j*shape[0]), arg1.dataset, arg1.stencil),
          NeighProxy<CellMatData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::SINGLE_MAT, i+j*shape[0]), arg2.dataset, arg2.stencil),
          NeighProxy<CellMatData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::SINGLE_MAT, i+j*shape[0]), arg3.dataset, arg3.stencil),
          data4[i+j*shape[0]]);
      } else {
        for (std::size_t structure_index = structure.structure[i+j*shape[0]].imat; 
            structure_index!=-1ul; structure_index = structure.mixed_storage[structure_index].nextfrac) {
          std::size_t mat_index = structure.mixed_storage[structure_index].material;
          anonymusAt3831(
            NeighProxy<CellData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::MULTIMAT, structure_index), arg0.dataset, arg0.stencil),
            NeighProxy<CellData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::MULTIMAT, structure_index), arg1.dataset, arg1.stencil),
            NeighProxy<CellMatData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::MULTIMAT, structure_index), arg2.dataset, arg2.stencil),
            NeighProxy<CellMatData<2>>(computation.data, Coords<2>(i, j), CellMatIndex(i+j*shape[0], mat_index), ValueIndex(ValueIndex::Type::MULTIMAT, structure_index), arg3.dataset, arg3.stencil),
            data4_list[structure_index]);
        }
      }
    }
  }
}
