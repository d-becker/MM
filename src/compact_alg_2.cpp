#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Computation.hpp"
#include <omp.h>

#define FUSED

using namespace MM;
using namespace MM::compressed_cell_centric;


void compact_alg_2(int sizex, int sizey,
    double *__restrict__ p_compact, double *__restrict__ p_compact_list,
    const double *__restrict__ n,
    const double *__restrict__ rho_compact, const double *__restrict__ rho_compact_list,
    const double *__restrict__ t_compact,
    const double *__restrict__ t_compact_list,
    const double *__restrict__ Vf,
    const double *__restrict__ Vf_compact_list,
          const CompressedDataStructure &structure) {

  double t1 = omp_get_wtime();

  #if defined(OMP)
  #pragma omp parallel for //collapse(2)
  #elif defined(ACC)
  #pragma acc parallel
  #pragma acc loop independent
  #endif
  for (int j = 0; j < sizey; j++) {
  #if defined(OMP)
  #pragma omp simd
  #elif defined(ACC)
  #pragma acc loop independent
  #endif
    for (int i = 0; i < sizex; i++) {

      if (structure.structure[i+sizex*j].nmats>1) {
#ifdef FUSED
        // NOTE: I think the paper describes this algorithm (Alg. 9) wrong.
        // The solution below is what I believe to good.

        // condition is 'ix >= 0', this is the equivalent of
        // 'until ix < 0' from the paper
#ifdef MM_LINKED
        for (std::size_t ix = structure.structure[i+sizex*j].imat;
             ix !=-1ul; ix = structure.mixed_storage[ix].nextfrac) {
          double nm = n[structure.mixed_storage[ix].material];
          p_compact_list[ix] = (nm * rho_compact_list[ix] * t_compact_list[ix]) / Vf_compact_list[ix];
        }
#else
        for (int idx = structure.structure[i+sizex*j].imat;
            idx < structure.structure[i+sizex*j].imat + structure.structure[i+sizex*j].nmats; idx++) {
          double nm = n[structure.mixed_storage[idx].material];
          p_compact_list[idx] = (nm * rho_compact_list[idx] * t_compact_list[idx]) / Vf_compact_list[idx];
        }
#endif
#endif
      }
      else {
        // NOTE: HACK: we index materials from zero, but zero can be a list index
        int mat = structure.structure[i+sizex*j].imat;
        // NOTE: There is no division by Vf here, because the fractional volume is 1.0 in the pure cell case.
        p_compact[i+sizex*j] = n[mat] * rho_compact[i+sizex*j] * t_compact[i+sizex*j];;
      }
    }
  }
#ifndef FUSED
  #if defined(OMP)
  #pragma omp parallel for simd
  #elif defined(ACC)
  #pragma acc parallel
  #pragma acc loop independent
  #endif
  for (int idx = 0; idx < mmc_index[mmc_cells]; idx++) {
    double nm = n[matids[idx]];
    p_compact_list[idx] = (nm * rho_compact_list[idx] * t_compact_list[idx]) / Vf_compact_list[idx];
  }
#endif

  printf("REFERENCE Compact matrix, cell centric, alg 2: %g sec\n", omp_get_wtime()-t1);
}
