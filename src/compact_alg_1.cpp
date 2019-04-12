#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Computation.hpp"
#include <omp.h>

#define FUSED

using namespace MM;
using namespace MM::compressed_cell_centric;


void compact_alg_1(int sizex, int sizey,
    const double *__restrict__ rho_compact, const double *__restrict__ rho_compact_list,
    const double *__restrict__ V, 
    const double *__restrict__ Vf, 
    const double *__restrict__ Vf_compact_list, 
          double *__restrict__ rho_ave_compact, 
          CompressedDataStructure &structure) {
// Computational loop 1 - average density in cell
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

#ifdef FUSED
      double ave = 0.0;
      if (structure.structure[i+sizex*j].nmats>1) {
#ifdef MM_LINKED
#pragma novector
        for (std::size_t ix = structure.structure[i+sizex*j].imat;
             ix !=-1ul; ix = structure.mixed_storage[ix].nextfrac) {
          ave += rho_compact_list[ix] * Vf_compact_list[ix];
        }
#else
        for (int idx = structure.structure[i+sizex*j].imat;
            idx < structure.structure[i+sizex*j].imat + structure.structure[i+sizex*j].nmats; idx++) {
          ave += rho_compact_list[idx] * Vf_compact_list[idx];
        }
#endif
        rho_ave_compact[i+sizex*j] = ave/V[i+sizex*j];
      }
      else {
#endif
        // We use a distinct output array for averages.
        // In case of a pure cell, the average density equals to the total.
        rho_ave_compact[i+sizex*j] = rho_compact[i+sizex*j] / V[i+sizex*j];
#ifdef FUSED
      }
#endif
    }
  }
#ifndef FUSED
  #if defined(OMP)
  #pragma omp parallel for simd
  #elif defined(ACC)
  #pragma acc parallel
  #pragma acc loop independent
  #endif
  for (int c = 0; c < mmc_cells; c++) {
    double ave = 0.0;
    for (int m = mmc_index[c]; m < mmc_index[c+1]; m++) {
      ave +=  rho_compact_list[m] * Vf_compact_list[m];
    }
    rho_ave_compact[mmc_i[c]+sizex*mmc_j[c]] = ave/V[mmc_i[c]+sizex*mmc_j[c]];
  }
#endif
  printf("REFERENCE Compact matrix, cell centric, alg 1: %g sec\n", omp_get_wtime()-t1);
}
