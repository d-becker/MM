#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Computation.hpp"
#include <omp.h>

#define FUSED

using namespace MM;
using namespace MM::compressed_cell_centric;


void compact_alg_3(int sizex, int sizey,
    double *__restrict__ rho_mat_ave_compact, double *__restrict__ rho_mat_ave_compact_list,
    const double *__restrict__ x,
    const double *__restrict__ y,
    const double *__restrict__ rho_compact, const double *__restrict__ rho_compact_list,
          const CompressedDataStructure &structure) {

  double t1 = omp_get_wtime();
	// Computational loop 3 - Average density of each material over neighborhood of each cell
  t1 = omp_get_wtime();
  #if defined(OMP)
  #pragma omp parallel for //collapse(2)
  #elif defined(ACC)
  #pragma acc parallel
  #pragma acc loop independent
  #endif
	for (int j = 1; j < sizey-1; j++) {
  #if defined(OMP)
  #pragma omp simd
  #elif defined(ACC)
  #pragma acc loop independent
  #endif
		for (int i = 1; i < sizex-1; i++) {
			// o: outer
			double xo = x[i+sizex*j];
			double yo = y[i+sizex*j];

			// There are at most 9 neighbours in 2D case.
			double dsqr[9];

			// for all neighbours
			for (int nj = -1; nj <= 1; nj++) {
				for (int ni = -1; ni <= 1; ni++) {

					dsqr[(nj+1)*3 + (ni+1)] = 0.0;

					// i: inner
					double xi = x[(i+ni)+sizex*(j+nj)];
					double yi = y[(i+ni)+sizex*(j+nj)];

					dsqr[(nj+1)*3 + (ni+1)] += (xo - xi) * (xo - xi);
					dsqr[(nj+1)*3 + (ni+1)] += (yo - yi) * (yo - yi);
				}
			}

      if (structure.structure[i+sizex*j].nmats>1) {
				// condition is 'ix >= 0', this is the equivalent of
				// 'until ix < 0' from the paper
        #ifdef MM_LINKED
#pragma novector
        for (std::size_t ix = structure.structure[i+sizex*j].imat;
             ix !=-1ul; ix = structure.mixed_storage[ix].nextfrac) {
				#else
        for (int ix = structure.structure[i+sizex*j].imat;
            ix < structure.structure[i+sizex*j].imat + structure.structure[i+sizex*j].nmats; ix++) {
				#endif

					int mat = structure.mixed_storage[ix].material;
					double rho_sum = 0.0;
					int Nn = 0;

					// for all neighbours
					for (int nj = -1; nj <= 1; nj++) {
						for (int ni = -1; ni <= 1; ni++) {
							int ci = i+ni, cj = j+nj;

              if (structure.structure[ci+sizex*cj].nmats>1) {
								// condition is 'jx >= 0', this is the equivalent of
								// 'until jx < 0' from the paper
#ifdef MM_LINKED
#pragma novector
                for (std::size_t jx = structure.structure[ci+sizex*cj].imat;
                    jx !=-1ul; jx = structure.mixed_storage[jx].nextfrac) {
#else
                for (int jx = structure.structure[ci+sizex*cj].imat;
                      jx < structure.structure[ci+sizex*cj].imat + structure.structure[ci+sizex*cj].nmats; jx++) {
#endif
									if (structure.mixed_storage[jx].material == mat) {
										rho_sum += rho_compact_list[jx] / dsqr[(nj+1)*3 + (ni+1)];
										Nn += 1;

										// The loop has an extra condition: "and not found".
										// This makes sense, if the material is found, there won't be any more of the same.
										break;
									}
								}
							}
							else {
								// NOTE: In this case, the neighbour is a pure cell, its material index is in jx.
								// In contrast, Algorithm 10 loads matids[jx] which I think is wrong.

								// NOTE: HACK: we index materials from zero, but zero can be a list index
								int mat_neighbour = structure.structure[ci+sizex*cj].imat;
								if (mat == mat_neighbour) {
									rho_sum += rho_compact[ci+sizex*cj] / dsqr[(nj+1)*3 + (ni+1)];
									Nn += 1;
								}
							} // end if (jx <= 0)
						} // end for (int ni)
					} // end for (int nj)

					rho_mat_ave_compact_list[ix] = rho_sum / Nn;
				} // end for (ix = -ix)
			} // end if (ix <= 0)
			else {
				// NOTE: In this case, the cell is a pure cell, its material index is in ix.
				// In contrast, Algorithm 10 loads matids[ix] which I think is wrong.

				// NOTE: HACK: we index materials from zero, but zero can be a list index
				int mat = structure.structure[i+sizex*j].imat;

				double rho_sum = 0.0;
				int Nn = 0;

				// for all neighbours
				for (int nj = -1; nj <= 1; nj++) {
					if ((j + nj < 0) || (j + nj >= sizey)) // TODO: better way?
						continue;

					for (int ni = -1; ni <= 1; ni++) {
						if ((i + ni < 0) || (i + ni >= sizex)) // TODO: better way?
							continue;

						int ci = i+ni, cj = j+nj;

            if (structure.structure[ci+sizex*cj].nmats>1) {
							// condition is 'jx >= 0', this is the equivalent of
							// 'until jx < 0' from the paper
#ifdef MM_LINKED
#pragma novector
                for (std::size_t jx = structure.structure[ci+sizex*cj].imat;
                    jx !=-1ul; jx = structure.mixed_storage[jx].nextfrac) {
#else
                for (int jx = structure.structure[ci+sizex*cj].imat;
                      jx < structure.structure[ci+sizex*cj].imat + structure.structure[ci+sizex*cj].nmats; jx++) {
#endif

                if (structure.mixed_storage[jx].material == mat) {
									rho_sum += rho_compact_list[jx] / dsqr[(nj+1)*3 + (ni+1)];
									Nn += 1;

									// The loop has an extra condition: "and not found".
									// This makes sense, if the material is found, there won't be any more of the same.
									break;
								}
							}
						}
						else {
							// NOTE: In this case, the neighbour is a pure cell, its material index is in jx.
							// In contrast, Algorithm 10 loads matids[jx] which I think is wrong.

							// NOTE: HACK: we index materials from zero, but zero can be a list index
              int mat_neighbour = structure.structure[ci+sizex*cj].imat;
							if (mat == mat_neighbour) {
								rho_sum += rho_compact[ci+sizex*cj] / dsqr[(nj+1)*3 + (ni+1)];
								Nn += 1;
							}
						} // end if (jx <= 0)
					} // end for (int ni)
				} // end for (int nj)

				rho_mat_ave_compact[i+sizex*j] = rho_sum / Nn;
			} // end else
		}
	}
  printf("REFERENCE Compact matrix, cell centric, alg 3: %g sec\n", omp_get_wtime()-t1);
}
