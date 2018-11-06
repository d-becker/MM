/*
* Open source copyright declaration based on BSD open source template:
* http://www.opensource.org/licenses/bsd-license.php
*
* Copyright (c) 2013, Istvan Reguly and others. 
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
* The name of Mike Giles may not be used to endorse or promote products
* derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY Istvan Reguly ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/** @brief initial version of mutli-material code with full dense matrix representaiton
  * @author Istvan Reguly
  */


#include <cmath>
#include <iostream>

#include "Arguments.hpp"
#include "Computation.hpp"
#include "Data.hpp"
#include "Datasets.hpp"
#include "IndexGenerator.hpp"

#include "CompressedDataStructure.hpp"

using namespace std;
using namespace MM;
using namespace MM::full_matrix;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

extern void full_matrix_cell_centric(unsigned int sizex, unsigned int sizey, int Nmats, Data<2> &d,
	CellMatData<2>& rho, CellMatData<2>& rho_mat_ave, CellMatData<2>& p, CellMatData<2>& Vf, CellMatData<2>& t,
	CellData<2>& V, CellData<2>& x, CellData<2>& y,
	MatData<>& n, CellData<2>& rho_ave);

extern void compact_matrix_cell_centric(unsigned int sizex, unsigned int sizey, int Nmats, Data<2> &d,
	CellMatData<2>& rho, CellMatData<2>& rho_mat_ave, CellMatData<2>& p, CellMatData<2>& Vf, CellMatData<2>& t,
	CellData<2>& V, CellData<2>& x, CellData<2>& y,
	MatData<>& n, CellData<2>& rho_ave, vector<vector<size_t>> &mats);

void initialise_field_rand(CellMatData<2>& rho, CellMatData<2>& t, CellMatData<2>& p, int Nmats, unsigned int sizex, unsigned int sizey, double prob2, double prob3, double prob4) {

  //let's use a morton space filling curve here
  srand(0);
  double prob1 = 1.0-prob2-prob3-prob4;
  printf("Random layout %g %g %g %g\n", prob1, prob2, prob3, prob4);

  for (int n = 0; n < sizex*sizey; n++) {
    unsigned int i = n%sizex;//n & 0xAAAA;
    unsigned int j = n/sizex;//n & 0x5555;

    double r = (double)rand()/(double)RAND_MAX;
    int m = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
    int m2, m3, m4;
    rho.at(Coords<2>(i,j),m) = 1.0;
    t.at(Coords<2>(i,j),m) = 1.0;
    p.at(Coords<2>(i,j),m) = 1.0;
    if (r >= prob1) {
      m2 = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
      while (m2 == m)
        m2 = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
      rho.at(Coords<2>(i,j), m2) = 1.0;
      t.at(Coords<2>(i,j), m2) = 1.0;
      p.at(Coords<2>(i,j), m2) = 1.0;
    }
    if (r >= 1.0-prob4-prob3) {
      m3 = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
      while (m3 == m && m3 == m2)
        m3 = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
      rho.at(Coords<2>(i,j),m3) = 1.0;
      t.at(Coords<2>(i,j),m3) = 1.0;
      p.at(Coords<2>(i,j),m3) = 1.0;
    }
    if (r >= 1.0-prob4) {
      m4 = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
      while (m4 == m && m4 == m2 && m4 == m3)
        m4 = (double)rand()/(double)RAND_MAX * Nmats/4 + (Nmats/4)*(i/(sizex*sizey/4));
      rho.at(Coords<2>(i,j),m4) = 1.0;
      t.at(Coords<2>(i,j),m4) = 1.0;
      p.at(Coords<2>(i,j),m4) = 1.0;
    }
  }
}

  void initialise_field_static(CellMatData<2>& rho, CellMatData<2>& t, CellMatData<2>& p, int Nmats, unsigned int sizex, unsigned int sizey) {
	//Pure cells and simple overlaps
	int width = sizex/Nmats;

  int overlap_i = std::max(0.0,ceil((double)sizey/1000.0)-1);
  int overlap_j = std::max(0.0,floor((double)sizex/1000.0)-1);

	//Top
	for (unsigned int mat = 0; mat < Nmats/2; mat++) {
#pragma omp parallel for
		for (unsigned int j = mat*width; j < sizey/2; j++) {
			for (unsigned int i = mat*width-(mat>0)-(mat>0)*overlap_i; i < (mat+1)*width; i++) { //+1 for overlap
				rho.at(Coords<2>(i,j),mat) = 1.0;
				t.at(Coords<2>(i,j),mat) = 1.0;
				p.at(Coords<2>(i,j),mat) = 1.0;
			}
			for (unsigned int i = sizex-mat*width-1+(mat>0)*overlap_i; i >= sizex-(mat+1)*width-1; i--) { //+1 for overlap
				rho.at(Coords<2>(i,j),mat) = 1.0;
				t.at(Coords<2>(i,j),mat) = 1.0;
				p.at(Coords<2>(i,j),mat) = 1.0;
			}
		}

#pragma omp parallel for
		for (unsigned int j = mat*width-(mat>0); j < (mat+1)*width; j++) { //+1 for overlap
			for (unsigned int i = mat*width-(mat>0)-(mat>0)*overlap_i; i < sizex-mat*width; i++) {
				rho.at(Coords<2>(i,j),mat) = 1.0;
				t.at(Coords<2>(i,j),mat) = 1.0;
				p.at(Coords<2>(i,j),mat) = 1.0;
			}
		}
	}
	
	//Bottom
	for (unsigned int mat = 0; mat < Nmats/2; mat++) {
#pragma omp parallel for
		for (unsigned int j = sizey/2-1; j < sizey-mat*width; j++) {
			for (unsigned int i = mat*width-(mat>0); i < (mat+1)*width; i++) { //+1 for overlap
				rho.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
				t.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
				p.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
			}
			for (unsigned int i = sizex-mat*width-1; i >= sizex-(mat+1)*width-1; i--) { //+1 for overlap
				rho.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
				t.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
				p.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
			}
		}
#pragma omp parallel for
		for (unsigned int j = sizey-mat*width-1+(mat>0)*overlap_j; j >= sizey-(mat+1)*width-(mat<(Nmats/2-1)); j--) { //+1 for overlap
			for (unsigned int i = mat*width; i < sizex-mat*width; i++) {
				rho.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
				t.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
				p.at(Coords<2>(i,j),mat+Nmats/2) = 1.0;
			}
		}
	}
	//Fill in corners
#pragma omp parallel for
	for (unsigned int mat = 1; mat < Nmats/2; mat++) {
		for (unsigned int j = sizey/2-3; j < sizey/2-1;j++)
			for (unsigned int i = 2; i < 5+overlap_i; i++) {
				//x neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-1,j),mat) = 1.0;t.at(Coords<2>(mat*width-1,j),mat) = 1.0;p.at(Coords<2>(mat*width-1,j),mat) = 1.0;
				//y neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-i,j),Nmats/2+mat) = 1.0;t.at(Coords<2>(mat*width-i,j),Nmats/2+mat) = 1.0;p.at(Coords<2>(mat*width-i,j),Nmats/2+mat) = 1.0;
				//x-y neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat) = 1.0;t.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat) = 1.0;p.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat) = 1.0;
				rho.at(Coords<2>(mat*width-i,j),Nmats/2+mat-1) = 1.0;t.at(Coords<2>(mat*width-i,j),Nmats/2+mat-1) = 1.0;p.at(Coords<2>(mat*width-i,j),Nmats/2+mat-1) = 1.0;
			}
		for (unsigned int j = sizey/2; j < sizey/2+2+overlap_j;j++)
			for (unsigned int i = 2; i < 5; i++) {
				//x neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),Nmats/2+mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-i,j),Nmats/2+mat) = 1.0;t.at(Coords<2>(mat*width-i,j),Nmats/2+mat) = 1.0;p.at(Coords<2>(mat*width-i,j),Nmats/2+mat) = 1.0;
				//y neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-1,j),mat) = 1.0;t.at(Coords<2>(mat*width-1,j),mat) = 1.0;p.at(Coords<2>(mat*width-1,j),mat) = 1.0;

			}
	}
	int only_8 = 0;
	for (unsigned int mat = Nmats/2+1; mat < Nmats; mat++) {
		for (unsigned int j = sizey/2-3; j < sizey/2-1;j++)
			for (unsigned int i = 2; i < 5; i++) {
				//x neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-i,j),-Nmats/2+mat) = 1.0;t.at(Coords<2>(mat*width-i,j),-Nmats/2+mat) = 1.0;p.at(Coords<2>(mat*width-i,j),-Nmats/2+mat) = 1.0;
				//y neighbour material
				rho.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-1,j),mat) = 1.0;t.at(Coords<2>(mat*width-1,j),mat) = 1.0;p.at(Coords<2>(mat*width-1,j),mat) = 1.0;
			}
		for (unsigned int j = sizey/2; j < sizey/2+2;j++)
			for (unsigned int i = 2; i < 4; i++) {
				if (i < 3 && only_8<6) {
					//y neighbour material
					rho.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat-1) = 1.0;
					rho.at(Coords<2>(mat*width-i,j),-Nmats/2+mat) = 1.0;t.at(Coords<2>(mat*width-i,j),-Nmats/2+mat) = 1.0;p.at(Coords<2>(mat*width-i,j),-Nmats/2+mat) = 1.0;
				}
				if (i==2 && only_8==0) {
					//x-y neighbour material
					rho.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat) = 1.0;t.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat) = 1.0;p.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat) = 1.0;
					rho.at(Coords<2>(mat*width-i,j),-Nmats/2+mat-1) = 1.0;t.at(Coords<2>(mat*width-i,j),-Nmats/2+mat-1) = 1.0;p.at(Coords<2>(mat*width-i,j),-Nmats/2+mat-1) = 1.0;
				}
				//x neighbour material
				if (mat >= Nmats-8 && j==sizey/2+1 && i==3) if (only_8++>=4) {
					break;
				}
				rho.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;t.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;p.at(Coords<2>(mat*width+i-2,j),mat-1) = 1.0;
				rho.at(Coords<2>(mat*width-1,j),mat) = 1.0;t.at(Coords<2>(mat*width-1,j),mat) = 1.0;p.at(Coords<2>(mat*width-1,j),mat) = 1.0;
			}
	}
#pragma omp parallel for
	for (unsigned int mat=Nmats/2+1; mat < Nmats/2+5; mat++) {
		unsigned int i = 2; unsigned int j = sizey/2+1;
		rho.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat) = 0.0;t.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat) = 0.0;p.at(Coords<2>(mat*width+i-2,j),-Nmats/2+mat) = 0.0;
	}
}
int main(int argc, char* argv[]) {
	unsigned int sizex = 1000;
  if (argc > 1)
    sizex = atoi(argv[1]);
	unsigned int sizey = 1000;
  if (argc > 2)
    sizey = atoi(argv[2]);
	unsigned int ncells = sizex*sizey;

	int Nmats = 40;

	const std::size_t COLS = sizex;
    const std::size_t ROWS = sizey;
    const std::size_t MAT_N = Nmats;

    Data<2> data({COLS, ROWS}, MAT_N);


	//Allocate the four state variables for all Nmats materials and all cells 
	//density
    CellMatData<2> rho = data.new_cell_mat_data();
    //density average in neighbourhood
    CellMatData<2> rho_mat_ave = data.new_cell_mat_data();
    //pressure
    CellMatData<2> p = data.new_cell_mat_data();
    //Fractional volume
    CellMatData<2> Vf = data.new_cell_mat_data();
    //temperature
    CellMatData<2> t = data.new_cell_mat_data();

	//Allocate per-cell only datasets
	CellData<2> V = data.new_cell_data();
	CellData<2> x = data.new_cell_data();
	CellData<2> y = data.new_cell_data();

	//Allocate per-material only datasets
	MatData<> n = data.new_mat_data(); // number of mats

	//Allocate output datasets
	CellData<2> rho_ave = data.new_cell_data();

	int imaterial_multi_cell;

	//Initialise arrays
	double dx = 1.0/sizex;
	double dy = 1.0/sizey;
	for (unsigned int j = 0; j < sizey; j++) {
		for (unsigned int i = 0; i < sizex; i++) {
			V.at(Coords<2>(i,j)) = dx*dy;
			x.at(Coords<2>(i,j)) = dx*i;
			y.at(Coords<2>(i,j)) = dy*j;
		}
	}

	for (unsigned int mat = 0; mat < Nmats; mat++) {
		n.at(mat) = 1.0; // dummy value
	}

  if (argc>=6) initialise_field_rand(rho, t, p, Nmats, sizex, sizey, atof(argv[3]), atof(argv[4]), atof(argv[5]));
  else initialise_field_static(rho, t, p, Nmats, sizex, sizey);

	FILE *f;
	int print_to_file = 0;

	if (print_to_file==1)
		FILE *f = fopen("map.txt","w");

  std::vector<std::vector<std::size_t>> mats(sizex*sizey);
	//Compute fractions and count cells
	int cell_counts_by_mat[4] = {0,0,0,0};
  int mmc_cells = 0;
	for (unsigned int j = 0; j < sizey; j++) {
		for (unsigned int i = 0; i < sizex; i++) {
			int count = 0;
			for (unsigned int mat = 0; mat < Nmats; mat++) {
				count += rho.at(Coords<2>(i,j),mat)!=0.0;
			}
			if (count == 0) {
				printf("Error: no materials in cell %d %d\n",i,j);
				if (print_to_file)
					fclose(f);
        exit(-1);
			}
      if (count > 1) mmc_cells++;
			cell_counts_by_mat[count-1]++;
      mats[j*sizex+i].resize(count);
      count = 0;
			for (unsigned int mat = 0; mat < Nmats; mat++) {
        if (rho.at(Coords<2>(i,j),mat)!=0.0) {
          mats[j*sizex+i][count++] = mat;
        }
			}


			if (print_to_file) {
				if (i!=0) fprintf(f,", %d",count);
				else fprintf(f,"%d",count);
			}

			for (unsigned int mat = 0; mat < Nmats; mat++) {
				if (rho.at(Coords<2>(i,j),mat)!=0.0) Vf.at(Coords<2>(i,j),mat)=1.0/count;
			}
		}
		if (print_to_file)
			fprintf(f,"\n");
	}
	printf("Pure cells %d, 2-materials %d, 3 materials %d, 4 materials %d: MMC cells %d\n",
		cell_counts_by_mat[0],cell_counts_by_mat[1],cell_counts_by_mat[2],cell_counts_by_mat[3], mmc_cells);

	if (print_to_file)
		fclose(f);


	full_matrix_cell_centric(sizex, sizey, Nmats, data, rho, rho_mat_ave, p, Vf, t, V, x, y, n, rho_ave);
  compact_matrix_cell_centric(sizex, sizey, Nmats, data, rho, rho_mat_ave, p, Vf, t, V, x, y, n, rho_ave, mats);
	
	return 0;
}
