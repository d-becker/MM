#include <cmath>
#include <iostream>

#include "IndexGenerator.hpp"
#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "full_matrix/Arguments.hpp"
#include "full_matrix/Computation.hpp"
#include "full_matrix/Data.hpp"
#include "full_matrix/Datasets.hpp"

using namespace std;
using namespace MM;
using namespace MM::full_matrix;

void print_2D(const CellData<2, double>& dataset) {
	for (std::size_t i = 0; i < dataset.get_size().at(0); ++i) {
		for (std::size_t j = 0; j < dataset.get_size().at(1); ++j) {
			cout << dataset[Coords<2>(i, j)] << "\t";
		}

		cout << endl;
	}
}

template<std::size_t N>
void print_array(const std::array<std::size_t, N> arr) {
	for (const std::size_t e : arr) {
		cout << e << "\t";
	}

	cout << endl;
}

void index_generator() {
	std::array<std::size_t, 3> begin {0, 1, 0};
	std::array<std::size_t, 3> end {2, 4, 2};

	IndexGenerator<3> generator(begin, end);

	while (generator.has_next()) {
		std::array<std::size_t, 3> array = generator.next();
		print_array(array);
	}
}

void datasets() {
	Data<2> data({2, 2}, 4);
	CellData<2> dataset = data.new_cell_data();
	// MatData<> mat_data = data.new_mat_data();
	// CellMatData<2> cell_mat_data = data.new_cell_mat_data();

        print_2D(dataset);
}

void fill_density(CellMatData<2>& density) {
	// density(:,:,0) =
	// 0.87093   0.57531
	// 0.69274   0.84571
	//
	// density(:,:,1) =
	// 0.98574   0.93665
	// 0.10462   0.96086

	density.at(Coords<2>(0u, 0u), 0) = 0.87093;
	density.at(Coords<2>(0u, 1u), 0) = 0.57531;
	density.at(Coords<2>(1u, 0u), 0) = 0.69274;
	density.at(Coords<2>(1u, 1u), 0) = 0.845713;

	density.at(Coords<2>(0u, 0u), 1) = 0.98574;
	density.at(Coords<2>(0u, 1u), 1) = 0.93665;
	density.at(Coords<2>(1u, 0u), 1) = 0.10462;
	density.at(Coords<2>(1u, 1u), 1) = 0.96086;
}

void fill_volume(CellMatData<2>& volume) {
	// volume(:,:,0) =
	// 0.052457   0.616349
	// 0.103264   0.413648
	//
	// volume(:,:,1) =
	// 0.51597   0.79952
	// 0.19611   0.10866

	volume.at(Coords<2>(0u, 0u), 0) = 0.052457;
	volume.at(Coords<2>(0u, 1u), 0) = 0.616349;
	volume.at(Coords<2>(1u, 0u), 0) = 0.103264;
	volume.at(Coords<2>(1u, 1u), 0) = 0.413648;

	volume.at(Coords<2>(0u, 0u), 1) = 0.51597;
	volume.at(Coords<2>(0u, 1u), 1) = 0.79952;
	volume.at(Coords<2>(1u, 0u), 1) = 0.19611;
	volume.at(Coords<2>(1u, 1u), 1) = 0.10866;
}

bool eq(double lhs, double rhs, double epsilon = 0.0001) {
	return abs(lhs - rhs) <= epsilon;
}

bool check_mass(const CellMatData<2>& mass) {
	// mass(:,:,0) =
	// 0.045686   0.354592
	// 0.071535   0.349827
	//
	// mass(:,:,1) =
	// 0.508612   0.748873
	// 0.020517   0.104403

	return
		eq(mass.at(Coords<2>(0u, 0u), 0), 0.045686) &&
		eq(mass.at(Coords<2>(0u, 1u), 0), 0.354592) &&
		eq(mass.at(Coords<2>(1u, 0u), 0), 0.071535) &&
		eq(mass.at(Coords<2>(1u, 1u), 0), 0.349827) &&

		eq(mass.at(Coords<2>(0u, 0u), 1), 0.508612) &&
		eq(mass.at(Coords<2>(0u, 1u), 1), 0.748873) &&
		eq(mass.at(Coords<2>(1u, 0u), 1), 0.020517) &&
		eq(mass.at(Coords<2>(1u, 1u), 1), 0.104403);
}

/* void mm_par_loop_kernel1(std::string name,
     Computation<2, double> &computation,
     IN<CellMatData<2>> arg0,
     IN<CellMatData<2>> arg1,
     OUT<CellMatData<2>> arg2);*/

void test1() {
	// Given a density and a volume state variable, calculate
	// the mass by cell and material (no reduction).

	const std::size_t COLS = 2;
	const std::size_t ROWS = 2;
	const std::size_t MAT_N = 2;
	
	Data<2> data({COLS, ROWS}, MAT_N);
	
	CellMatData<2> density = data.new_cell_mat_data();
	CellMatData<2> volume = data.new_cell_mat_data();

	fill_density(density);
	fill_volume(volume);

	// Output
	CellMatData<2> mass = data.new_cell_mat_data();

	IndexGenerator<2> index_generator({0, 0}, {2, 2});
	Computation<2> computation(data, index_generator);
	
	computation.compute([] (double density,
				double volume,
				double& mass) {
		mass = density * volume;
	},
			    IN<CellMatData<2>>(density),
			    IN<CellMatData<2>>(volume),
			    OUT<CellMatData<2>>(mass));
  /*m_par_loop_kernel1("kernel1", computation,
           IN<CellMatData<2>>(density),
           IN<CellMatData<2>>(volume),
           OUT<CellMatData<2>>(mass));*/

	bool ok = check_mass(mass);
	cout << "Mass ok: " << ok << "." << endl;
}

bool check_mass_by_cell(const CellData<2>& mass) {
	// mass_by_cell =
	// 0.554299   1.103465
	// 0.092052   0.454230

	return
		eq(mass[Coords<2>(0u, 0u)], 0.554299) &&
		eq(mass[Coords<2>(0u, 1u)], 1.103465) &&
		eq(mass[Coords<2>(1u, 0u)], 0.092052) &&
		eq(mass[Coords<2>(1u, 1u)], 0.454230);
}


void test2() {
	// Given a density and a volume state variable, calculate
	// the mass by cell (with reduction).

	const std::size_t COLS = 2;
	const std::size_t ROWS = 2;
	const std::size_t MAT_N = 2;

	Data<2> data({COLS, ROWS}, MAT_N);
	
	CellMatData<2> density = data.new_cell_mat_data();
	CellMatData<2> volume = data.new_cell_mat_data();

	fill_density(density);
	fill_volume(volume);

	// Output
	CellData<2> mass_by_cell = data.new_cell_data();

	auto INC = [] (double left, double right) {
		return left + right;
	};
	
	IndexGenerator<2> index_generator({0, 0}, {2, 2});
	Computation<2> computation(data, index_generator);
	computation.compute([] (double density,
				double volume,
				ReduceProxy<double> mass_by_cell) {
		mass_by_cell << density * volume;
	},
			    IN<CellMatData<2>>(density),
			    IN<CellMatData<2>>(volume),
			    REDUCE<CellData<2>>(INC, mass_by_cell));

	bool ok = check_mass_by_cell(mass_by_cell);
	cout << "Mass by cell ok: " << ok << "." << endl;
}

void test3() {
	const std::size_t COLS = 128;
	const std::size_t ROWS = 128;
	const std::size_t MAT_N = 1;

	Data<2> data({COLS, ROWS}, MAT_N);

	CellData<2> x = data.new_cell_data();
	CellData<2> y = data.new_cell_data();

	Stencil<2> s9pt({{1,1},  {1,0},  {1,-1},
			 {0,1},  {0,0},  {0,-1},
			 {-1,1}, {-1,0}, {-1,-1}});

        // Fill the datasets with data.
	
	IndexGenerator<2> index_generator({1, 1}, {127, 127});
	Computation<2> computation(data, index_generator);
	
	computation.compute([] (NeighProxy<CellData<2>> x,
				double& y) {
				    y = -x[{1,1}]  - x[{1,0}]   - x[{1,-1}]
					-x[{0,1}]  + 8*x[{0,0}] - x[{0,-1}]
					-x[{-1,1}] - x[{-1,0}]  -x[{-1,-1}];
			    },
			    NEIGH<CellData<2>>(x, s9pt),
			    OUT<CellData<2>>(y));
}

void test4_index() {
	const std::size_t COLS = 2;
	const std::size_t ROWS = 2;
	const std::size_t MAT_N = 1;

	Data<2> data({COLS, ROWS}, MAT_N);

	CellData<2> y = data.new_cell_data();

        // Fill the datasets with data.
	
	IndexGenerator<2> index_generator({0, 0}, {2, 2});
	Computation<2> computation(data, index_generator);
	
	computation.compute([] (Coords<2> index,
				double& y) {
      std::cout << index[0] << ", " << index[1] << std::endl;
			    },
			    INDEX<2>(),
			    OUT<CellData<2>>(y));
}

void test5_free_data() {
	const std::size_t COLS = 128;
	const std::size_t MAT_N = 50;

	Data<1> data({COLS}, MAT_N);

	const double dt = 1e-2;

	CellData<1> input = data.new_cell_data();
	CellData<1> output = data.new_cell_data();
	
	Stencil<1> neighs({{-1}, {0}, {1}});

        // Fill the dataset with data.
	
	IndexGenerator<1> index_generator({1}, {COLS - 1});
	Computation<1> computation(data, index_generator);
	
	computation.compute([] (
		const double& dt,
		const NeighProxy<CellData<1>> input,
		double& output) {
		    output = dt * (-input[{-1}] + 2*input[{0}] - input[{1}]);
	        },
	        FREE_SCALAR<>(dt),
                NEIGH<CellData<1>>(input, neighs),
	        OUT<CellData<1>>(output));
}

int main() {
	// index_generator();
	test1();
	test2();
	test3();
	test4_index();
	test5_free_data();

	return 0;
}
