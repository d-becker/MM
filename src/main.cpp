#include <iostream>

#include "Arguments.hpp"
#include "Data.hpp"
#include "Datasets.hpp"
#include "IndexGenerator.hpp"

using namespace std;
using namespace MM;

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
	MatData<> mat_data = data.new_mat_data();
	CellMatData<2> cell_mat_data = data.new_cell_mat_data();

        print_2D(dataset);
}

int main() {
	index_generator();
	return 0;
}
