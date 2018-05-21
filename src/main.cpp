#include <iostream>

#include "Data.hpp"
#include "Datasets.hpp"

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

int main() {
	Data<2> data({2, 2}, 4);
	CellData<2> dataset = data.new_cell_data();
	MatData<> mat_data = data.new_mat_data();

        print_2D(dataset);
	return 0;
}
