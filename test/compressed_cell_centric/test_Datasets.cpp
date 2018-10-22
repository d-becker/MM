#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "gtest/gtest.h"

#include "Coords.hpp"
#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "compressed_cell_centric/Datasets.hpp"

namespace MM::compressed_cell_centric {

namespace {

// TODO: Test accessor functions.

// class CompressedCellMatData : public ::testing::Test {
// public:
// 	CompressedCellMatData()
// 		: cols(2),
// 		  rows(2),
// 		  size({cols, rows}),
// 		  mat_no(4),
// 		  materials({
// 				  std::vector<std::size_t>{0},
// 				  std::vector<std::size_t>{0, 1},
// 				  std::vector<std::size_t>{2, 3},
// 				  std::vector<std::size_t>{1, 3},
// 			    }),
// 		  structure(materials),
// 		  cell_values({2.8, 0.0, 0.0, 0.0}),
// 		  mixed_storage_values({2.6, 8.2, 5.82, 6.38, 1.24, 3.4})
// 	{
// 	}
// protected:
// 	const std::size_t cols;
// 	const std::size_t rows;
// 	const std::array<std::size_t, 2> size;
// 	const std::size_t mat_no;
// 	const std::vector<std::vector<std::size_t>> materials;
// 	CompressedDataStructure structure;
// 	std::vector<double> cell_values;
// 	std::vector<double> mixed_storage_values;
// };

// TEST_F(CompressedCellMatData, get_size_ok) {
// 	CellMatData<2> cell_mat_data(size,
// 				     structure,
// 				     cell_values,
// 				     mixed_storage_values);
// 	ASSERT_EQ(cell_mat_data.get_size(), size);
// }

// TEST_F(CompressedCellMatData, access_single_mat_cell) {
// 	CellMatData<2> cell_mat_data(size,
// 				     structure,
// 				     cell_values,
// 				     mixed_storage_values);
// 	const Coords<2> index = {0u, 0u};
// 	const double value = cell_mat_data.at(index, 0);
// 	ASSERT_EQ(value, cell_values.at(0));

// 	const double new_value = value + 1.0;
// 	cell_mat_data.at(index, 0) = new_value;
// 	ASSERT_EQ(new_value, cell_mat_data.at(index, 0));
// }

// TEST_F(CompressedCellMatData, access_multimat_cell) {
// 	CellMatData<2> cell_mat_data(size,
// 				     structure,
// 				     cell_values,
// 				     mixed_storage_values);
// 	const Coords<2> index = {1u, 0u};
// 	const std::size_t mat_index = 0;
// 	const double value = cell_mat_data.at(index, mat_index);
// 	ASSERT_EQ(value, mixed_storage_values.at(0));

// 	const double new_value = value + 1.0;
// 	cell_mat_data.at(index, mat_index) = new_value;
// 	ASSERT_EQ(new_value, cell_mat_data.at(index, mat_index));
// }

// TEST_F(CompressedCellMatData, single_mat_iteration) {
// 	CellMatData<2> cell_mat_data(size,
// 				     structure,
// 				     cell_values,
// 				     mixed_storage_values);
// 	const Coords<2> index = {0u, 0u};

// 	const std::vector<double> expected = {cell_values.at(0)};

// 	auto iteration_proxy = cell_mat_data.mat_iteration(index);
// 	std::vector<double> actual;
// 	for (double value : iteration_proxy) {
// 		actual.push_back(value);
// 	}

// 	ASSERT_EQ(actual, expected);
// }

// TEST_F(CompressedCellMatData, multimat_iteration) {
// 	CellMatData<2> cell_mat_data(size,
// 				     structure,
// 				     cell_values,
// 				     mixed_storage_values);
// 	const Coords<2> index = {1u, 0u};

// 	const std::vector<double> expected = {mixed_storage_values.at(0),
// 					      mixed_storage_values.at(1)};

// 	auto iteration_proxy = cell_mat_data.mat_iteration(index);
// 	std::vector<double> actual;
// 	for (double value : iteration_proxy) {
// 		actual.push_back(value);
// 	}

// 	ASSERT_EQ(actual, expected);

// 	std::vector<double> new_expected;
// 	std::transform(expected.begin(),
// 		       expected.end(),
// 		       std::back_inserter(new_expected),
// 		       [](const double value) {return value + 1.0;});

// 	for (double& value : iteration_proxy) {
// 		value += 1.0;
// 	}

// 	std::vector<double> new_actual;
// 	for (double value : iteration_proxy) {
// 		new_actual.push_back(value);
// 	}

// 	ASSERT_EQ(new_actual, new_expected);
// }

} // anonymous namespace

} // MM::compressed_cell_centric
