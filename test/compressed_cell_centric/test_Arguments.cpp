#include <array>

#include "gtest/gtest.h"

#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Datasets.hpp"

namespace MM::compressed_cell_centric {

namespace {

TEST(test_IN, cell_data) {
	MultidimArray<2, double> arr({2, 2});
	const double value = 2.4;
	const std::size_t cell_index = 1;
	arr.at_raw_index(cell_index) = value;
	
	const CellData<2> cell_data(arr);
	const IN<CellData<2>> in(cell_data);

	const CellMatIndex cell_mat_index(cell_index, 0);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(in.get(cell_mat_index, value_index), value);
}

TEST(test_IN, mat_data) {
	std::vector<double> arr(4, 0.0);
	const double value = 2.4;
	arr.at(1) = value;
	
	const MatData<> mat_data(arr);
	const IN<MatData<>> in(mat_data);

	const CellMatIndex cell_mat_index(0, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(in.get(cell_mat_index, value_index), value);
}

TEST(test_IN, cell_mat_data) {
	std::vector<double> cell_values{2.2, 1.4, 3.1, 4.6};
	std::vector<double> mixed_values{5.6, 8.2, 6.8, 4.3};
	
	const CellMatData<2> cell_mat_data(cell_values, mixed_values);
	const IN<CellMatData<2>> in(cell_mat_data);

	const CellMatIndex cell_mat_index(1, 1);
	const std::size_t mixed_index = 2;
	const ValueIndex value_index(ValueIndex::Type::MULTIMAT, mixed_index);

	ASSERT_EQ(in.get(cell_mat_index, value_index),
		  mixed_values.at(mixed_index));
}

TEST(test_OUT, cell_data) {
	MultidimArray<2, double> arr({2, 2});
	
	const CellData<2> cell_data(arr);
	OUT<CellData<2>> out(cell_data);

	const double value = 2.4;
        const std::size_t cell_index = 1;
	const CellMatIndex cell_mat_index(cell_index, 0);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	out.get(cell_mat_index, value_index) = value;
	ASSERT_EQ(arr.at_raw_index(cell_index), value);
}

TEST(test_OUT, mat_data) {
	std::vector<double> arr(4, 0.0);
	
	const MatData<> mat_data(arr);
	OUT<MatData<>> out(mat_data);

	const double value = 2.4;
        const std::size_t mat_index = 1;
	const CellMatIndex cell_mat_index(0, mat_index);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	out.get(cell_mat_index, value_index) = value;
	ASSERT_EQ(arr.at(mat_index), value);
}

TEST(test_OUT, cell_mat_data) {
	std::vector<double> cell_values{2.2, 1.4, 3.1, 4.6};
	std::vector<double> mixed_values{5.6, 8.2, 6.8, 4.3};
	
	const CellMatData<2> cell_mat_data(cell_values, mixed_values);
	OUT<CellMatData<2>> out(cell_mat_data);

	const CellMatIndex cell_mat_index(1, 1);
	const std::size_t mixed_index = 2;
	const ValueIndex value_index(ValueIndex::Type::MULTIMAT, mixed_index);

	const double value = 12.3;
	out.get(cell_mat_index, value_index) = value;
	ASSERT_EQ(mixed_values.at(mixed_index), value);
}

TEST(test_REDUCE, cell_data) {
	MultidimArray<2, double> arr({2, 2});
	
	const CellData<2> cell_data(arr);
	auto reducer = [](const double a, const double b) {return a + b;};
	REDUCE<CellData<2>> reduce(reducer, cell_data);

	const double value1 = 2.4;
	const double value2 = 8.2;
        const std::size_t cell_index = 1;
	const CellMatIndex cell_mat_index1(cell_index, 0);
	const CellMatIndex cell_mat_index2(cell_index, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, cell_index);

	reduce.get(cell_mat_index1, value_index) << value1;
	reduce.get(cell_mat_index2, value_index) << value2;
	ASSERT_EQ(arr.at_raw_index(cell_index), reducer(value1, value2));
}

TEST(test_REDUCE, mat_data) {
	std::vector<double> arr(4, 0.0);
	const MatData<> mat_data(arr);
	
	auto reducer = [](const double a, const double b) {return a + b;};
	REDUCE<MatData<>> reduce(reducer, mat_data);

	const double value1 = 2.4;
	const double value2 = 8.2;
        const std::size_t mat_index = 1;
	const CellMatIndex cell_mat_index1(0, mat_index);
	const CellMatIndex cell_mat_index2(1, mat_index);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	reduce.get(cell_mat_index1, value_index) << value1;
	reduce.get(cell_mat_index2, value_index) << value2;
	ASSERT_EQ(arr.at(mat_index), reducer(value1, value2));
}

TEST(test_FREE_SCALAR, free_scalar) {
	const double value = 2.5;
	const FREE_SCALAR<double> free_scalar(value);

	const CellMatIndex cell_mat_index(0, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(free_scalar.get(cell_mat_index, value_index), value);
}

TEST(test_FREE_ARRAY, free_array) {
	const std::vector<double> value = {2.5, 3.0};
	const FREE_ARRAY<double> free_array(value);

	const CellMatIndex cell_mat_index(0, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(free_array.get(cell_mat_index, value_index), value);
}

} // anonymous namespace

} // MM::compressed_cell_centric
