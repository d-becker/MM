#include <array>

#include "gtest/gtest.h"

#include "Coords.hpp"
#include "MultidimArray.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/Data.hpp"
#include "compressed_cell_centric/Datasets.hpp"

namespace MM::compressed_cell_centric {

namespace {

Data<2, double> dummy_data() {
  return Data<2, double>({2, 2}, {{}, {}, {}, {}});
}

TEST(test_IN, cell_data) {
	MultidimArray<2, double> arr({2, 2});
	const double value = 2.4;
	const std::size_t cell_index = 1;
	arr.at_raw_index(cell_index) = value;

	const CellData<2> cell_data(arr);
	const IN<CellData<2>> in(cell_data);

	const CellMatIndex cell_mat_index(cell_index, 0);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(in.get(dummy_data(), cell_mat_index, value_index), value);
}

TEST(test_IN, mat_data) {
	std::vector<double> arr(4, 0.0);
	const double value = 2.4;
	arr.at(1) = value;

	constexpr std::size_t N = 2;
	const MatData<N> mat_data(arr);
	const IN<MatData<N>> in(mat_data);

	const CellMatIndex cell_mat_index(0, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(in.get(dummy_data(), cell_mat_index, value_index), value);
}

TEST(test_IN, cell_mat_data) {
	std::vector<double> cell_values{2.2, 1.4, 3.1, 4.6};
	std::vector<double> mixed_values{5.6, 8.2, 6.8, 4.3};

	const CellMatData<2> cell_mat_data(cell_values, mixed_values);
	const IN<CellMatData<2>> in(cell_mat_data);

	const CellMatIndex cell_mat_index(1, 1);
	const std::size_t mixed_index = 2;
	const ValueIndex value_index(ValueIndex::Type::MULTIMAT, mixed_index);

  ASSERT_EQ(in.get(dummy_data(), cell_mat_index, value_index),
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

	out.get(dummy_data(), cell_mat_index, value_index) = value;
  ASSERT_EQ(arr.at_raw_index(cell_index), value);
}

TEST(test_OUT, mat_data) {
	std::vector<double> arr(4, 0.0);

	constexpr std::size_t N = 2;
	const MatData<N> mat_data(arr);
	OUT<MatData<N>> out(mat_data);

	const double value = 2.4;
        const std::size_t mat_index = 1;
	const CellMatIndex cell_mat_index(0, mat_index);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	out.get(dummy_data(), cell_mat_index, value_index) = value;
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
	out.get(dummy_data(), cell_mat_index, value_index) = value;
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

	reduce.get(dummy_data(), cell_mat_index1, value_index) << value1;
	reduce.get(dummy_data(), cell_mat_index2, value_index) << value2;
	ASSERT_EQ(arr.at_raw_index(cell_index), reducer(value1, value2));
}

TEST(test_REDUCE, mat_data) {
	constexpr std::size_t N = 2;
	std::vector<double> arr(4, 0.0);
	const MatData<N> mat_data(arr);

	auto reducer = [](const double a, const double b) {return a + b;};
	REDUCE<MatData<N>> reduce(reducer, mat_data);

	const double value1 = 2.4;
	const double value2 = 8.2;
  const std::size_t mat_index = 1;
	const CellMatIndex cell_mat_index1(0, mat_index);
	const CellMatIndex cell_mat_index2(1, mat_index);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	reduce.get(dummy_data(), cell_mat_index1, value_index) << value1;
	reduce.get(dummy_data(), cell_mat_index2, value_index) << value2;
	ASSERT_EQ(arr.at(mat_index), reducer(value1, value2));
}

TEST(test_NEIGH, has_neigh) {
  Stencil<2> stencil({{0, 0}, {0, 1}});

  MultidimArray<2, double> arr({2, 2});
  const CellData<2> cell_data(arr);

	NEIGH<CellData<2>> neigh(cell_data, stencil);

	NeighProxy<CellData<2>> proxy = neigh.get(
		dummy_data(),
		CellMatIndex(0, 0),
		ValueIndex(ValueIndex::Type::SINGLE_MAT, 0));

	ASSERT_TRUE(proxy.has_neigh({0, 0}));
	ASSERT_TRUE(proxy.has_neigh({0, 1}));
	ASSERT_FALSE(proxy.has_neigh({1, 0}));
}

TEST(test_NEIGH, get_neigh) {
	Stencil<2> stencil({{0, 0}, {0, 1}});

  MultidimArray<2, double> arr({2, 2});
	const double value1 = 2.0;
	const double value2 = 4.0;
	arr[Coords<2>(0u, 0u)] = value1;
	arr[Coords<2>(0u, 1u)] = value2;
	const CellData<2> cell_data(arr);

	NEIGH<CellData<2>> neigh(cell_data, stencil);

	NeighProxy<CellData<2>> proxy = neigh.get(
		dummy_data(),
		CellMatIndex(0, 0),
		ValueIndex(ValueIndex::Type::SINGLE_MAT, 0));

	Data<2, double> data = dummy_data();
	ASSERT_EQ(value1, proxy.get_neigh({0, 0}));
	ASSERT_EQ(value2, proxy.get_neigh({0, 1}));
}

TEST(test_INDEX, index) {
	const Coords<2> coords(0u, 1u);
  const Data<2>& data = dummy_data();
  const std::size_t cell_index = coords_to_flat_index(coords, data.get_size());
	INDEX<2> index;

  const Coords<2> result = index.get(data, CellMatIndex(cell_index, 0),
      ValueIndex(ValueIndex::Type::SINGLE_MAT, 0));
	ASSERT_EQ(coords, result);
}

TEST(test_FREE_SCALAR, free_scalar) {
	const double value = 2.5;
	const FREE_SCALAR<double> free_scalar(value);

	const CellMatIndex cell_mat_index(0, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(free_scalar.get(dummy_data(), cell_mat_index, value_index), value);
}

TEST(test_FREE_ARRAY, free_array) {
	const std::vector<double> value = {2.5, 3.0};
	const FREE_ARRAY<double> free_array(value);

	const CellMatIndex cell_mat_index(0, 1);
	const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, 0);

	ASSERT_EQ(free_array.get(dummy_data(), cell_mat_index, value_index), value);
}

} // anonymous namespace

} // MM::compressed_cell_centric
