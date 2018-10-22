#include <array>
#include <vector>

#include "gtest/gtest.h"

#include "compressed_cell_centric/Data.hpp"

#include "Util.hpp"

namespace MM::compressed_cell_centric {

namespace {

TEST(CompressedData, size_correct) {
	const std::vector<std::vector<std::size_t>> arr
		= get_raw_data(10, 4);
	const std::array<std::size_t, 2> size = {4, 8};

	Data<2> data(size, arr);
	ASSERT_EQ(data.get_size(), size);
}

TEST(CompressedData, material_number_correct) {
	constexpr std::size_t cell_no = 4;
	constexpr std::size_t mat_no = 3;

	const std::vector<std::vector<std::size_t>> arr
		= get_raw_data(cell_no, mat_no);

	Data<2> data({2, 2}, arr);
	ASSERT_EQ(data.get_mat_number(), mat_no);
}

TEST(CompressedData, new_cell_data) {
	constexpr std::size_t cell_no = 4;
	constexpr std::size_t mat_no = 3;

	const std::vector<std::vector<std::size_t>> arr
		= get_raw_data(cell_no, mat_no);

	Data<2> data({2, 2}, arr);
	CellData<2> cell_data = data.new_cell_data();
	ASSERT_EQ(cell_data.get_size(), data.get_size());
}

TEST(CompressedData, new_mat_data) {
	constexpr std::size_t cell_no = 4;
	constexpr std::size_t mat_no = 3;

	const std::vector<std::vector<std::size_t>> arr
		= get_raw_data(cell_no, mat_no);

	Data<2> data({2, 2}, arr);
	MatData<double> mat_data = data.new_mat_data();
	ASSERT_EQ(mat_data.get_size(), data.get_mat_number());
}

// TODO: Decide if we delete it or add get_size to CellMatData.
// TEST(CompressedData, new_cell_mat_data) {
// 	constexpr std::size_t cell_no = 4;
// 	constexpr std::size_t mat_no = 3;

// 	const std::vector<std::vector<std::size_t>> arr
// 		= get_raw_data(cell_no, mat_no);

// 	Data<2> data({2, 2}, arr);
// 	CellMatData<2, double> cell_mat_data = data.new_cell_mat_data();
// 	ASSERT_EQ(cell_mat_data.get_size(), data.get_size());
// }

} // anonymous namespace

} // MM::compressed_cell_centric
