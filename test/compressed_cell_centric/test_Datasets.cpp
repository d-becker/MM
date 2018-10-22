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

class CompressedCellMatData : public ::testing::Test {
public:
	CompressedCellMatData()
		: cell_values({2.8, 0.0, 0.0, 0.0}),
		  mixed_storage_values({2.6, 8.2, 5.82, 6.38, 1.24, 3.4})
	{
	}
protected:
	std::vector<double> cell_values;
	std::vector<double> mixed_storage_values;
};

TEST_F(CompressedCellMatData, get_cell_number_ok) {
	CellMatData<2> cell_mat_data(cell_values,
				     mixed_storage_values);
	ASSERT_EQ(cell_mat_data.cell_number(), cell_values.size());
}


TEST_F(CompressedCellMatData, get_mixed_storage_size_ok) {
	CellMatData<2> cell_mat_data(cell_values,
				     mixed_storage_values);
	ASSERT_EQ(cell_mat_data.mixed_storage_size(),
		  mixed_storage_values.size());
}

TEST_F(CompressedCellMatData, access_and_modify_cell_value) {
	CellMatData<2> cell_mat_data(cell_values,
				     mixed_storage_values);
	const std::size_t index = 0;
	const double value = cell_values.at(index);
	ASSERT_EQ(cell_mat_data.cell_value_at(index), value);

	const double new_value = value + 1.0;
	cell_mat_data.cell_value_at(index) = new_value;
	ASSERT_EQ(cell_mat_data.cell_value_at(index), new_value);
}

TEST_F(CompressedCellMatData, access_and_modify_mixed_storage_value) {
	CellMatData<2> cell_mat_data(cell_values,
				     mixed_storage_values);
	const std::size_t index = 0;
	const double value = mixed_storage_values.at(index);
	ASSERT_EQ(cell_mat_data.mixed_storage_value_at(index), value);

	const double new_value = value + 1.0;
	cell_mat_data.mixed_storage_value_at(index) = new_value;
	ASSERT_EQ(cell_mat_data.mixed_storage_value_at(index), new_value);
}

} // anonymous namespace

} // MM::compressed_cell_centric
