#include <array>
#include <cstddef>
#include <cstdlib>
#include <set>
#include <utility>
#include <vector>


#include "gtest/gtest.h"

#include "compressed_cell_centric/CompressedDataStructure.hpp"

namespace MM::compressed_cell_centric {

namespace {

void check_structure_iteration(
	const CompressedDataStructure& structure,
	const std::vector<std::vector<std::size_t>>& raw_data)
{
	for (std::size_t cell_index = 0;
	     cell_index < raw_data.size();
	     ++cell_index) {
		const std::vector<std::size_t> materials
			= raw_data.at(cell_index);
		CompressedDataStructure::CellIteration iteration
			= structure.cell_iteration(cell_index);
		CompressedDataStructure::CellIterator it = iteration.begin();
		for (std::size_t mat_index = 0;
		     mat_index < materials.size();
		     ++mat_index) {
			const std::pair<CellMatIndex, ValueIndex> pair = *it;
			const CellMatIndex& cell_mat_index = pair.first;
			ASSERT_EQ(cell_mat_index.cell_index, cell_index);
			ASSERT_EQ(cell_mat_index.mat_index,
				  materials.at(mat_index));

			const ValueIndex& value_index = pair.second;
			const Cell& cell = structure.cell_at(cell_index);

			if (cell.nmats > 1) {
				ASSERT_EQ(value_index.type,
					  ValueIndex::Type::MULTIMAT);
				ASSERT_EQ(structure.mixed_cell_at(
						  value_index.index).material,
					  materials.at(mat_index));
			} else {
				ASSERT_EQ(value_index.type,
					  ValueIndex::Type::SINGLE_MAT);
				ASSERT_EQ(value_index.index, cell_index);
			}

			++it;
		}

		ASSERT_EQ(it, iteration.end());
	}
}

std::vector<std::vector<std::size_t>> get_raw_data() {
	return {
		std::vector<std::size_t>{1},
		std::vector<std::size_t>{1, 2, 3},
		std::vector<std::size_t>{4},
		std::vector<std::size_t>{2, 4},
		std::vector<std::size_t>{3}
	};
}

TEST(CompressedDataStructure, create_and_check_data_structure) {
	const std::vector<std::vector<std::size_t>> arr = get_raw_data();
	const CompressedDataStructure structure(arr);

	check_structure_iteration(structure, arr);
}

} // anonymous namespace

} // MM::compressed_cell_centric
