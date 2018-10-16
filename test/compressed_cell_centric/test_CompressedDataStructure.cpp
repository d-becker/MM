#include <array>
#include <cstddef>
#include <cstdlib>
#include <set>
#include <vector>


#include "gtest/gtest.h"

#include "compressed_cell_centric/CompressedDataStructure.hpp"

#include "Util.hpp"

namespace MM::compressed_cell_centric {

namespace {

void check_cell(const std::size_t cell_index,
	        const CompressedDataStructure& structure,
		const std::vector<std::size_t>& materials) {
	const Cell& cell = structure.cell_at(cell_index);
	ASSERT_EQ(cell.nmats, materials.size());
	
	if (materials.size() <= 1) {
		ASSERT_EQ(cell.imat, materials.at(0));
	} else {
		std::set<std::size_t> materials_from_structure;
		for (std::size_t material_index : structure.mixed_mat_iteration(cell_index)) {
			const MixedStorageCell& mixed_cell
				= structure.mixed_cell_at(material_index);
			ASSERT_EQ(mixed_cell.frac2cell, cell_index);
			materials_from_structure.insert(mixed_cell.material);
		}

		const std::set<std::size_t> raw_materials_as_set(materials.begin(),
								 materials.end());
		ASSERT_EQ(materials_from_structure, raw_materials_as_set);
	}
}

void check_cells(const std::vector<std::vector<std::size_t>>& raw_data,
		 const CompressedDataStructure& structure) {
	for (std::size_t i = 0; i < raw_data.size(); ++i) {
		const std::vector<std::size_t>& materials
			= raw_data.at(i);
		check_cell(i, structure, materials);
	}
}

TEST(CompressedDataStructure, create_and_check_data_structure) {
	constexpr std::size_t N = 5;
	constexpr std::size_t M = 5;
	const std::vector<std::vector<std::size_t>> arr = get_raw_data(N, M);
	const CompressedDataStructure structure(arr);

	check_cells(arr, structure);
}

TEST(CompressedDataStructure,
     mixed_mat_iteration_on_single_mat_cell_empty_iterator) {
	const std::vector<std::vector<std::size_t>> materials = {
		std::vector<std::size_t>{1}
	};

	const CompressedDataStructure structure(materials);

	auto iteration = structure.mixed_mat_iteration(0);

	ASSERT_EQ(iteration.begin(), iteration.end());
}

} // anonymous namespace

} // MM::compressed_cell_centric
