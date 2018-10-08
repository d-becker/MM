#include <array>
#include <cstddef>
#include <cstdlib>
#include <set>
#include <vector>


#include "gtest/gtest.h"

#include "compressed_cell_centric/CompressedDataStructure.hpp"

namespace MM::compressed_cell_centric {

namespace {

template <std::size_t N>
std::array<std::vector<std::size_t>, N> get_raw_data() {
	constexpr std::size_t M = 5;

	srand(0); // Deterministic pseudo-random number generation.

	std::array<std::vector<std::size_t>, N> res;

	for (std::size_t i = 0; i < N; ++i) {
		const std::size_t length = rand() % M;

		res[i] = std::vector<std::size_t>(length);

		for (std::size_t& material : res[i]) {
			material = rand() % M; // TODO: Avoid duplicate materials.
		}
	}

	return res;
}

template <std::size_t N>
void check_cell(const std::size_t cell_index,
	        const CompressedDataStructure<N>& structure,
		const std::vector<std::size_t>& materials) {
	const Cell& cell = structure.cell_at(cell_index);
	
	ASSERT_EQ(cell.nmats, materials.size());

	std::set<std::size_t> materials_from_structure;
	for (std::size_t material_index : structure.iteration(cell_index)) {
		const MixedStorageCell& mixed_cell
			= structure.mixed_cell_at(material_index);
		ASSERT_EQ(mixed_cell.frac2cell, cell_index);
		materials_from_structure.insert(mixed_cell.material);
	}

	const std::set<std::size_t> raw_materials_as_set(materials.begin(),
							 materials.end());
	ASSERT_EQ(materials_from_structure, raw_materials_as_set);
}

template <std::size_t N>
void check_cells(const std::array<std::vector<std::size_t>, N>& raw_data,
		 const CompressedDataStructure<N>& structure) {
	for (std::size_t i = 0; i < N; ++i) {
		const std::vector<std::size_t>& materials
			= raw_data.at(i);
		check_cell(i, structure, materials);
	}
}

TEST(CompressedDataStructure, create_data_structure) {
	constexpr std::size_t N = 5;
	
	const std::array<std::vector<std::size_t>, N> arr = get_raw_data<N>();
	const CompressedDataStructure<N> structure(arr);

	check_cells(arr, structure);
}

} // anonymous namespace

} // MM::compressed_cell_centric
