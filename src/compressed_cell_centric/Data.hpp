#ifndef MM_COMPRESSED_CELL_DATA_HPP
#define MM_COMPRESSED_CELL_DATA_HPP

#include <array>
#include <cstddef>
#include <unordered_set>
#include <vector>

#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "compressed_cell_centric/Datasets.hpp"

namespace MM {

namespace compressed_cell_centric {

template <std::size_t N, typename dtype = double>
class Data {
public:
	Data(const std::array<std::size_t, N> p_size,
	     const std::vector<std::vector<std::size_t>>& materials)
		: mat_number(count_materials(materials)),
		  size(p_size),
		  structure(CompressedDataStructure(materials))
	{
		
	}

	const std::array<std::size_t, N>& get_size() const {
		return size;
	}

	std::size_t get_mat_number() const {
		return mat_number;
	}

	// CellData<N, dtype> new_cell_data() {
	// }

	// MatData<dtype> new_mat_data() {
	// }

	// CellMatData<N, dtype> new_cell_mat_data() {
	// }

private:
	static
	std::size_t count_materials(
		const std::vector<std::vector<std::size_t>>& materials) {
		std::unordered_set<std::size_t> set;

		for (const std::vector<std::size_t>& in_cell : materials) {
			for (const std::size_t mat_id : in_cell) {
				set.insert(mat_id);
			}
		}

		return set.size();
	}
	
	const std::size_t mat_number;
	const std::array<std::size_t, N> size;
	CompressedDataStructure structure;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATA_HPP
