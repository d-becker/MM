#ifndef MM_COMPRESSED_CELL_DATA_HPP
#define MM_COMPRESSED_CELL_DATA_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <numeric>
#include <unordered_set>
#include <vector>

#include "mm_assert.hpp"
#include "MultidimArray.hpp"
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
		  structure(CompressedDataStructure(materials)),
		  cell_only_buffers(),
		  mat_only_buffers(),
		  cell_values(),
		  mixed_storage_values()
	{
    const std::size_t total_cell_number = std::accumulate(
        size.begin(),
        size.end(),
        1,
        std::multiplies<std::size_t>());
    MM_ASSERT(total_cell_number == materials.size());
	}

	const std::array<std::size_t, N>& get_size() const {
		return size;
	}

	std::size_t get_mat_number() const {
		return mat_number;
	}

  const CompressedDataStructure get_structure() const {
    return structure;
  }

	CellData<N, dtype> new_cell_data() {
		cell_only_buffers.emplace_back(size);
		return CellData<N, dtype>(cell_only_buffers.back());
	}

	MatData<N, dtype> new_mat_data() {
		mat_only_buffers.emplace_back(mat_number, 0.0);
		return MatData<N, dtype>(mat_only_buffers.back());
	}

	CellMatData<N, dtype> new_cell_mat_data() {
		cell_values.emplace_back(structure.cell_number(), 0.0);
		mixed_storage_values.emplace_back(
			structure.mixed_storage_size(), 0.0);

		return CellMatData<N, dtype>(cell_values.back(),
					     mixed_storage_values.back());
	}

	CompressedDataStructure::CellIteration
	cell_iteration(const Coords<N>& coords) const {
		const std::size_t flat_index
			= coords_to_flat_index(coords, size);
		return structure.cell_iteration(flat_index);
	}

//private:
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

	std::list<MultidimArray<N, dtype>> cell_only_buffers;
  std::list<std::vector<dtype>> mat_only_buffers;
	std::list<std::vector<dtype>> cell_values;
	std::list<std::vector<dtype>> mixed_storage_values;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATA_HPP
