#ifndef MM_COMPRESSED_CELL_COMPUTATION_HPP
#define MM_COMPRESSED_CELL_COMPUTATION_HPP

#include <type_traits>

#include "IndexGenerator.hpp"
#include "compressed_cell_centric/Arguments.hpp"
#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "compressed_cell_centric/Data.hpp"

namespace MM {

namespace compressed_cell_centric {

template <std::size_t N, typename dtype = double>
class Computation {
public:
	Computation(Data<N, dtype>& p_data,
		    IndexGenerator<N> p_index_generator)
		: data(p_data),
		  index_generator(p_index_generator)
	{
	}

	template<typename FUNC, typename ...ARGS>
	void compute(FUNC func, ARGS ...args) {
		// TODO: check all args belong to this->data.
#ifdef MM_FUSED
  const CompressedDataStructure structure = data.get_structure();
  const std::array<std::size_t, N>& size = data.get_size();

  // Iterate over single-material cells.
  for (std::size_t i = 0; i < structure.cell_number(); i++) {
    const Cell& cell = structure.cell_at(i);
    const Coords<N> coords = flat_index_to_coords(i, size);

    if (index_generator.is_in_range(coords) && cell.nmats == 1) {
      const std::size_t mat_id = cell.imat;
      const CellMatIndex cell_mat_index(i, mat_id);
      const ValueIndex value_index(ValueIndex::Type::SINGLE_MAT, i);

      func(args.get(coords, data, cell_mat_index, value_index)...);
    }
  }

  // Iterate through the mixed cells.
  for (std::size_t i = 0; i < structure.mixed_storage_size(); i++) {
    const MixedStorageCell& mixed_cell = structure.mixed_cell_at(i);
    const std::size_t cell_id = mixed_cell.frac2cell;
    const Coords<N> coords = flat_index_to_coords(cell_id, size);

    if (index_generator.is_in_range(coords)) {
      const CellMatIndex cell_mat_index(cell_id, mixed_cell.material);
      const ValueIndex value_index(ValueIndex::Type::MULTIMAT, i);

      func(args.get(coords, data, cell_mat_index, value_index)...);
    }
  }
#else
  while (index_generator.has_next()) {
    const Coords<N> coords = index_generator.next();
    for (std::pair<CellMatIndex, ValueIndex> pair
        : data.cell_iteration(coords)) {
      const CellMatIndex& cell_mat_index
        = pair.first;
      const ValueIndex& value_index
        = pair.second;
      func(args.get(coords,
            data,
            cell_mat_index,
            value_index)...);
    }
  }
#endif
	}

//private:
	Data<N, dtype>& data;
	IndexGenerator<N> index_generator;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_COMPUTATION_HPP
