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

    while (index_generator.has_next()) {
      const Coords<N> coords
        = Coords<N>::from_array(index_generator.next());
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
	}

//private:
	Data<N, dtype>& data;
	IndexGenerator<N> index_generator;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_COMPUTATION_HPP
