#ifndef MM_COMPRESSED_DATA_STRUCTURE_HPP
#define MM_COMPRESSED_DATA_STRUCTURE_HPP

#include <vector>

#include "MultidimArray.hpp"

namespace MM {

namespace compressed_cell_centric {

template <std::size_t N, typename dtype>
class CompressedDataStructure {
public:

private:
	struct Cell {
		long long nmats;
		long long imat;
	};

	struct MixedStorageCell {
		long long nextfrac;
		long long material;
	};

	MultidimArray<N, Cell> structure;
	std::vector<MixedStorageCell> mixed_storage;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_DATA_STRUCTURE_HPP
