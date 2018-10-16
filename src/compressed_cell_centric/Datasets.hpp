#ifndef MM_COMPRESSED_CELL_DATASETS_HPP
#define MM_COMPRESSED_CELL_DATASETS_HPP

#include <cstddef>

#include "full_matrix/Datasets.hpp"

namespace MM {

namespace compressed_cell_centric {

template<std::size_t _N, typename _dtype = double>
using CellData = full_matrix::CellData<_N, _dtype>;

template<typename _dtype = double>
using MatData = full_matrix::MatData<_dtype>;

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATASETS_HPP
