#ifndef MM_DATA_HPP
#define MM_DATA_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <list>
#include <sstream>
#include <vector>

#include "Coords.hpp"
#include "DimensionError.hpp"
#include "mm_assert.hpp"

#include "full_matrix/Datasets.hpp"

namespace MM {

namespace full_matrix {

template <std::size_t N, typename dtype = double>
class Data {
public:
  Data(const std::array<std::size_t, N> p_size,
       const std::size_t p_mat_number)
    : size(p_size),
      mat_number(p_mat_number)
  {
    // The number of materials cannot be zero.
    MM_ASSERT(p_mat_number != 0);

    // The size cannot contain ones.
    for (const std::size_t d : p_size) {
      MM_ASSERT(d != 1);
    }

    assert_no_zero_size(p_size);
  }

  const std::array<std::size_t, N>& get_size() const {
    return size;
  }

  std::size_t get_mat_number() const {
    return mat_number;
  }

  CellData<N, dtype> new_cell_data() {
          return new_cell_data(size);
  }

  CellData<N, dtype>
  new_cell_data(const std::array<std::size_t, N>& arr_size) {
    assert_no_zero_size(arr_size);
    cell_buffers.emplace_back(arr_size);
    return CellData<N, dtype>(cell_buffers.back());
  }

  MatData<N, dtype> new_mat_data() {
    mat_buffers.emplace_back(mat_number, 0.0);
    return MatData<N, dtype>(mat_buffers.back());
  }

  CellMatData<N, dtype> new_cell_mat_data() {
    return new_cell_mat_data(size);
  }

  CellMatData<N, dtype>
  new_cell_mat_data(const std::array<std::size_t, N>& arr_size) {
    assert_no_zero_size(arr_size);

    std::array<std::size_t, N + 1> extended_arr_size;
    std::copy(arr_size.begin(),
        arr_size.end(),
        extended_arr_size.begin());
    extended_arr_size[N] = mat_number;

    cell_mat_buffers.emplace_back(extended_arr_size);

    return CellMatData<N, dtype>(cell_mat_buffers.back());
  }

private:
  static void assert_no_zero_size(const std::array<std::size_t, N> size) {
    for (const std::size_t d : size) {
      // The size cannot contain zeros.
      MM_ASSERT(d != 1);
    }
  }

  const std::array<std::size_t, N> size;
  const std::size_t mat_number;

  std::list<MultidimArray<N, dtype>> cell_buffers;
  std::list<std::vector<dtype>> mat_buffers;
  std::list<MultidimArray<N + 1, dtype>> cell_mat_buffers;
};

} // namespace full_matrix

} // namespace MM

#endif // MM_DATA_HPP
