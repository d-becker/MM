#ifndef MM_MULTIDIM_ARRAY_HPP
#define MM_MULTIDIM_ARRAY_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

#include "Coords.hpp"

namespace MM {

template <std::size_t N>
std::size_t multidim_index_to_raw(const Coords<N>& coords,
                                  const std::array<std::size_t, N> size) {
  std::size_t multiplier = 1;
  std::size_t index = 0;

  for (std::size_t i = 0; i < N; ++i) {
    index += coords[i] * multiplier;
    multiplier *= size[i];
  }

  return index;
}

template <std::size_t N, typename dtype>
class MultidimArray {
public:
  MultidimArray(const std::array<std::size_t, N> p_size)
    : size(p_size),
      buffer(calculate_buffer_length(p_size), 0.0)
  {
  }

  const std::array<std::size_t, N>& get_size() const {
    return size;
  }

  std::size_t get_flat_size() const {
    return buffer.size();
  }

  const dtype& operator[](const Coords<N>& coords) const {
    assert(bounds_check(coords));
    const std::size_t index = multidim_index_to_raw(coords, size);
    return buffer[index];
  }

  dtype& operator[](const Coords<N>& coords) {
    const MultidimArray<N, dtype>& const_this
      = const_cast<const MultidimArray<N, dtype>&>(*this);
    const dtype& const_result = const_this[coords];
    return const_cast<dtype&>(const_result);
  }

  const dtype& at_raw_index(std::size_t index) const {
    assert(index < buffer.size());
    return buffer[index];
  }

  dtype& at_raw_index(std::size_t index) {
    assert(index < buffer.size());
    return buffer[index];
  }
  
private:
  static std::size_t
  calculate_buffer_length(const std::array<std::size_t, N> size) {
    std::size_t res = 1;

    for (const std::size_t d : size) {
      res *= d;
    }

    return res;
  }
  
  
  bool bounds_check(const Coords<N>& coords) const {
    for (std::size_t i = 0; i < N; ++i) {
      if (coords[i] >= size[i]) {
        return false;
      }
    }

    return true;
  }
  
  const std::array<std::size_t, N> size;
  std::vector<dtype> buffer;
};

} // namespace MM

#endif // MM_MULTIDIM_ARRAY_HPP
