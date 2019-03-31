#ifndef MM_COORDS_HPP
#define MM_COORDS_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "mm_assert.hpp"

namespace MM {

template <typename dtype, std::size_t N>
class VecN {
public:
	static VecN from_array(std::array<dtype, N> p_coordinates)
	{
		VecN res;
		res.coordinates = p_coordinates;
		return res;
	}

	template <typename ...Coordinates>
	VecN(Coordinates... coords)
		: coordinates({coords...})
	{
		static_assert(sizeof...(coords) == N, "Wrong number of coordinates.");
	}

	dtype operator[](const std::size_t n) const {
    MM_ASSERT(n < N);
		return coordinates[n];
	}

  dtype& operator[](const std::size_t n) {
    MM_ASSERT(n < N);
		return coordinates[n];
	}

	bool operator==(const VecN<dtype, N>& other) const {
		return this->coordinates == other.coordinates;
	}
private:
	VecN()
    : coordinates()
  {
    coordinates.fill(0);
  }

	std::array<dtype, N> coordinates;
};

template <std::size_t N>
using Coords = VecN<std::size_t, N>;

template <std::size_t N>
using Offsets = VecN<long long int, N>;

template <std::size_t N>
Coords<N> operator+(const Coords<N>& coords, const Offsets<N>& offsets) {
	std::array<std::size_t, N> c;

	for (std::size_t i = 0; i < N; ++i) {
		c[i] = coords[i] + offsets[i];
	}

	return Coords<N>::from_array(c);
}

template <std::size_t N>
std::size_t coords_to_flat_index(const Coords<N>& coords,
                                 const std::array<std::size_t, N>& size) {
  std::size_t multiplier = 1;
  std::size_t index = 0;

  for (std::size_t i = 0; i < N; ++i) {
    index += coords[i] * multiplier;
    multiplier *= size[i];
  }

  return index;
}

template <std::size_t N>
Coords<N> flat_index_to_coords(const std::size_t index,
                               const std::array<std::size_t, N>& size) {
  std::array<std::size_t, N> partial_products;
  partial_products[0] = 1;
  for (std::size_t i = 0; i < N - 1; i++) {
    partial_products[i + 1] = partial_products[i] * size[i];
  }

  std::size_t remaining_index = index;
  std::array<std::size_t, N> res;
  for (std::size_t i = N - 1; (i + 1) > 0; i--) {
    res[i] = remaining_index / partial_products[i];
    remaining_index = remaining_index % partial_products[i];
  }


  return Coords<N>::from_array(res);
}

} // namespace MM

#endif // MM_COORDS_HPP
