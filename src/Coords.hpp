#ifndef MM_COORDS_HPP
#define MM_COORDS_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>


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
    // assert(n < N);
    return coordinates.at(n);
		return coordinates[n];
	}

  dtype& operator[](const std::size_t n) {
    assert(n < N);
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

} // namespace MM

#endif // MM_COORDS_HPP
