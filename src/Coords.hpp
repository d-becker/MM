#ifndef MM_COORDS_HPP
#define MM_COORDS_HPP

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace MM {

template <std::size_t N>
class Coords {
public:
	Coords(std::array<std::size_t, N> p_coordinates)
		: coordinates(p_coordinates)
		{
		}

	template <typename ...Coordinates>
	Coords(Coordinates... coords)
		: coordinates({coords...})
		{
			static_assert(sizeof...(coords) == N);
		}

	std::size_t operator[](const std::size_t n) const {
		return coordinates[n];
	}
	
	std::size_t& operator[](const std::size_t n) {
		return coordinates[n];
	}

	std::size_t at(const std::size_t n) const {
		return coordinates.at(n);
	}

	std::size_t& at(const std::size_t n) {
		return coordinates.at(n);
	}
private:
	std::array<std::size_t, N> coordinates;
};

} // namespace MM

#endif // MM_COORDS_HPP
