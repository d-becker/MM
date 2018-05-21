#ifndef MM_DATA_HPP
#define MM_DATA_HPP

#include <array>
#include <cstddef>
#include <list>
#include <vector>

#include "Coords.hpp"
#include "Datasets.hpp"

namespace MM {

template <std::size_t N, typename dtype = double>
class Data {
public:
	Data(const std::array<std::size_t, N> p_size,
	     const std::size_t p_mat_number)
		: size(p_size),
		  mat_number(p_mat_number)
	{
	}
	
	const std::array<std::size_t, N>& get_size() const {
		return size;
	}

	CellData<N, dtype> new_cell_data() {
		const std::size_t buffer_length = calculate_buffer_length();
		cell_buffers.emplace_back(buffer_length, 0.0);
		return CellData<N, dtype>(size, cell_buffers.back());
	}

	MatData<dtype> new_mat_data() {
		mat_buffers.emplace_back(mat_number, 0.0);
		return MatData<dtype>(mat_buffers.back());
	}

	CellMatData<N, dtype> new_cell_mat_data() {
		
	}
	
private:
	std::size_t calculate_buffer_length() {
		std::size_t res = 1;

		for (const std::size_t d : size) {
			res *= d;
		}

		return res;
	}
	
	const std::array<std::size_t, N> size;
	const std::size_t mat_number;
	
	std::list<std::vector<dtype>> cell_buffers;
	std::list<std::vector<dtype>> mat_buffers;
	std::list<std::vector<std::vector<dtype>> cell_mat_buffers;
};

} // namespace MM

#endif // MM_DATA_HPP
