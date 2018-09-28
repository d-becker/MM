#ifndef MM_DATA_HPP
#define MM_DATA_HPP

#include <array>
#include <cstddef>
#include <list>
#include <sstream>
#include <vector>

#include "Coords.hpp"
#include "Datasets.hpp"
#include "DimensionError.hpp"

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
		if (p_mat_number == 0) {
			throw DimensionError(
				"The number of materials cannot be zero.");
		}

		for (const std::size_t d : p_size) {
			if (d == 1) {
				throw DimensionError(
					"The size cannot contain ones.");
			}
		}

		throw_on_zero_size(p_size);
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
		throw_on_zero_size(arr_size);
		cell_buffers.emplace_back(arr_size);
		return CellData<N, dtype>(cell_buffers.back());
	}

	MatData<dtype> new_mat_data() {
		mat_buffers.emplace_back(mat_number, 0.0);
		return MatData<dtype>(mat_buffers.back());
	}

	CellMatData<N, dtype> new_cell_mat_data() {
		return new_cell_mat_data(size);
	}
	
	CellMatData<N, dtype>
	new_cell_mat_data(const std::array<std::size_t, N>& arr_size) {
		throw_on_zero_size(arr_size);
		cell_mat_buffers.emplace_back(
			mat_number,
			MultidimArray<N, dtype>(arr_size));
		
		return CellMatData<N, dtype>(cell_mat_buffers.back());
	}
	
private:
	static void throw_on_zero_size(const std::array<std::size_t, N> size) {
		for (const std::size_t d : size) {
			if (d == 0) {
				throw DimensionError(
					"The size cannot contain zeros.");
			}
		}
	}
	
	const std::array<std::size_t, N> size;
	const std::size_t mat_number;
	
	std::list<MultidimArray<N, dtype>> cell_buffers;
	std::list<std::vector<dtype>> mat_buffers;
	std::list<std::vector<MultidimArray<N, dtype>>> cell_mat_buffers;
};

} // namespace full_matrix

} // namespace MM

#endif // MM_DATA_HPP
