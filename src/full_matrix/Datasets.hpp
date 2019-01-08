#ifndef MM_DATASETS_HPP
#define MM_DATASETS_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <list>
#include <vector>

#include "Coords.hpp"
#include "MultidimArray.hpp"

namespace MM {

namespace full_matrix {

template<std::size_t _N, typename _dtype = double>
class CellData {
public:
	constexpr static std::size_t N = _N;
	using dtype = _dtype;
	
	CellData(MultidimArray<N, dtype>& p_buffer)
		: buffer(p_buffer)
	{
	}

	const std::array<std::size_t, N>& get_size() const {
		return buffer.get_size();
	}

	std::size_t get_flat_size() const {
		return buffer.get_flat_size();
	}
	
	const dtype& operator[](const Coords<N>& index) const {
		return buffer[index];
	}

	dtype& operator[](const Coords<N>& index) {
		return buffer[index];
	}

	const dtype& at(const Coords<N>& index) const {
		return buffer.at(index);
	}

	dtype& at(const Coords<N>& index) {
		return buffer.at(index);
	}

	const dtype& at_raw_index(const std::size_t index) const {
		return buffer.at_raw_index(index);
	}

	dtype& at_raw_index(const std::size_t index) {
		return buffer.at_raw_index(index);
	}

  dtype* get_raw(std::array<std::size_t, N> &shape) {
    shape = buffer.get_size();
    return &buffer.at_raw_index(0);
  }

private:
	MultidimArray<N, dtype>& buffer;
};

template<std::size_t _N, typename _dtype = double>
class MatData {
public:
	using dtype = _dtype;
	constexpr static std::size_t N = _N;
	
	MatData(std::vector<dtype>& p_material_data)
	        : material_data(p_material_data)
	{
	}

	std::size_t get_size() const {
		return material_data.size();
	}
	
	const dtype& operator[](const std::size_t index) const {
		return material_data[index];
	}

	dtype& operator[](const std::size_t index) {
		return material_data[index];
	}

	const dtype& at(const std::size_t index) const {
		return material_data.at(index);
	}

	dtype& at(const std::size_t index) {
		return material_data.at(index);
	}

  dtype* get_raw(std::array<std::size_t,N> &shape) {
    std::array<std::size_t,N> _shape = {material_data.size()};
    shape = _shape;
    return &material_data[0];
  }

private:
	std::vector<dtype>& material_data;
};

template<std::size_t _N, typename _dtype = double>
class CellMatData {
public:
	constexpr static std::size_t N = _N;
	using dtype = _dtype;
	
	CellMatData(MultidimArray<N + 1, dtype>& p_data)
		: data(p_data),
		  size_without_mat_dimension(
			  cut_last_dimension(p_data.get_size()))
	{
	}

	const std::array<std::size_t, N>& get_size() const {
		return size_without_mat_dimension;
	}

	const dtype& get_unchecked(const Coords<N> cell_index,
				   const std::size_t mat_index) const {
		return data[extend_coords(cell_index, mat_index)];
	}

	dtype& get_unchecked(const Coords<N> cell_index,
			     const std::size_t mat_index) {
		return data[extend_coords(cell_index, mat_index)];
	}

	const dtype& at(const Coords<N> cell_index,
			const std::size_t mat_index) const {
		return data.at(extend_coords(cell_index, mat_index));
	}

	dtype& at(const Coords<N> cell_index,
		  const std::size_t mat_index) {
		return data.at(extend_coords(cell_index, mat_index));
	}

  dtype* get_raw(std::array<std::size_t, N> &shape) {
    shape = size_without_mat_dimension;
    return &data.at_raw_index(0);
  }

private:
	static std::array<std::size_t, N>
	cut_last_dimension(const std::array<std::size_t, N + 1>& arr) {
		std::array<std::size_t, N> res;
		std::copy(arr.begin(), arr.end() - 1, res.begin());

		return res;
	}

	static Coords<N + 1>
	extend_coords(const Coords<N> coords, const std::size_t mat_index) {
		std::array<std::size_t, N + 1> arr;

		for (std::size_t i = 0; i < N; ++i) {
			arr[i] = coords[i];
		}

		arr[N] = mat_index;

		return Coords<N + 1>::from_array(arr);
	}
	
	// The last dimension is the material.
	MultidimArray<N + 1, dtype>& data;

	std::array<std::size_t, N> size_without_mat_dimension;
};

} // namespace full_matrix

} // namespace MM

#endif // MM_DATASETS_HPP
