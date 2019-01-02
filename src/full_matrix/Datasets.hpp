#ifndef MM_DATASETS_HPP
#define MM_DATASETS_HPP

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

  std::vector<dtype*> get_raw(std::array<std::size_t, N> &shape) {
    shape = buffer.get_size();
    std::vector<dtype*> arr = {&buffer.at_raw_index(0)};
    return arr;
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

  std::vector<dtype*> get_raw(std::array<std::size_t,N> &shape) {
    std::array<std::size_t,1> _shape = {material_data.size()};
    shape = _shape;
    std::vector<dtype*> arr = {&material_data[0]};
    return arr;
  }

private:
	std::vector<dtype>& material_data;
};

template<std::size_t _N, typename _dtype = double>
class CellMatData {
public:
	constexpr static std::size_t N = _N;
	using dtype = _dtype;
	
	CellMatData(std::vector<MultidimArray<N, dtype>>& p_data)
		: data(p_data)
	{
	}

	const std::array<std::size_t, N>& get_size() const {
		return data.at(0).get_size();
	}

	const dtype& get_unchecked(const Coords<N> cell_index,
				   const std::size_t mat_index) const {
		return data[mat_index][cell_index];
	}

	dtype& get_unchecked(const Coords<N> cell_index,
			     const std::size_t mat_index) {
		return data[mat_index][cell_index];
	}

	const dtype& at(const Coords<N> cell_index,
			const std::size_t mat_index) const {
		return data.at(mat_index).at(cell_index);
	}

	dtype& at(const Coords<N> cell_index,
		  const std::size_t mat_index) {
		return data.at(mat_index).at(cell_index);
	}


  std::vector<dtype*> get_raw(std::array<std::size_t, N> &shape) {
    shape = data.at(0).get_size();
    std::vector<dtype*> arr(data.size());
    for (size_t i = 0; i < data.size(); i++)
      arr[i] = &(data[i].at_raw_index(0)); 
    return arr;
  }
	
private:
	std::vector<MultidimArray<N, dtype>>& data;
	
};

} // namespace full_matrix

} // namespace MM

#endif // MM_DATASETS_HPP
