#ifndef MM_DATASETS_HPP
#define MM_DATASETS_HPP

#include <array>
#include <cstddef>
#include <list>
#include <vector>

#include "Coords.hpp"
#include "MultidimArray.hpp"

namespace MM {

template<std::size_t N, typename dtype = double>
class CellData : private MultidimArray<N, dtype> {
public:
	CellData(const std::array<std::size_t, N> p_size,
		 std::vector<dtype>& p_buffer)
		: MultidimArray<N, dtype>(p_size, p_buffer)
	{
	}

	using MultidimArray<N, dtype>::get_size;
	using MultidimArray<N, dtype>::at;
	using MultidimArray<N, dtype>::operator[];
};

template<typename dtype = double>
class MatData {
public:
	MatData(std::vector<dtype>& p_material_data)
	        : material_data(p_material_data)
	{
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
	
private:
	std::vector<dtype>& material_data;
};

template<std::size_t N, typename dtype = double>
class CellMatData {
public:
	CellMatData(std::vector<MultidimArray<N, dtype>> p_data)
		: data(p_data)
		{
		}
private:
	std::vector<MultidimArray<N, dtype>> data;
	
};

} // namespace MM

#endif // MM_DATASETS_HPP
