#ifndef MM_ARGUMENTS_HPP
#define MM_ARGUMENTS_HPP

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

#include "Coords.hpp"
#include "Datasets.hpp"

namespace MM {

template<std::size_t N>
void fit_index_to_reduced_data(Coords<N>& index,
			       std::array<std::size_t, N> size) {
	for (std::size_t i = 0; i < N; ++i) {
		if (size[i] == 1) {
			index[i] = 0;
		}
	}
}

template<std::size_t N, typename dtype>
double& unified_data_get(CellData<N, dtype> data,
			 Coords<N> cell_index,
			 const std::size_t mat_index) {
	fit_index_to_reduced_data(cell_index, data.get_size());
	return data.at(cell_index);
}

template<std::size_t N, typename dtype>
double& unified_data_get(MatData<dtype> data,
			 const Coords<N>& cell_index,
			 const std::size_t mat_index) {
	return data.at(mat_index);
}

template<std::size_t N, typename dtype>
double& unified_data_get(CellMatData<N, dtype> data,
			 Coords<N> cell_index,
			 const std::size_t mat_index) {
	fit_index_to_reduced_data(cell_index, data.get_size());
	return data.at(cell_index, mat_index);
}

template<typename T, std::size_t N, typename dtype>
class IN {
public:
	IN(T p_data) : data(p_data)
	{
		static_assert(std::is_same<T, CellData<N, dtype>>::value
			      || std::is_same<T, MatData<dtype>>::value
			      || std::is_same<T, CellMatData<N, dtype>>::value);
	}
	
	double get(const Coords<N>& cell_index, const std::size_t mat_index) {
		return unified_data_get(data, cell_index, mat_index);
	}
	
private:
	T data;
};

template<typename T, std::size_t N, typename dtype>
class OUT {
public:
        OUT(T p_data) : data(p_data)
	{
		static_assert(std::is_same<T, CellData<N, dtype>>::value
			      || std::is_same<T, MatData<dtype>>::value
			      || std::is_same<T, CellMatData<N, dtype>>::value);
	}
	
	double& get(const Coords<N>& cell_index, const std::size_t mat_index) {
		return unified_data_get(data, cell_index, mat_index);
	}
	
private:
	T data;
};


class ReduceProxy {
public:
	ReduceProxy(double& p_reduced_value,
		    std::function<double(double, double)>& p_reducer)
		: reduced_value(p_reduced_value),
		  reducer(p_reducer)
		{
		}
		
	void operator<<(double value) {
		reduced_value = reducer(reduced_value, value);
	}
		
private:
	double& reduced_value;
	std::function<double(double, double)>& reducer;
};

template<typename T, std::size_t N, typename dtype>
class REDUCE {
public:	
	REDUCE(std::function<double(double, double)> p_reducer, T p_data)
		: reducer(p_reducer),
		  data(p_data)
	{
		static_assert(std::is_same<T, CellData<N, dtype>>::value
			      || std::is_same<T, MatData<dtype>>::value
			      || std::is_same<T, CellMatData<N, dtype>>::value);
	}

	ReduceProxy get(const Coords<N>& cell_index,
			const std::size_t mat_index) {
		double& value = unified_data_get(data, cell_index, mat_index);
		return ReduceProxy(value, reducer);
	}

private:
	std::function<double(double, double)> reducer; // A commutative and associative function.
	T data;
};

template <std::size_t N>
class Stencil {
public:
	Stencil(std::vector<Offsets<N>> p_offsets)
		: offsets(p_offsets)
	{
	}

	bool contains_offset(const Offsets<N>& offset) const {
		auto it = std::find(offsets.begin(), offsets.end(), offset);
		return it != offsets.end();
	}
	
private:
	std::vector<Offsets<N>> offsets;
};

template<typename T, std::size_t N, typename dtype>
class NeighProxy {
public:
	NeighProxy(const Coords<N> p_cell_coords, const std::size_t p_mat_index,
		   T p_data, const Stencil<N> p_stencil)
		: cell_coords(p_cell_coords), mat_index(p_mat_index),
		  data(p_data), stencil(p_stencil)
	{
	}

	dtype get_neigh(const Offsets<N>& offset) const {
		if (!stencil.contains_offset(offset)) {
			throw "No such offset in stencil.";
		}
		
		const Coords<N> neighbour_coords = cell_coords + offset;
		return unified_data_get(data, neighbour_coords, mat_index);
	}
	
private:
	const Coords<N> cell_coords;
	const std::size_t mat_index;

	T data;
	const Stencil<N> stencil;	
};

template<typename T, std::size_t N, typename dtype>
class NEIGH {
public:
	NEIGH(T p_data, const Stencil<N> p_stencil)
		: data(p_data), stencil(p_stencil)
	{
		static_assert(std::is_same<T, CellData<N, dtype>>::value
			      || std::is_same<T, MatData<dtype>>::value
			      || std::is_same<T, CellMatData<N, dtype>>::value);
	}

	NeighProxy<T, N, dtype> get(const Coords<N>& cell_index,
				    const std::size_t mat_index) {
		return NeighProxy<T, N, dtype>(cell_index,
					       mat_index,
					       data,
					       stencil);
	}

private:
	const Stencil<N> stencil;
	T data;
};

} // namespace MM

#endif // MM_ARGUMENTS_HPP
