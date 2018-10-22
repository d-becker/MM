#ifndef MM_ARGUMENTS_HPP
#define MM_ARGUMENTS_HPP

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

#include "Coords.hpp"
#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "compressed_cell_centric/Datasets.hpp"

namespace MM {

namespace compressed_cell_centric {

template<std::size_t N, typename dtype>
dtype& unified_data_get(CellData<N, dtype> data,
		        const CellMatIndex& cell_mat_index,
			const ValueIndex& _value_index) {
	return data.at_raw_index(cell_mat_index.cell_index);
}

template<typename dtype>
dtype& unified_data_get(MatData<dtype> data,
		        const CellMatIndex& cell_mat_index,
			const ValueIndex& _value_index) {
	return data.at(cell_mat_index.mat_index);
}

template<std::size_t N, typename dtype>
dtype& unified_data_get(CellMatData<N, dtype> data,
		        const CellMatIndex& cell_mat_index,
			const ValueIndex& value_index) {
	if (value_index.type == ValueIndex::Type::MULTIMAT) {
		return data.mixed_storage_value_at(value_index.index);
	}
	
	return data.cell_value_at(value_index.index);
}

template<typename T>
class IN {
public:
	using dtype = typename T::dtype;
	
	IN(T p_data) : data(p_data)
	{
                // TODO: Solve this.
		// constexpr std::size_t N = T::N;
		
		// static_assert(std::is_same<T, CellData<N, dtype>>::value
		// 	      || std::is_same<T, MatData<dtype>>::value
		// 	      || std::is_same<T, CellMatData<N, dtype>>::value);
	}
	
	const dtype& get(const CellMatIndex& cell_mat_index,
			 const ValueIndex& value_index) const {
		return unified_data_get(data, cell_mat_index, value_index);
	}
	
private:
	T data;
};

template<typename T>
class OUT {
public:
        OUT(T p_data) : data(p_data)
	{
		// TODO: Solve this.
		// constexpr std::size_t N = T::N;
		// using dtype = typename T::dtype;
		
		// static_assert(std::is_same<T, CellData<N, dtype>>::value
		// 	      || std::is_same<T, MatData<dtype>>::value
		// 	      || std::is_same<T, CellMatData<N, dtype>>::value);
	}
	
	typename T::dtype& get(const CellMatIndex& cell_mat_index,
			       const ValueIndex& value_index) {
		return unified_data_get(data, cell_mat_index, value_index);
	}
	
private:
	T data;
};

template <typename dtype = double>
class ReduceProxy {
public:
	ReduceProxy(dtype& p_reduced_value,
		    std::function<dtype(dtype, dtype)>& p_reducer)
		: reduced_value(p_reduced_value),
		  reducer(p_reducer)
		{
		}
		
	void operator<<(dtype value) {
		reduced_value = reducer(reduced_value, value);
	}
		
private:
	dtype& reduced_value;
	std::function<dtype(dtype, dtype)>& reducer;
};

template<typename T>
class REDUCE {
public:
	using dtype = typename T::dtype;
	
	REDUCE(std::function<dtype(dtype, dtype)> p_reducer, T p_data)
		: reducer(p_reducer),
		  data(p_data)
	{
		// TODO: Solve this.
		// constexpr std::size_t N = T::N;
		
		// static_assert(std::is_same<T, CellData<N, dtype>>::value
		// 	      || std::is_same<T, MatData<dtype>>::value
		// 	      || std::is_same<T, CellMatData<N, dtype>>::value);
	}

	ReduceProxy<typename T::dtype> get(const CellMatIndex& cell_mat_index,
					   const ValueIndex& value_index) {
		dtype& value
			= unified_data_get(data, cell_mat_index, value_index);
		return ReduceProxy<dtype>(value, reducer);
	}

private:
	// A commutative and associative function.
	std::function<dtype(dtype, dtype)> reducer;
	T data;
};

// template <std::size_t N>
// class Stencil {
// public:
// 	Stencil(std::vector<Offsets<N>> p_offsets)
// 		: offsets(p_offsets)
// 	{
// 	}

// 	bool contains_offset(const Offsets<N>& offset) const {
// 		auto it = std::find(offsets.begin(), offsets.end(), offset);
// 		return it != offsets.end();
// 	}
	
// private:
// 	std::vector<Offsets<N>> offsets;
// };

// template<typename T>
// class NeighProxy {
// public:
// 	using dtype = typename T::dtype;

// 	NeighProxy(const Coords<T::N> p_cell_coords,
// 		   const std::size_t p_mat_index,
// 		   T p_data, const Stencil<T::N> p_stencil)
// 		: cell_coords(p_cell_coords), mat_index(p_mat_index),
// 		  data(p_data), stencil(p_stencil)
// 	{
// 	}

// 	const dtype& get_neigh(const Offsets<T::N>& offset) const {
// 		if (!stencil.contains_offset(offset)) {
// 			throw "No such offset in stencil.";
// 		}
		
// 		const Coords<T::N> neighbour_coords = cell_coords + offset;
// 		return unified_data_get(data, neighbour_coords, mat_index);
// 	}

// 	const dtype& operator[](const Offsets<T::N>& offset) const {
// 		return get_neigh(offset);
// 	}
	
// private:
// 	const Coords<T::N> cell_coords;
// 	const std::size_t mat_index;

// 	T data;
// 	const Stencil<T::N> stencil;	
// };

// template<typename T>
// class NEIGH {
// public:	
// 	using dtype = typename T::dtype;
	
// 	NEIGH(T p_data, const Stencil<T::N> p_stencil)
// 		: data(p_data), stencil(p_stencil)
// 	{
// 		constexpr std::size_t N = T::N;
		
// 		static_assert(std::is_same<T, CellData<N, dtype>>::value
// 			      || std::is_same<T, MatData<dtype>>::value
// 			      || std::is_same<T, CellMatData<N, dtype>>::value);
// 	}

// 	NeighProxy<T>
// 	get(const Coords<T::N>& cell_index, const std::size_t mat_index) {
// 		return NeighProxy<T>(cell_index,
// 				     mat_index,
// 				     data,
// 				     stencil);
// 	}

// private:
// 	T data;
// 	const Stencil<T::N> stencil;
// };

// template <std::size_t N>
// class INDEX {
// public:
// 	INDEX() {};

// 	Coords<N>
// 	get(const Coords<N>& cell_index, const std::size_t mat_index) {
// 		return cell_index;
// 	}
// };


// template<typename dtype = double>
// class FREE_SCALAR {
// public:
// 	FREE_SCALAR(const dtype p_value)
// 		: value(p_value)
// 	{
// 	}

// 	template<std::size_t N>
// 	const dtype& get(const Coords<N>& cell_index,
// 			 const std::size_t mat_index) const {
// 		return value;
// 	}
// private:
// 	const dtype value;
// };

// template<typename dtype = double>
// class FREE_ARRAY {
// public:
// 	using array_type = std::vector<dtype>;
	
// 	FREE_ARRAY(const array_type& p_values)
// 		: values(p_values)
// 	{
// 	}

// 	template<std::size_t N>
// 	const array_type& get(const Coords<N>& cell_index,
// 			      const std::size_t mat_index) const {
// 		return values;
// 	}
	
// private:
// 	const array_type& values;
// };

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_ARGUMENTS_HPP
