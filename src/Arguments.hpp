#ifndef MM_ARGUMENTS_HPP
#define MM_ARGUMENTS_HPP

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

/*
namespace MM {

class Stencil {
public:
	template<typename ...Offsets>
	Stencil(Offsets ...offsets) : offsets({offsets...}) {}
		
	bool contains_offset(const Offset& offset) const {
		return offsets.count(offset) > 0;
	}
	
private:
	std::unordered_set<Offset> offsets;
};

template <typename T>
class NeighProxy {
public:
	NeighProxy(const CELL_ID p_cell, const MAT_ID p_mat,
		   T p_data, const Stencil p_stencil)
		: cell(p_cell), mat(p_mat),
		  data(p_data), stencil(p_stencil)
		{}

	double get_neigh(int x, int y) {
		const Offset offset(x, y);
		if (!stencil.contains_offset(offset)) {
			throw "No offset.";
		}
		
		const std::size_t cell_index = cell.id;

		const std::size_t COLS = data.get_data_structure().get_n_of_columns();
		const std::size_t ROWS = data.get_data_structure().get_n_of_rows();
		
		const std::size_t cell_col = cell_index % COLS;
		const std::size_t cell_row = cell_index / COLS;

		const std::size_t new_x = cell_col + x;
		const std::size_t new_y = cell_row + y;

		// Bounds checking.
		if (new_x >= COLS) {
			std::stringstream s;
			s << "The column index is out of bounds. "
			  << "Provided index: "  << new_x
			  << ", valid range: [0 - " << COLS << ").";
			throw std::out_of_range(s.str());
		}

		if (new_y >= ROWS) {
			std::stringstream s;
			s << "The row index is out of bounds. "
			  << "Provided index: "  << new_y
			  << ", valid range: [0 - " << ROWS << ").";
			throw std::out_of_range(s.str());
		}

		// TODO: check stencil.
		
		const std::size_t new_index = new_y * COLS + new_x;
		const CELL_ID new_id {new_index};
		return unified_data_get(data, new_id, mat);
	}
	
private:
	const CELL_ID cell;
	const MAT_ID mat;

	T data;
	const Stencil stencil;
};

template <typename T>
class NEIGH {
public:
	NEIGH(T p_data, const Stencil stencil)
		: data(p_data)
		{}

	NeighProxy<T> get(const CELL_ID& cell_id, const MAT_ID& mat_id) {
		return NeighProxy<T>(cell_id, mat_id, data, stencil);
	}

private:
	const Stencil stencil;
	T data;
};

*/

} // namespace MM

#endif // MM_ARGUMENTS_HPP
