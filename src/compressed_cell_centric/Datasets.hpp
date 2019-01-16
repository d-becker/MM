#ifndef MM_COMPRESSED_CELL_DATASETS_HPP
#define MM_COMPRESSED_CELL_DATASETS_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "MultidimArray.hpp"
#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "full_matrix/Datasets.hpp"

namespace MM {

namespace compressed_cell_centric {

template<std::size_t _N, typename _dtype = double>
using CellData = full_matrix::CellData<_N, _dtype>;

template<std::size_t N, typename _dtype = double>
using MatData = full_matrix::MatData<N, _dtype>;

template<std::size_t _N, typename _dtype = double>
class CellMatData {
public:
	constexpr static std::size_t N = _N;
	using dtype = _dtype;

	CellMatData(std::vector<dtype>& p_cell_values,
		    std::vector<dtype>& p_mixed_storage_values)
		: cell_values(p_cell_values),
		  mixed_storage_values(p_mixed_storage_values)
	{
	}
		
	dtype& cell_value_at(const std::size_t index) {
		return cell_values.at(index);
	}

	const dtype& cell_value_at(const std::size_t index) const {
		return cell_values.at(index);
	}

	dtype& mixed_storage_value_at(const std::size_t index) {
		return mixed_storage_values.at(index);
	}

	const dtype& mixed_storage_value_at(const std::size_t index) const {
		return mixed_storage_values.at(index);
	}

	std::size_t cell_number() const {
		return cell_values.size();
	}

	std::size_t mixed_storage_size() const {
		return mixed_storage_values.size();
	}

	dtype* get_raw() {
    	return &cell_values[0];
  	}

  	dtype* get_raw_list() {
    	return &mixed_storage_values[0];
  	}
private:
	// const std::array<std::size_t, N>& size;
	// const CompressedDataStructure& structure;
	std::vector<dtype>& cell_values;
	std::vector<dtype>& mixed_storage_values;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATASETS_HPP
