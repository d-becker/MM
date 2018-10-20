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

template<typename _dtype = double>
using MatData = full_matrix::MatData<_dtype>;

template<std::size_t _N, typename _dtype = double>
class CellMatData {
public:
	constexpr static std::size_t N = _N;
	using dtype = _dtype;

	CellMatData(const std::array<std::size_t, N>& p_size,
		    const CompressedDataStructure& p_structure,
		    std::vector<dtype> p_cell_values,
		    std::vector<dtype> p_mixed_storage_values)
		: size(p_size),
		  structure(p_structure),
		  cell_values(p_cell_values),
		  mixed_storage_values(p_mixed_storage_values)
	{
		// TODO: Check the vector lengths.
	}

	const std::array<std::size_t, N>& get_size() const {
		return size;
	}

	const dtype& at(const Coords<N> cell_index,
			const std::size_t mat_index) const {
		return const_cast<CellMatData*>(this)
			->at(cell_index, mat_index);
	}

	dtype& at(const Coords<N> cell_index,
		  const std::size_t mat_index) {
	        const std::size_t flat_index
			= multidim_index_to_raw(cell_index, size);
		const Cell& cell = structure.cell_at(flat_index);

		if (cell.nmats == 0) {
		        throw std::out_of_range(
				"Non-existant material requested.");
		}
		
		if (cell.nmats == 1) {
			if (cell.imat == mat_index) {
				return cell_values.at(flat_index);
			}

			throw std::out_of_range(
				"Non-existant material requested.");
		}

		auto iteration = structure.mixed_mat_iteration(flat_index);
		for (const std::size_t mixed_mat_index : iteration) {
			const MixedStorageCell& mixed_cell
				= structure.mixed_cell_at(mixed_mat_index);
			if (mixed_cell.material == mat_index) {
				return mixed_storage_values.at(mixed_mat_index);
			}
		}

		throw std::out_of_range("Non-existant material requested.");
	}

	class MaterialIterator {
	public:
		MaterialIterator(const bool p_is_end,
				 const bool p_is_multimat,
				 const std::size_t p_cell_index,
				 CompressedDataStructure::MixedStorageIterator p_mixed_iterator,
				 CompressedDataStructure::MixedStorageIterator p_mixed_end_iterator,
				 std::vector<dtype>& p_cell_values,
				 std::vector<dtype>& p_mixed_storage_values)
			: is_end(p_is_end),
			  is_multimat(p_is_multimat),
			  cell_index(p_cell_index),
			  mixed_iterator(p_mixed_iterator),
			  mixed_end_iterator(p_mixed_end_iterator),
			  cell_values(p_cell_values),
			  mixed_storage_values(p_mixed_storage_values)
		{
		}

		dtype& operator*() {
			if (is_multimat) {
				const std::size_t flat_index = *mixed_iterator;
				return mixed_storage_values.at(flat_index);
			}	

			return cell_values.at(cell_index);
		}

		MaterialIterator operator++() {
			if (is_multimat) {
				++mixed_iterator;
				if (mixed_iterator == mixed_end_iterator) {
					is_end = true;
				}
			} else {
				is_end = true;
			}

			return *this;
		}

		MaterialIterator operator++(int) {
			MaterialIterator old_value = (*this);
			++(*this);
			return old_value;
		}

		bool operator==(const MaterialIterator& other) const {
			return is_end == other.is_end
				&& is_multimat == other.is_multimat
				&& cell_index == other.cell_index
				&& mixed_iterator == other.mixed_iterator
				&& mixed_end_iterator == other.mixed_end_iterator
				&& (&cell_values) == (&other.cell_values)
				&& (&mixed_storage_values) == (&other.mixed_storage_values);
		}

		bool operator!=(const MaterialIterator& other) const {
			return !(*this == other);
		}

	private:
		bool is_end;
		bool is_multimat;
		std::size_t cell_index;
		CompressedDataStructure::MixedStorageIterator mixed_iterator;
		CompressedDataStructure::MixedStorageIterator mixed_end_iterator;
		std::vector<dtype>& cell_values;
		std::vector<dtype>& mixed_storage_values;
	};

	class MaterialIterationProxy {
	public:
		MaterialIterationProxy(const std::size_t p_cell_index,
				       const CompressedDataStructure& p_structure,
				       std::vector<dtype>& p_cell_values,
				       std::vector<dtype>& p_mixed_storage_values)
			: cell_index(p_cell_index),
			  structure(p_structure),
			  cell_values(p_cell_values),
			  mixed_storage_values(p_mixed_storage_values)
		{
		}

		MaterialIterator begin() {
			const Cell& cell = structure.cell_at(cell_index);
			bool is_multimat = cell.nmats > 1;
		        auto iteration_proxy = structure.mixed_mat_iteration(cell_index);

			return MaterialIterator(
				false,
			        is_multimat,
				cell_index,
				iteration_proxy.begin(),
			        iteration_proxy.end(),
				cell_values,
				mixed_storage_values
				);
		}

		MaterialIterator end() {
			const Cell& cell = structure.cell_at(cell_index);
			bool is_multimat = cell.nmats > 1;
		        auto iteration_proxy = structure.mixed_mat_iteration(cell_index);

			return MaterialIterator(
				true,
			        is_multimat,
				cell_index,
				iteration_proxy.end(),
			        iteration_proxy.end(),
				cell_values,
				mixed_storage_values
				);
		}

	private:
		const std::size_t cell_index;
		const CompressedDataStructure& structure;
		std::vector<dtype>& cell_values;
		std::vector<dtype>& mixed_storage_values;
	};

	MaterialIterationProxy mat_iteration(const Coords<N> cell_index) {
		const std::size_t flat_index = multidim_index_to_raw(cell_index,
								     size);
		return MaterialIterationProxy(flat_index,
					      structure,
					      cell_values,
					      mixed_storage_values);
	}
private:
	const std::array<std::size_t, N>& size;
	const CompressedDataStructure& structure;
	std::vector<dtype> cell_values;
	std::vector<dtype> mixed_storage_values;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATASETS_HPP
