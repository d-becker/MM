#ifndef MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP
#define MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP

#include <cstddef>
#include <vector>

namespace MM {

namespace compressed_cell_centric {

struct Cell {
	// The number of materials in the cell.
	std::size_t nmats;

	// If nmats is 1, this is the material number,
	// otherwise an index into mixed_storage.
	std::size_t imat;
};

struct MixedStorageCell {
	// The index of the next material in the same cell
	// or -1 if this is the last one.
	std::size_t nextfrac;

	// The index of the cell this material belongs to.
	std::size_t frac2cell;

	std::size_t material;
};

class CompressedDataStructure {
public:
	CompressedDataStructure(
		const std::vector<std::vector<std::size_t>>& materials)
		: structure(materials.size())
	{
		const std::size_t num_of_cells = materials.size();
		for (std::size_t i = 0; i < num_of_cells; ++i) {
			const std::vector<std::size_t>& materials_in_cell
				= materials[i];
			
			handle_cell(i, materials_in_cell);
		}
	}

	class MixedStorageIterator {
	public:
		MixedStorageIterator(const std::vector<MixedStorageCell>& p_mixed_data,
				     const std::size_t p_ms_index)
			: mixed_data(p_mixed_data),
			  ms_index(p_ms_index)
		{
		}

		std::size_t operator*() const {
			return ms_index;
		}

		MixedStorageIterator operator++() {
		        const std::size_t next = mixed_data.at(ms_index).nextfrac;
			ms_index = next;
			return *this;
		}

		bool operator==(const MixedStorageIterator& other) const {
			return &mixed_data == &other.mixed_data
				&& ms_index == other.ms_index;
		}

		bool operator!=(const MixedStorageIterator& other) const {
			return !(*this == other);
		}

	private:
		const std::vector<MixedStorageCell>& mixed_data;
		std::size_t ms_index;
	};

	class MixedStorageIterationProxy {
	public:
		MixedStorageIterationProxy(const CompressedDataStructure& p_data,
					   const std::size_t p_cell_index)
			: data(p_data),
			  cell_index(p_cell_index)
		{
		}

		const MixedStorageIterator begin() const {
			const Cell& cell = data.cell_at(cell_index);
			if (cell.nmats <= 1) {
				return MixedStorageIterator(data.mixed_storage, -1);
			} else {
				return MixedStorageIterator(data.mixed_storage, cell.imat);
			}
		}

		const MixedStorageIterator end() const {
			return MixedStorageIterator(data.mixed_storage, -1);
		}
	private:
		const CompressedDataStructure& data;
		const std::size_t cell_index;
	};

	const Cell& cell_at(const std::size_t index) const {
		return structure.at(index);
	}

	const MixedStorageCell& mixed_cell_at(const std::size_t index) const {
		return mixed_storage.at(index);
	}

	const MixedStorageIterationProxy mixed_mat_iteration(const std::size_t cell_index) const {
		return MixedStorageIterationProxy(*this, cell_index);
	}

	std::size_t mixed_storage_size() const {
		return mixed_storage.size();
	}
private:
	void handle_cell(std::size_t cell_index,
			 const std::vector<std::size_t>& materials_in_cell) {
		structure[cell_index].nmats = materials_in_cell.size();

		if (structure[cell_index].nmats == 1) {
			structure[cell_index].imat = materials_in_cell[0];
		} else if (structure[cell_index].nmats > 1) {
			structure[cell_index].imat = mixed_storage.size();

			add_materials_to_mixed_storage(cell_index,
						       materials_in_cell);
		}
	}

	void add_materials_to_mixed_storage(
		std::size_t cell_index,
		const std::vector<std::size_t>& materials_in_cell)
	{
		for (const std::size_t material : materials_in_cell) {
			MixedStorageCell ms;
			ms.nextfrac = mixed_storage.size() + 1;
			ms.frac2cell = cell_index;
			ms.material = material;
			mixed_storage.emplace_back(ms);
		}

		mixed_storage.back().nextfrac = -1;
	}

	std::vector<Cell> structure;
	std::vector<MixedStorageCell> mixed_storage;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP
