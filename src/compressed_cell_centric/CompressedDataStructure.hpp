#ifndef MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP
#define MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP

#include <cstddef>
#include <utility>
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

struct CellMatIndex {
	CellMatIndex(const std::size_t p_cell_index,
		     const std::size_t p_mat_index)
		: cell_index(p_cell_index),
		  mat_index(p_mat_index)
	{
	}
	
	std::size_t cell_index;
	std::size_t mat_index;
};

struct ValueIndex {
	enum struct Type {
		SINGLE_MAT,
		MULTIMAT
	};

	ValueIndex(const Type p_type, const std::size_t p_index)
		: type(p_type),
		  index(p_index)
	{
	}

	Type type;
	std::size_t index;
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

	class Iterator {
	public:
		Iterator(const CompressedDataStructure& p_structure,
			 const std::size_t p_cell_index,
			 const std::size_t p_mixed_storage_index)
			: structure(p_structure),
			  cell_index(p_cell_index),
			  mixed_storage_index(p_mixed_storage_index)
		{
		}
		
		std::pair<CellMatIndex, ValueIndex> operator*() {
			const Cell& cell = structure.cell_at(cell_index);
			if (cell.nmats > 1) {
				const MixedStorageCell& mixed_cell
					= structure.mixed_cell_at(
						mixed_storage_index);
				const std::size_t mat_id = mixed_cell.material;

				const CellMatIndex cell_mat_index(cell_index,
								  mat_id);
				const ValueIndex value_index(
					ValueIndex::Type::MULTIMAT,
					mixed_storage_index);
				return std::make_pair(cell_mat_index,
						      value_index);
			}

			const CellMatIndex cell_mat_index(cell_index,
							  cell.imat);
			const ValueIndex value_index(
				ValueIndex::Type::SINGLE_MAT, cell_index);
			return std::make_pair(cell_mat_index, value_index);
		}
		
		Iterator& operator++() {
			const Cell& cell = structure.cell_at(cell_index);

			if (cell.nmats > 1) {
				const MixedStorageCell& mixed_cell
					= structure.mixed_cell_at(
						mixed_storage_index);
				const std::size_t next = mixed_cell.nextfrac;

				if (next != -1u) {
					mixed_storage_index = next;
					return *this;
				}
			}

			// Either this is a single material cell or we
			// have gone past the last material in the cell.
			++cell_index;

			if (cell_index < structure.cell_number()) {
				const Cell& new_cell
					= structure.cell_at(cell_index);
				mixed_storage_index = new_cell.imat;
			} else {
				mixed_storage_index = -1;
			}
			
			return *this;
		}

		bool operator==(const Iterator& other) const {
			return &structure == &other.structure
				&& cell_index == other.cell_index
				&& mixed_storage_index
				    == other.mixed_storage_index;
		}

		bool operator!=(const Iterator& other) const {
			return !(*this == other);
		}
	private:
		const CompressedDataStructure& structure;
		std::size_t cell_index;
		std::size_t mixed_storage_index;
	};

		Iterator begin() const {
		if (cell_number() == 0) {
			return Iterator(*this, 0, -1);
		}
		
		const Cell& first_cell = cell_at(0);
		if (first_cell.nmats <= 1) {
			return Iterator(*this, 0, -1);
		}
		
		return Iterator(*this, 0, first_cell.imat);
	}

	Iterator end() const {
		return Iterator(*this, cell_number(), -1);
	}

	const Cell& cell_at(const std::size_t index) const {
		return structure.at(index);
	}

	const MixedStorageCell& mixed_cell_at(const std::size_t index) const {
		return mixed_storage.at(index);
	}

	std::size_t cell_number() const {
		return structure.size();
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
