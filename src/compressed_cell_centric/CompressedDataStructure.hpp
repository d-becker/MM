#ifndef MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP
#define MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

#include "mm_assert.hpp"

#include "Coords.hpp"

namespace MM {

namespace compressed_cell_centric {

struct Cell {
	// The number of materials in the cell.
	std::size_t nmats;

	// If nmats is 1, this is the material number,
	// otherwise an index into mixed_storage.
	std::size_t imat;

  // TODO:
  // Only for measurements.
  std::size_t a, b;
};

struct MixedStorageCell {
#ifdef MM_LINKED
	// The index of the next material in the same cell
	// or -1 if this is the last one.
	std::size_t nextfrac;
#endif

	// The index of the cell this material belongs to.
	std::size_t frac2cell;

  // The number of the material.
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

	class CellIterator {
	public:
		CellIterator(const CompressedDataStructure& p_structure,
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

		CellIterator& operator++() {
#ifdef MM_LINKED
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

			mixed_storage_index = -1;
			return *this;
#else
      const Cell& cell = structure.cell_at(cell_index);

      if (cell.nmats > 1) {
        mixed_storage_index++;

        const std::size_t end_index = cell.imat + cell.nmats;
        if (mixed_storage_index < end_index) {
          return *this;
        }
      }

      mixed_storage_index = -1;
      return *this;
#endif
		}

		bool operator==(const CellIterator& other) const {
			return &structure == &other.structure
				&& cell_index == other.cell_index
				&& mixed_storage_index
				    == other.mixed_storage_index;
		}

		bool operator!=(const CellIterator& other) const {
			return !(*this == other);
		}
	private:
		const CompressedDataStructure& structure;
		const std::size_t cell_index;
		std::size_t mixed_storage_index;
	};

	class CellIteration {
	public:
		CellIteration(const CompressedDataStructure& p_structure,
			      const std::size_t p_cell_index)
			: structure(p_structure),
			  cell_index(p_cell_index)
		{
		}

		CellIterator begin() const {
			if (structure.cell_number() == 0) {
				return CellIterator(structure, 0, -1);
			}

			const Cell& cell = structure.cell_at(cell_index);
			if (cell.nmats <= 1) {
				return CellIterator(structure, cell_index, 0);
			}

			return CellIterator(structure,
					    cell_index,
					    cell.imat);
		}

		CellIterator end() const {
			return CellIterator(structure, cell_index, -1);
		}
	private:
		const CompressedDataStructure& structure;
		const std::size_t cell_index;
	};

	CellIteration cell_iteration(const std::size_t cell_index) const {
		return CellIteration(*this, cell_index);
	}

	const Cell& cell_at(const std::size_t index) const {
    MM_ASSERT(index < structure.size());
    return structure[index];
	}

	const MixedStorageCell& mixed_cell_at(const std::size_t index) const {
    MM_ASSERT(index < mixed_storage.size());
    return mixed_storage[index];
	}

	std::size_t cell_number() const {
		return structure.size();
	}

	std::size_t mixed_storage_size() const {
		return mixed_storage.size();
	}
//private:
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

  void add_materials_to_mixed_storage(std::size_t cell_index,
		const std::vector<std::size_t>& materials_in_cell)
	{
		for (const std::size_t material : materials_in_cell) {
			MixedStorageCell ms;
#ifdef MM_LINKED
			ms.nextfrac = mixed_storage.size() + 1;
#endif
			ms.frac2cell = cell_index;
			ms.material = material;
			mixed_storage.emplace_back(ms);
		}

#ifdef MM_LINKED
		mixed_storage.back().nextfrac = -1;
#endif
	}

	std::vector<Cell> structure;
	std::vector<MixedStorageCell> mixed_storage;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_DATA_STRUCTURE_HPP
