#ifndef MM_COMPRESSED_CELL_ARGUMENTS_HPP
#define MM_COMPRESSED_CELL_ARGUMENTS_HPP

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

#include "Coords.hpp"
#include "compressed_cell_centric/CompressedDataStructure.hpp"
#include "compressed_cell_centric/Datasets.hpp"
#include "compressed_cell_centric/Data.hpp"
#include "mm_assert.hpp"

namespace MM {

namespace compressed_cell_centric {

template<std::size_t N, typename dtype>
dtype& unified_data_get(CellData<N, dtype> data,
                        const CellMatIndex& cell_mat_index,
                        const ValueIndex& /* value_index */) {
  return data.at_raw_index(cell_mat_index.cell_index);
}

template<std::size_t N, typename dtype>
dtype& unified_data_get(MatData<N, dtype> data,
                        const CellMatIndex& cell_mat_index,
                        const ValueIndex& /* _value_index */) {
  return data[cell_mat_index.mat_index];
}

template<std::size_t N, typename dtype>
dtype& unified_data_get(CellMatData<N, dtype> data,
                        const CellMatIndex& /* cell_mat_index */,
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
    constexpr std::size_t N = T::N;

    static_assert(std::is_same<T, CellData<N, dtype>>::value
            || std::is_same<T, MatData<N, dtype>>::value
            || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  template<std::size_t N>
  const dtype& get(const Coords<N>&,
                   const Data<N, dtype>&,
                   const CellMatIndex& cell_mat_index,
                   const ValueIndex& value_index) const {
    return unified_data_get(data, cell_mat_index, value_index);
  }

  typename T::dtype *get_raw() {
    return data.get_raw();
  }

  typename T::dtype *get_raw_list() {
    return data.get_raw_list();
  }

private:
  T data;
};

template<typename T>
class OUT {
public:
  using dtype = typename T::dtype;

  OUT(T p_data) : data(p_data)
  {
    constexpr std::size_t N = T::N;
    using dtype = typename T::dtype;

    static_assert(std::is_same<T, CellData<N, dtype>>::value
            || std::is_same<T, MatData<N, dtype>>::value
            || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  template <std::size_t N>
  typename T::dtype& get(const Coords<N>&,
                         const Data<N, dtype>&,
                         const CellMatIndex& cell_mat_index,
                         const ValueIndex& value_index) {
    return unified_data_get(data, cell_mat_index, value_index);
  }

  typename T::dtype *get_raw() {
    return data.get_raw();
  }

  typename T::dtype *get_raw_list() {
    return data.get_raw_list();
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
    constexpr std::size_t N = T::N;

    static_assert(std::is_same<T, CellData<N, dtype>>::value
            || std::is_same<T, MatData<N, dtype>>::value
            || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  template<std::size_t N>
  ReduceProxy<typename T::dtype> get(const Coords<N>&,
                                     const Data<N, dtype>&,
                                     const CellMatIndex& cell_mat_index,
                                     const ValueIndex& value_index) {
    dtype& value = unified_data_get(data, cell_mat_index, value_index);
    return ReduceProxy<dtype>(value, reducer);
  }

  // A commutative and associative function.
  std::function<dtype(dtype, dtype)> reducer;
  T data;
private:
};

template <std::size_t N>
class Stencil {
public:
  Stencil(std::vector<Offsets<N>> p_offsets) noexcept
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

template<typename T>
class NeighProxyDirect {
public:
  using dtype = typename T::dtype;

  NeighProxyDirect(const Data<T::N, dtype>& p_data,
             const Coords<T::N> p_cell_coords,
             const CellMatIndex& p_cell_mat_index,
             const ValueIndex& p_value_index,
             T p_dataset)
    : data(p_data),
      cell_coords(p_cell_coords),
      cell_mat_index(p_cell_mat_index),
      value_index(p_value_index),
      dataset(p_dataset)
  {
  }

  bool has_neigh(const Offsets<T::N>& offset) const {
    const Coords<T::N> neighbour_coords = cell_coords + offset;
    if (std::is_same<T, CellData<T::N, dtype>>::value) {
      return true;
    } else {
      return has_neigh_cell_mat_data(neighbour_coords);
    }
  }

  const dtype& get_neigh(const Offsets<T::N>& offset) const {
    MM_ASSERT(has_neigh(offset));

    const Coords<T::N> neighbour_coords = cell_coords + offset;

    if (std::is_same<T, CellData<T::N, dtype>>::value) {
      return get_neigh_cell_data(neighbour_coords);
    } else {
      return get_neigh_cell_mat_data(neighbour_coords);
    }
  }

  const dtype& operator[](const Offsets<T::N>& offset) const {
    return get_neigh(offset);
  }

  struct Token {
    friend class NeighProxyDirect;
    static_assert(std::is_same<T, CellMatData<T::N, dtype>>::value);


    // Returns whether the token is valid. An invalid token can be returned if
    // the user asks for a token for a material in a cell where the material is
    // not present.
    bool is_valid() const {
      return value_index.index != -1;
    }

  private:
    static Token invalid_token() {
      return Token(ValueIndex(ValueIndex::Type::SINGLE_MAT, -1));
    }

    Token(ValueIndex p_value_index) : value_index(p_value_index) {}

    ValueIndex value_index;
  };

  Token get_cell_mat_token(const Offsets<T::N>& neighbour_offset) const {
    static_assert(std::is_same<T, CellMatData<T::N, dtype>>::value);

    const Coords<T::N> neighbour_coords = cell_coords + neighbour_offset;

    for (const std::pair<CellMatIndex, ValueIndex> pair
           : data.cell_iteration(neighbour_coords)) {
      const CellMatIndex& n_cell_mat_index = pair.first;
      if (n_cell_mat_index.mat_index == cell_mat_index.mat_index) {
        const ValueIndex& n_value_index = pair.second;
        return Token(n_value_index);
      }
    }

    // Material not found, returning an invalid token.
    return Token::invalid_token();
  }

  const dtype& get_with_token(const Token& token) const {
    MM_ASSERT(token.is_valid());

    // CellMatIndex is not used.
    return unified_data_get(dataset, CellMatIndex(-1, -1), token.value_index);
  }

private:
  const dtype&
  get_neigh_cell_data(const Coords<T::N>& neighbour_coords) const {
    const std::size_t flat_index = multidim_index_to_raw(
      neighbour_coords,
      data.get_size());
    const CellMatIndex n_cell_mat_index(flat_index, -1);
    const ValueIndex n_value_index(ValueIndex::Type::SINGLE_MAT, -1);
    return unified_data_get(dataset,
                            n_cell_mat_index,
                            n_value_index);
  }

  const dtype&
  get_neigh_cell_mat_data(const Coords<T::N>& neighbour_coords) const {
    for (const std::pair<CellMatIndex, ValueIndex> pair
           : data.cell_iteration(neighbour_coords)) {
      const CellMatIndex& n_cell_mat_index = pair.first;
      if (n_cell_mat_index.mat_index == cell_mat_index.mat_index) {
        const ValueIndex& n_value_index = pair.second;
        return unified_data_get(dataset,
                                n_cell_mat_index,
                                n_value_index);
      }
    }

    MM_ASSERT(false);
  }

  bool has_neigh_cell_mat_data(const Coords<T::N>& neighbour_coords) const {
    for (const std::pair<CellMatIndex, ValueIndex> pair
           : data.cell_iteration(neighbour_coords)) {
      const CellMatIndex& n_cell_mat_index = pair.first;
      if (n_cell_mat_index.mat_index == cell_mat_index.mat_index) {
        return true;
      }
    }

    return false;
  }

  const Data<T::N, dtype>& data;
  const Coords<T::N> cell_coords;
  const CellMatIndex& cell_mat_index;
  const ValueIndex& value_index;

  T dataset;
};

template <typename T>
class NeighProxy : public NeighProxyDirect<T> {
public:
  using dtype = typename T::dtype;

  NeighProxy(const Data<T::N, dtype>& p_data,
             const Coords<T::N> p_cell_coords,
             const CellMatIndex& p_cell_mat_index,
             const ValueIndex& p_value_index,
             T p_dataset,
             const Stencil<T::N> p_stencil)
    : NeighProxyDirect<T>(p_data,
                 p_cell_coords,
                 p_cell_mat_index,
                 p_value_index,
                 p_dataset),
      stencil(p_stencil)
  {
  }

  bool has_neigh(const Offsets<T::N>& offset) const {
    if (!stencil.contains_offset(offset)) {
      return false;
    }

    return NeighProxyDirect<T>::has_neigh(offset);
  }

private:
  const Stencil<T::N> stencil;
};

template<typename T>
class NEIGH {
public:
  using dtype = typename T::dtype;

  NEIGH(T p_dataset, const Stencil<T::N> p_stencil)
    : dataset(p_dataset), stencil(p_stencil)
  {
    constexpr std::size_t N = T::N;

    static_assert(std::is_same<T, CellData<N, dtype>>::value
            || std::is_same<T, MatData<N, dtype>>::value
            || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  NeighProxy<T>
  get(const Coords<T::N>& coords,
      const Data<T::N, dtype>& data,
      const CellMatIndex& cell_mat_index,
      const ValueIndex& value_index) {
    return NeighProxy<T>(
        data,
        coords,
        cell_mat_index,
        value_index,
        dataset,
        stencil);
  }

  dtype* get_raw() {
    return dataset.get_raw();
  }

  dtype* get_raw_list() {
    return dataset.get_raw_list();
  }

// TODO: private.
// private:
  T dataset;
  const Stencil<T::N> stencil;
};

template <std::size_t N>
class INDEX {
public:
  INDEX() {};

  template<typename dtype>
  Coords<N>
  get(const Coords<N>& cell_index,
      const Data<N, dtype>&,
      const CellMatIndex&,
      const ValueIndex&) {
    return cell_index;
  }
};

template<typename dtype = double>
class FREE_SCALAR {
public:
  FREE_SCALAR(const dtype p_value)
    : value(p_value)
  {
  }

  template<std::size_t N>
  const dtype& get(const Coords<N>&,
                   const Data<N, dtype>&,
                   const CellMatIndex& /* cell_mat_index */,
                   const ValueIndex& /* value_index */) const {
    return value;
  }

private:
  const dtype value;
};

template<typename dtype = double>
class FREE_ARRAY {
public:
  using array_type = std::vector<dtype>;

  FREE_ARRAY(const array_type& p_values)
    : values(p_values)
  {
  }

  template<std::size_t N>
  const array_type& get(const Coords<N>&,
                        const Data<N, dtype>&,
                        const CellMatIndex& /* cell_mat_index */,
                        const ValueIndex& /* value_index */) const {
    return values;
  }

private:
  const array_type& values;
};

} // namespace compressed_cell_centric

} // namespace MM

#endif // MM_COMPRESSED_CELL_ARGUMENTS_HPP
