#ifndef MM_ARGUMENTS_HPP
#define MM_ARGUMENTS_HPP

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>

#include "Coords.hpp"
#include "mm_assert.hpp"

#include "full_matrix/Datasets.hpp"

namespace MM {

namespace full_matrix {

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
dtype& unified_data_get(CellData<N, dtype> data,
                        Coords<N> cell_index,
                        const std::size_t /* mat_index */) {
  fit_index_to_reduced_data(cell_index, data.get_size());
  return data[cell_index];
}

template<std::size_t N, typename dtype>
dtype& unified_data_get(MatData<N, dtype> data,
                        const Coords<N>& cell_index,
                        const std::size_t mat_index) {
  return data[mat_index];
}

template<std::size_t N, typename dtype>
dtype& unified_data_get(CellMatData<N, dtype> data,
                        Coords<N> cell_index,
                        const std::size_t mat_index) {
  fit_index_to_reduced_data(cell_index, data.get_size());
  return data.at(cell_index, mat_index);
}

template<typename T>
class IN {
public:
  using dtype = typename T::dtype;

  IN(T p_data) : data(p_data)
  {
    // constexpr std::size_t N = T::N;

    // static_assert(std::is_same<T, CellData<N, dtype>>::value
    //         || std::is_same<T, MatData<dtype>>::value
    //         || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  const dtype& get(const Coords<T::N>& cell_index,
                   const std::size_t mat_index) {
    return unified_data_get(data, cell_index, mat_index);
  }

  typename T::dtype *get_raw(std::array<size_t, T::N> &shape) {
    return data.get_raw(shape);
  }

private:
  T data;
};

template<typename T>
class OUT {
public:
        OUT(T p_data) : data(p_data)
  {
    // constexpr std::size_t N = T::N;
    // using dtype = typename T::dtype;

    // static_assert(std::is_same<T, CellData<N, dtype>>::value
    //         || std::is_same<T, MatData<dtype>>::value
    //         || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  typename T::dtype& get(const Coords<T::N>& cell_index,
                         const std::size_t mat_index) {
    return unified_data_get(data, cell_index, mat_index);
  }

  typename T::dtype *get_raw(std::array<size_t, T::N> &shape) {
    return data.get_raw(shape);
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
    // constexpr std::size_t N = T::N;

    // static_assert(std::is_same<T, CellData<N, dtype>>::value
    //         || std::is_same<T, MatData<dtype>>::value
    //         || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  ReduceProxy<typename T::dtype> get(const Coords<T::N>& cell_index,
             const std::size_t mat_index) {
    dtype& value = unified_data_get(data, cell_index, mat_index);
    return ReduceProxy<dtype>(value, reducer);
  }

  typename T::dtype * get_raw(std::array<size_t, T::N> &shape) {
    return data.get_raw(shape);
  }

private:
  std::function<dtype(dtype, dtype)> reducer; // A commutative and associative function.
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

template<typename T>
class NeighProxy {
public:
  using dtype = typename T::dtype;

  NeighProxy(const Coords<T::N> p_cell_coords,
             const std::size_t p_mat_index,
             T p_data,
             const Stencil<T::N> p_stencil)
    : cell_coords(p_cell_coords), mat_index(p_mat_index),
      data(p_data), stencil(p_stencil)
  {
  }

  bool has_neigh(const Offsets<T::N>& offset) const {
    return true;
  }

  const dtype& get_neigh(const Offsets<T::N>& offset) const {
    MM_ASSERT(stencil.contains_offset(offset));

    const Coords<T::N> neighbour_coords = cell_coords + offset;
    return unified_data_get(data, neighbour_coords, mat_index);
  }

  const dtype& operator[](const Offsets<T::N>& offset) const {
    return get_neigh(offset);
  }

  struct Token {
    friend class NeighProxy;

    static_assert(std::is_same<T, CellMatData<T::N, dtype>>::value);
    bool is_valid() const {
      for (std::size_t i = 0; i < T::N; i++) {
        if (coords[i] == -1) {
          return false;
        }
      }

      return true;
    }

  private:
    static Token invalid_token() {
      Coords<T::N> res;
      res[0] = -1;
      return res;
    }

    Token(Coords<T::N> p_coords) : coords(p_coords) {}

    Coords<T::N> coords;
  };

  Token get_cell_mat_token(const Offsets<T::N>& neighbour_offset) const {
    const Coords<T::N> neighbour_coords = cell_coords + neighbour_offset;
    return Token(neighbour_coords);
  }

  const dtype get_with_token(const Token& token) const {
    MM_ASSERT(token.is_valid());

    return unified_data_get(data, token.coords, mat_index);
  }

private:
  const Coords<T::N> cell_coords;
  const std::size_t mat_index;

  T data;
  const Stencil<T::N> stencil;
};

template <class T>
class NeighProxyDirect {
  using dtype = typename T::dtype;
  public:
  NeighProxyDirect(std::array<std::size_t, T::N> _shape, const typename T::dtype * __restrict__ _ptr) : shape(_shape), ptr(_ptr) {}

  const typename T::dtype& operator[](const std::array<int, T::N> offsets) const {
   return ptr[offsets[0] + offsets[1]*shape[0]];
  }

  // TODO: Decide if we would like to use these.
  struct Token {
    friend class NeighProxyDirect;

    static_assert(std::is_same<T, CellMatData<T::N, dtype>>::value);

    bool is_valid() const {
      return true;
    }

  private:
    Token(Offsets<T::N> p_offsets) : offsets(p_offsets) {}

    Offsets<T::N> offsets;
  };

  Token get_cell_mat_token(const Offsets<T::N>& neighbour_offset) const {
    return Token(neighbour_offset);
  }

  const dtype get_with_token(const Token& token) const {
    MM_ASSERT(token.is_valid());


    return ptr[token.offsets[0] + token.offsets[1]*shape[0]];
  }


  private:
    std::array<std::size_t, T::N> shape;
    const typename T::dtype * __restrict__ ptr;
};

template<typename T>
class NEIGH {
public:
  using dtype = typename T::dtype;

  NEIGH(T p_data, const Stencil<T::N> p_stencil)
    : data(p_data), stencil(p_stencil)
  {
    // constexpr std::size_t N = T::N;

    // static_assert(std::is_same<T, CellData<N, dtype>>::value
    //         || std::is_same<T, MatData<dtype>>::value
    //         || std::is_same<T, CellMatData<N, dtype>>::value);
  }

  NeighProxy<T>
  get(const Coords<T::N>& cell_index, const std::size_t mat_index) {
    return NeighProxy<T>(cell_index,
             mat_index,
             data,
             stencil);
  }

  typename T::dtype * get_raw(std::array<size_t, T::N> &shape) {
    return data.get_raw(shape);
  }
private:
  T data;
  const Stencil<T::N> stencil;
};

template <std::size_t N>
class INDEX {
public:
  INDEX() {};

  Coords<N>
  get(const Coords<N>& cell_index, const std::size_t /* mat_index */) {
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
  const dtype& get(const Coords<N>& /* cell_index */,
                   const std::size_t /* mat_index */) const {
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
  const array_type& get(const Coords<N>& cell_index,
                        const std::size_t mat_index) const {
    return values;
  }

private:
  const array_type& values;
};

} // namespace full_matrix

} // namespace MM

#endif // MM_ARGUMENTS_HPP
