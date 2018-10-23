#ifndef MM_INDEX_GENERATOR_HPP
#define MM_INDEX_GENERATOR_HPP

#include <array>
#include <cstddef>

namespace MM {

template<std::size_t N>
class IndexGenerator {
public:
	IndexGenerator(const std::array<std::size_t, N> p_begin,
		       const std::array<std::size_t, N> p_end)
		: begin(p_begin),
		  end(p_end),
		  current(p_begin),
		  _has_next(true)
	{
		for (std::size_t i = 0; i < N; ++i) {
			if (begin[i] >= end[i]) {
				_has_next = false;
				break;
			}
		}
	}

	const std::array<std::size_t, N>& get_begin() const {
		return begin;
	}
	
	const std::array<std::size_t, N>& get_end() const {
		return end;
	}

	bool has_next() const {
		return _has_next;
	}
	
	std::array<std::size_t, N> next() {
		std::array<std::size_t, N> res = current;
		
		for (std::size_t i = 0; i < N; ++i) {
			++current[i];

			if (current[i] >= end[i]) {
				if (i == N - 1) {
					_has_next = false;
					break;
				}
				
				current[i] = begin[i];
			} else {
				break;
			}
		}

		return res;
	}
	
private:
	const std::array<std::size_t, N> begin;
	const std::array<std::size_t, N> end; // Exclusive.

	std::array<std::size_t, N> current;
	bool _has_next;
};

} // namespace MM

#endif // MM_INDEX_GENERATOR_HPP
