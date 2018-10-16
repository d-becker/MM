#ifndef MM_TEST_COMPRESSED_UTIL_HPP
#define MM_TEST_COMPRESSED_UTIL_HPP

#include <vector>

namespace MM::compressed_cell_centric {

namespace {

std::vector<std::vector<std::size_t>> get_raw_data(const std::size_t N, const std::size_t M) {
	srand(0); // Deterministic pseudo-random number generation.
		  // TODO: It is not guaranteed that the given number of materials will actually be generated as it is "random".

	std::vector<std::vector<std::size_t>> res(N);

	std::vector<std::size_t> first_cell{0};
	res[0] = first_cell;

	for (std::size_t i = 1; i < N; ++i) {
		const std::size_t length = rand() % M;

		res[i] = std::vector<std::size_t>(length);

		for (std::size_t& material : res[i]) {
			material = rand() % M; // TODO: Avoid duplicate materials.
		}
	}

	return res;
}

} // anonymous namespace

} // namespace MM::compressed_cell_centric

#endif // MM_TEST_COMPRESSED_UTIL_HPP
