#include "Coords.hpp"

#include "gtest/gtest.h"

namespace MM {

namespace {

TEST(CoordsTest, ConvertTOFlatIndex) {
  const Coords<2> coords{2u, 3u};
  const std::array<std::size_t, 2> size{4u, 5u};
  const std::size_t result = coords_to_flat_index(coords, size);

  ASSERT_EQ(3*4 + 2, result);
}

TEST(CoordsTest, ConvertFromFlatIndex) {
  const std::size_t flat_index = 14;
  const std::array<std::size_t, 2> size{4u, 5u};

  const Coords<2> expected{2u, 3u};
  const Coords<2> result = flat_index_to_coords(flat_index, size);

  ASSERT_EQ(expected, result);
}

TEST(CoordsTest, CoordsToFlatRoundtrip) {
  constexpr std::array<std::size_t, 3> size{2, 3, 4};
  for(std::size_t i = 0; i < size[0]; i++) {
    for (std::size_t j = 0; j < size[1]; j++) {
      for (std::size_t k = 0; k < size[2]; k++) {
        const Coords<3> coords(i, j, k);
        const std::size_t flat_index = coords_to_flat_index(coords, size);
        const Coords<3> result = flat_index_to_coords(flat_index, size);
        EXPECT_EQ(coords, result);
      }
    }
  }
}

} // anonymous namespace

} // namespace MM
