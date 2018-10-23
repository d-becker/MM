#include <array>
#include <cstddef>
#include <numeric>
#include <vector>

#include "gtest/gtest.h"

#include "compressed_cell_centric/Computation.hpp"

#include "compressed_cell_centric/Arguments.hpp"

namespace MM::compressed_cell_centric {

namespace {

class CompressedComputation : public ::testing::Test {
protected:
	CompressedComputation()
		: size({2, 2}),
		  raw_materials {
		        std::vector<std::size_t>{0},
			std::vector<std::size_t>{1, 3},
			std::vector<std::size_t>{2},
			std::vector<std::size_t>{3},
	          },
		  data(size, raw_materials)
	{
	}

	void fill_density(CellMatData<2> density) {
		density.cell_value_at(0) = 2.3;
		density.cell_value_at(1) = 0.0; // Ignored as multimat cell.
		density.cell_value_at(2) = 4.8;
		density.cell_value_at(3) = 3.6;

		density.mixed_storage_value_at(0) = 0.3;
		density.mixed_storage_value_at(1) = 0.7;
	}

	void fill_volume(CellMatData<2> volume) {
		volume.cell_value_at(0) = 12.3;
		volume.cell_value_at(1) = 0.0; // Ignored as multimat cell.
		volume.cell_value_at(2) = 14.8;
		volume.cell_value_at(3) = 13.6;

		volume.mixed_storage_value_at(0) = 10.3;
		volume.mixed_storage_value_at(1) = 10.7;
	}

	void check_mass(CellMatData<2> mass,
			CellMatData<2> density,
			CellMatData<2> volume) {
		for (std::size_t cell_index = 0;
		     cell_index < mass.cell_number();
		     ++cell_index) {
			const double mass_value
				= mass.cell_value_at(cell_index);
			const double density_value
				= density.cell_value_at(cell_index);
			const double volume_value
				= volume.cell_value_at(cell_index);

			ASSERT_EQ(mass_value, density_value * volume_value);
		}

		for (std::size_t mixed_index = 0;
		     mixed_index < mass.mixed_storage_size();
		     ++mixed_index) {
			const double mass_value
				= mass.mixed_storage_value_at(mixed_index);
			const double density_value
				= density.mixed_storage_value_at(mixed_index);
			const double volume_value
				= volume.mixed_storage_value_at(mixed_index);
			
			ASSERT_EQ(mass_value, density_value * volume_value);
		}
	}

	void check_mass_by_cell(CellData<2> mass,
				CellMatData<2> density,
				CellMatData<2> volume) {
		std::vector<double> cell_products;
		for (std::size_t i = 0; i < density.cell_number(); ++i) {
			cell_products.push_back(
				density.cell_value_at(i)
				* volume.cell_value_at(i));
		}

		std::vector<double> mixed_products;
		for (std::size_t i = 0; i < density.mixed_storage_size(); ++i) {
			mixed_products.push_back(
				density.mixed_storage_value_at(i)
				* volume.mixed_storage_value_at(i));
		}

		double mixed_sum = std::accumulate(mixed_products.begin(),
						   mixed_products.end(),
						   0.0);

		// This is the multimaterial cell.
		cell_products.at(1) = mixed_sum;

		for (std::size_t i = 0; i < cell_products.size(); ++i) {
			ASSERT_EQ(mass.at_raw_index(i), cell_products.at(i));
		}
	}

	void check_multiplication(const CellMatData<2> result,
				  const CellMatData<2> original,
				  const double factor) {
		for (std::size_t i = 0; i < result.cell_number(); ++i) {
			ASSERT_EQ(result.cell_value_at(i),
				  original.cell_value_at(i) * factor);
		}

		for (std::size_t i = 0; i < result.mixed_storage_size(); ++i) {
			ASSERT_EQ(result.mixed_storage_value_at(i),
				  original.mixed_storage_value_at(i) * factor);
		}
	}

	const std::array<std::size_t, 2> size;
	const std::vector<std::vector<std::size_t>> raw_materials;
	Data<2> data;
};

TEST_F(CompressedComputation, in_out) {
	CellMatData<2> density = data.new_cell_mat_data();
	CellMatData<2> volume = data.new_cell_mat_data();

	fill_density(density);
	fill_volume(volume);

	CellMatData<2> mass = data.new_cell_mat_data();

	auto kernel = [] (double density,
			  double volume,
			  double& mass) {
		mass = density * volume;
	};

	IndexGenerator<2> index_generator({0, 0}, {2, 2});
	Computation<2> computation(data, index_generator);
	
        computation.compute(kernel,
			    IN<CellMatData<2>>(density),
			    IN<CellMatData<2>>(volume),
			    OUT<CellMatData<2>>(mass));

	check_mass(mass, density, volume);
}

TEST_F(CompressedComputation, in_reduce) {
	CellMatData<2> density = data.new_cell_mat_data();
	CellMatData<2> volume = data.new_cell_mat_data();

	fill_density(density);
	fill_volume(volume);

	CellData<2> mass_by_cell = data.new_cell_data();

	auto sum = [] (double left, double right) {
			   return left + right;
		   };
	
	auto kernel = [] (double density,
			  double volume,
			  ReduceProxy<double> mass_by_cell) {
			      mass_by_cell << density * volume;
		      };
	
	IndexGenerator<2> index_generator({0, 0}, {2, 2});
	Computation<2> computation(data, index_generator);
	computation.compute(kernel,
			    IN<CellMatData<2>>(density),
			    IN<CellMatData<2>>(volume),
			    REDUCE<CellData<2>>(sum, mass_by_cell));

	check_mass_by_cell(mass_by_cell, density, volume);
}

TEST_F(CompressedComputation, free_scalar) {
	CellMatData<2> density = data.new_cell_mat_data();

	fill_density(density);

	CellMatData<2> double_density = data.new_cell_mat_data();
	const double factor = 2.0;

	auto kernel = [] (double density,
			  const double& free_scalar,
			  double& double_density) {
                double_density = density * free_scalar;
	};

	IndexGenerator<2> index_generator({0, 0}, {2, 2});
	Computation<2> computation(data, index_generator);
	
        computation.compute(kernel,
			    IN<CellMatData<2>>(density),
			    FREE_SCALAR<>(factor),
			    OUT<CellMatData<2>>(double_density));

        check_multiplication(double_density, density, factor);
}


} // anonymous namespace

} // MM::compressed_cell_centric
