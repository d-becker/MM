#ifndef MM_COMPUTATION_HPP
#define MM_COMPUTATION_HPP

#include "Data.hpp"
#include "IndexGenerator.hpp"

namespace MM {

template <std::size_t N, typename dtype = double>
class Computation {
public:
	Computation(Data<N, dtype>& p_data,
		    IndexGenerator<N> p_index_generator)
		: data(p_data),
		  index_generator(p_index_generator)
	{
	}
		
	template<typename FUNC, typename ...ARGS>
	void compute(FUNC func, ARGS ...args) {
		// TODO: check all args belong to this->data.
		
		while (index_generator.has_next()) {
			const Coords<N> coords = index_generator.next();
			for (std::size_t mat_index = 0;
			     mat_index < data.get_mat_number();
			     ++mat_index) {
				func(args.get(coords, mat_index)...);
			}
		}
	
	}

private:
	Data<N, dtype>& data;
	IndexGenerator<N> index_generator;
};

} // namespace MM

#endif // MM_COMPUTATION_HPP
