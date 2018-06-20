#ifndef MM_DIMENSION_ERROR_HPP
#define MM_DIMENSION_ERROR_HPP

#include <stdexcept>
#include <string>

namespace MM {

class DimensionError : public std::runtime_error {
public:
	DimensionError(const std::string& what)
		: std::runtime_error(what)
	{
	}
};

} // namespace MM

#endif // MM_DIMENSION_ERROR_HPP
