#ifndef MM_ASSERT_HPP
#define MM_ASSERT_HPP

#include <cassert>

// Wrapping the assert macro to avoid unused variable warnings in relase
// builds.
#ifdef NDEBUG
#define MM_ASSERT(condition) ((void) (condition))
#else
#define MM_ASSERT(condition) assert(condition)
#endif

#endif // MM_ASSERT_HPP
