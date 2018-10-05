SET(SRC_FILES
        IndexGenerator.cpp
   )

include_directories("full_matrix")
include_directories("compressed_cell_centric")
include_directories(".")

add_executable(MM main.cpp)
