SET(SRC_FILES
        test_Coords.cpp
        compressed_cell_centric/test_Arguments.cpp
        compressed_cell_centric/test_Computation.cpp
        compressed_cell_centric/test_CompressedDataStructure.cpp
        compressed_cell_centric/test_Data.cpp
        compressed_cell_centric/test_Datasets.cpp
   )

# include_directories(".")

add_executable(MM_Test ${SRC_FILES})
target_link_libraries(MM_Test gtest_main)
