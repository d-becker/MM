SET(SRC_FILES
        compressed_cell_centric/test_CompressedDataStructure.cpp
   )

# include_directories(".")

add_executable(MM_Test ${SRC_FILES})
target_link_libraries(MM_Test gtest_main)
