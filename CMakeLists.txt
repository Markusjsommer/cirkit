cmake_minimum_required(VERSION 3.16)
project(cirkit)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(include)
include_directories(include/parallel_hashmap)
include_directories(include/cereal)

add_executable(cirkit main.cpp include/gzip include/tqdm.h include/FastaReader.cpp include/FastaReader.h include/build_index.cpp include/build_index.h include/FileReader.cpp include/FileReader.h include/kmer_flathash.cpp include/kmer_flathash.h include/BitsetManager.cpp include/BitsetManager.h)


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(cirkit ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

install(
        TARGETS cirkit
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
)