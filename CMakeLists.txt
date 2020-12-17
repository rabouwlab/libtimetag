cmake_minimum_required(VERSION 3.9)

project(libtimetag 
		VERSION 0.8
		DESCRIPTION "Library for storage and processing of time-correlated single-photon counting (TCSPC) data."
)

set(CMAKE_CXX_STANDARD 11)

include(GNUInstallDirs)

add_library(libtimetag SHARED)

target_sources(libtimetag
	PRIVATE src/getline.cpp
	PRIVATE src/sstt_file.cpp
	PRIVATE src/algos.cpp
	PRIVATE src/sstt_file2.cpp
)

set_target_properties(libtimetag PROPERTIES PUBLIC_HEADER "include/algos.h;include/sstt_file.h;include/sstt_file2.h")

add_compile_definitions(BUILDING_LIBTIMETAG)

target_include_directories(libtimetag
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
	PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

install(TARGETS libtimetag
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libtimetag)