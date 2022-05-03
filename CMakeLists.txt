cmake_minimum_required(VERSION 3.10)

############################################################
# set the project name
############################################################

project(ChIMES_LSQ
	VERSION 1.0.0
	DESCRIPTION "Utilities to generate ChIMES models"
	LANGUAGES CXX)

# Set default build type and build flags

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "") # Can be updated if we want to adde extra flags 
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_VERBOSE_MAKEFILE on) # Uncomment if you want to see compiler flags
    

# set default preprocessor flag values
 
set(DEBUG_CHEBY "0" CACHE STRING "Debug verbosity level")
set(VERBOSITY   "1" CACHE STRING "General verbosity level")
set(USE_MPI     "1" CACHE STRING "Use MPI")

# Set preprocessor flags to user-specified values, if available

add_definitions(-DDEBUG=${DEBUG})
add_definitions(-DVERBOSITY=${VERBOSITY})
add_definitions(-DUSE_MPI=${USE_MPI})

# Add support for MPI

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(imports/chimes_calculator/chimesFF/src)

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

############################################################
# Executables for chimes_lsq
############################################################

# Define an executable target

add_executable(chimes_lsq
    ${CMAKE_CURRENT_SOURCE_DIR}/src/chimes_lsq.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/functions.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Ewald.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ClassDefs.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Fcut.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Cluster.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Cheby.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io_styles.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/feexcept.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/A_Matrix.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/input.C
)

# Define properties for the executable target

target_include_directories(chimes_lsq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/)
target_compile_features   (chimes_lsq PRIVATE cxx_std_11)
target_link_libraries     (chimes_lsq ${MPI_LIBRARIES})
target_link_directories   (chimes_lsq PUBLIC imports/chimes_calculator/build)
target_link_libraries     (chimes_lsq libchimescalc.so)
set_target_properties(chimes_lsq PROPERTIES INSTALL_RPATH ${CMAKE_CURRENT_SOURCE_DIR}/imports/chimes_calculator/build/lib)

IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) 
	install(TARGETS chimes_lsq DESTINATION ${CMAKE_INSTALL_PREFIX}/src OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) 

############################################################
# Executables for chimes_md
############################################################

# Define an executable target

add_executable(chimes_md
    ${CMAKE_CURRENT_SOURCE_DIR}/src/chimes_md.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/functions.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Ewald.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ClassDefs.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Fcut.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Cluster.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/util.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Cheby.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io_styles.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/feexcept.C 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/input.C
)

# Define properties for the executable target

target_include_directories(chimes_md PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/)
target_compile_features   (chimes_md PRIVATE cxx_std_11)
target_link_libraries     (chimes_md ${MPI_LIBRARIES})
target_link_directories   (chimes_md PUBLIC imports/chimes_calculator/build)
target_link_libraries     (chimes_md libchimescalc.so)
set_target_properties(chimes_md PROPERTIES INSTALL_RPATH ${CMAKE_CURRENT_SOURCE_DIR}/imports/chimes_calculator/build/lib)

IF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) 
	install(TARGETS chimes_md DESTINATION ${CMAKE_INSTALL_PREFIX}/src OPTIONAL)
ENDIF (NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
