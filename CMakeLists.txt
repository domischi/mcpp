PROJECT(mc++)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

enable_language(C CXX)
#-------------------------------------------------------------------
# Options
#-------------------------------------------------------------------
#if you don't want the full compiler output, set to "OFF"
option(CMAKE_VERBOSE_MAKEFILE "Verbose Makefile" OFF)
option(USE_MPI "Use distributed memory parallization" ON)
option(USE_OPENMP "Use shared memory parallelization" ON)
option(USE_FFTW "Use FFTW" ON)
option(USE_ARMADILLO "Use Armadillo" ON)
option(OPTIMIZE "Set optimization flags" ON)
option(DEBUG "Set debug flags" OFF)
option(DEBUG_VECTORIZATION "Set vectorization debug flags (attention, verbose)" OFF)
option(ALPS_HAS_VALARRAY_SET_ZERO "ALPS has valarray set_zero" ON)
#-------------------------------------------------------------------
# ALPS
#-------------------------------------------------------------------
# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

set(LIBS "${ALPS_LIBRARIES}")
if(ALPS_HAS_VALARRAY_SET_ZERO)
    add_definitions(-DALPS_HAS_VALARRAY_SET_ZERO)
endif(ALPS_HAS_VALARRAY_SET_ZERO)
#-------------------------------------------------------------------
# FFTW
#-------------------------------------------------------------------
if(USE_FFTW)
    find_path (FFTW_INCLUDES fftw3.h)
    find_library (FFTW_LIBRARIES NAMES fftw3)
    mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
    set(FFTW_LIBS "${FFTW_LIBRARIES}")
    if(FFTW_LIBRARIES)
        message(STATUS "Found FFTW")
    else(FFTW_LIBRARIES)
        message(STATUS "Did not find FFTW, build with trivial implementation")
        add_definitions(-DNFFTW) 
    endif(FFTW_LIBRARIES)
else(USE_FFTW)
    add_definitions(-DNFFTW) 
endif(USE_FFTW)

#-------------------------------------------------------------------
# Armadillo
#-------------------------------------------------------------------
if(USE_ARMADILLO)
    find_package(Armadillo)
    if(${ARMADILLO_FOUND})
        message(STATUS "Found Armadillo: v${ARMADILLO_VERSION_STRING}")
        include_directories(${ARMADILLO_INCLUDE_DIRS})
else(${ARMADILLO_FOUND})
    add_definitions(-DNARMADILLO)
endif(${ARMADILLO_FOUND})
else(USE_ARMADILLO)
    add_definitions(-DNARMADILLO)
endif(USE_ARMADILLO)
#-------------------------------------------------------------------
# Flags
#-------------------------------------------------------------------
SET(COMMON_FLAGS "-w")

IF(DEBUG)
    SET(COMMON_FLAGS "${COMMON_FLAGS} -pg -g -Wall")
    ELSE(DEBUG)
    SET(COMMON_FLAGS "${COMMON_FLAGS} -DNDEBUG")
ENDIF(DEBUG)
IF(DEBUG_VECTORIZATION)
    SET(COMMON_FLAGS "${COMMON_FLAGS} -fopt-info-vec-missed")
ENDIF(DEBUG_VECTORIZATION)

IF(OPTIMIZE)
    SET(COMMON_FLAGS "${COMMON_FLAGS} -O3 -ffast-math -ftree-vectorize -march=native")
ENDIF(OPTIMIZE)

IF(USE_MPI)
    SET(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK")
ENDIF(USE_MPI)

IF(USE_OPENMP)
    SET(COMMON_FLAGS "${COMMON_FLAGS} -fopenmp -DUSE_OPENMP")
ENDIF(USE_OPENMP)

SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS ,"cpp")

#-------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------
SET(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Single output directory for building all libraries")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Single output directory for building all executables")
SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}" CACHE PATH "Default install to home directory")
SET(PROJECT_CMAKE "${CMAKE_SOURCE_DIR}/CMake")
SET(BUILD_DIR "${CMAKE_SOURCE_DIR}/build")

#-------------------------------------------------------------------
# GIT VERSION (from pmirshad, github: https://github.com/pmirshad/cmake-with-git-metadata/)
#-------------------------------------------------------------------
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(GIT_COMMIT_HASH "")
endif(EXISTS "${CMAKE_SOURCE_DIR}/.git")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version/version.h.in
  ${CMAKE_SOURCE_DIR}/src/version/version.h
)
#-------------------------------------------------------------------
# mc++
#-------------------------------------------------------------------
SET(SRCS ${SRCS} src/main.cpp src/parapack_register.cpp)
ADD_EXECUTABLE(mc++ ${SRCS})
TARGET_LINK_LIBRARIES(mc++ ${LIBS})
TARGET_LINK_LIBRARIES(mc++ ${FFTW_LIBS})
TARGET_LINK_LIBRARIES(mc++ ${ARMADILLO_LIBRARIES})
INSTALL(TARGETS mc++ DESTINATION $ENV{HOME}/bin)
