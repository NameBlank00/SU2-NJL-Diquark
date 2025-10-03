cmake_minimum_required(VERSION 3.1)

# set(CMAKE_C_COMPILER "gcc")
# set(CMAKE_CXX_COMPILER "g++")

project(TMQGP)

include(GenerateExportHeader)
# find_package(libInterpolate REQUIRED)
# find_package(libGSL REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_STANDARD 14)

if (UNIX)
set(SWIG_EXEC swig)
endif(UNIX)

find_package(Python3 COMPONENTS Interpreter Development NumPy)


if(UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g -Ofast -std=c++2a ")
endif(UNIX)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi /MP3")
endif()

set(VERBOSE TRUE)

file (GLOB_RECURSE libTMQGP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP/*.h
)

file (GLOB_RECURSE libTMQGP_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP/*.h
)

if (UNIX)
# include_directories("/usr/include/levmar")
include_directories(${Python3_INCLUDE_DIRS})
include_directories("/opt/anaconda3/include/")
include_directories(${Python3_NumPy_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/TMQGP)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)

# include_directories("/opt/homebrew/Cellar/gsl/2.8/include")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

include_directories("/usr/local/include")
# include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/")

# include_directories("/opt/homebrew/Cellar/boost/1.86.0/include/")
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eosstudio)

link_directories("/usr/lib/levmar")
link_directories("/usr/lib/x86_64-linux-gnu/")
# link_directories("/opt/homebrew/Cellar/gsl/2.8/lib")

# link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eosstudio/lib)

endif (UNIX)

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_library(libTMQGP SHARED ${libTMQGP_SRC})

set_target_properties(libTMQGP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# set_property(TARGET EosStudio
#    PROPERTY VS_PACKAGE_REFERENCES "boost_1.7.0"
# )

add_library(TMQGP SHARED ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i)

# add_executable(DirectUse ${CMAKE_CURRENT_SOURCE_DIR}/DirectUse.cpp)

if (APPLE)
set_property(TARGET libTMQGP APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
set_property(TARGET DirectUse APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")

set_property(TARGET TMQGP APPEND PROPERTY LINK_FLAGS "-flat_namespace -undefined suppress")
endif(APPLE)
# add_executable(TMQGP_mpi ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP_mpi.cpp)


set_target_properties(TMQGP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

set_target_properties(TMQGP PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set_target_properties(TMQGP PROPERTIES PREFIX _)
set_target_properties(TMQGP PROPERTIES IMPORT_PREFIX _)
#set_target_properties(TMQGP PROPERTIES RELW "_rd")
if (WIN32)
set_target_properties(TMQGP PROPERTIES SUFFIX .pyd)
endif (WIN32)
set_target_properties(TMQGP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
add_dependencies(TMQGP libTMQGP)
# add_dependencies(DirectUse libTMQGP)





#set_property(TARGET EosStudio
#    PROPERTY VS_PACKAGE_REFERENCES "stacktrace_0.1"
#)

#set_property(TARGET EosStudio
#    PROPERTY VS_PACKAGE_REFERENCES "stacktrace_0.1"
#)

if (UNIX)
target_link_libraries(TMQGP ${GSL_LIBRARIES} lapack ${Python3_LIBRARIES} libTMQGP ${MPI_CXX_LIBRARIES} Eigen3::Eigen)
target_link_libraries(libTMQGP ${GSL_LIBRARIES} lapack ${Python3_LIBRARIES} ${MPI_CXX_LIBRARIES} Eigen3::Eigen)
# target_link_libraries(DirectUse ${GSL_LIBRARIES} lapack ${Python3_LIBRARIES} libTMQGP ${MPI_CXX_LIBRARIES} Eigen3::Eigen)
# target_link_libraries(TMQGP_mpi gsl gslcblas lapack ${Python3_LIBRARIES} libTMQGP ${MPI_CXX_LIBRARIES})
endif(UNIX)


if (WIN32)
target_link_libraries(TMQGP gsl cblas blas blas_nowrap levmar lapack python38 vcf2c libTMQGP gslbase Tdbase StackTrace)
target_link_libraries(libTMQGP gsl cblas blas blas_nowrap levmar lapack python38 vcf2c gslbase Tdbase StackTrace)
endif(WIN32)

# WES stands for wrapEosStudio
get_target_property(WES_RUNTIME_OUTPUT_DIRECTORY TMQGP RUNTIME_OUTPUT_DIRECTORY)

if(WIN32)
link_directories(${LIBRARY_OUTPUT_PATH})
link_directories(${libTMQGP_RUNTIME_OUTPUT_DIRECTORY})
endif(WIN32)

#set (SWIGCOMMAND "${SWIG_EXEC} -DSWIG_PYTHON_LEGACY_BOOL -I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio -c++ -python -o ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp 
#-outdir ${EosStudio_RUNTIME_OUTPUT_DIRECTORY}/ ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i" )
set(CMAKE_VERBOSE_MAKEFILE on)
#-I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio 
# -IC:/Users/const.maslov/Numerics/EosStudio/packages/stacktrace.0.1/build/native/include/


find_package(Doxygen)

# set(SWIG_DOXYGEN "")

if (DOXYGEN_FOUND)
    set(DOXYGEN_GENERATE_HTML YES)


    doxygen_add_docs(
        doxygen
        ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP
        ${CMAKE_CURRENT_SOURCE_DIR}/EosStudio
        ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i
        COMMENT "Generating doxygen documentation"
    )

    # set(SWIG_DEGN)
endif(DOXYGEN_FOUND)



add_custom_command(
OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp
	#message(WARNING "Swig_exec=${SWIG_EXEC}")
COMMAND ${SWIG_EXEC} -DSWIG_PYTHON_LEGACY_BOOL -I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio -c++ -python -keyword -o ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp -outdir ${CMAKE_CURRENT_SOURCE_DIR}/lib/ ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i
DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i
DEPENDS ${libTMQGP_HEADERS}
)



#COMMAND ${SWIG_EXEC} -DSWIG_PYTHON_LEGACY_BOOL -I${GOVNODIR} -I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio -c++ -python -o ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp 
 #       -outdir ${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_CFG_INTDIR} ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i

# add_executable(DirectUse DirectUse.cpp)

# if (UNIX)
# target_link_libraries(DirectUse gsl gslcblas levmar lapack EosStudio python3.10)
# endif(UNIX)

# if(WIN32)
# target_link_libraries(DirectUse gsl cblas blas blas_nowrap levmar lapack libTMQGP python38 vcf2c gslbase Tdbase)
# endif(WIN32)
