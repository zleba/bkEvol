cmake_minimum_required (VERSION 2.6)
project (ChebMethod)

execute_process(COMMAND root-config --incdir OUTPUT_VARIABLE ROOTINC)
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOTLIB)
execute_process(COMMAND root-config --libs   OUTPUT_VARIABLE ROOTLIBS)

#Remove new line
string(REGEX REPLACE "\n$" "" ROOTINC "${ROOTINC}")
string(REGEX REPLACE "\n$" "" ROOTLIB "${ROOTLIB}")
string(REGEX REPLACE "\n$" "" ROOTLIBS "${ROOTLIBS}")

message(STATUS "inc:${ROOTINC}:")
message(STATUS "lib:${ROOTLIB}:")
message(STATUS "lis:${ROOTLIBS}:")

set(SRC src/cheb.cpp src/integration.cpp src/nodes.cpp src/Spline.cpp src/alphaSpline.cpp src/kernels.cpp)

if("$ENV{ARMADIR}" STREQUAL "")
    set(ARMAinc ../arma/install/include)
else()
    set(ARMAinc "$ENV{ARMADIR}/install/include")
endif()

add_definitions(-DARMA_DONT_USE_WRAPPER)

include_directories(inc  ${ROOTINC} ${ARMAinc})
#link_directories(${ROOTLIB})

add_library(base STATIC ${SRC})
target_compile_features(base PRIVATE cxx_range_for cxx_final cxx_attribute_deprecated)
add_executable(cheb src/main.cpp)
add_executable(one src/one.cpp)
add_executable(solveOne src/solveOne.cpp)
target_link_libraries(cheb base  -lopenblas ${ROOTLIBS})
target_link_libraries(one base  -lopenblas ${ROOTLIBS})
target_link_libraries(solveOne base  -lopenblas ${ROOTLIBS})
