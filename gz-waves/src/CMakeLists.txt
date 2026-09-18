add_subdirectory(gui)
add_subdirectory(systems)

# Collect source files into the "sources" variable and unit test files into the
# "gtest_sources" variable.
# gz_get_libsources_and_unittests(sources gtest_sources)

# Collect source and test files manually

set(sources
  Algorithm.cc
  geom/RayMeshQuery.cc
  Convert.cc
  Geometry.cc
  Grid.cc
  MeshTools.cc
  OceanTile.cc
  PhysicalConstants.cc
  Physics.cc
  SubMeshWithTangents.cc
  TriangulatedGrid.cc
  Utilities.cc
  WaterCurrentGrid.cc
  Wavefield.cc
  WavefieldSampler.cc
  WaveParameters.cc
  WaveSimulation.cc
  LinearRandomFFTWaveSimulation.cc
  LinearRandomFFTWaveSimulationRef.cc
  fft/Fft.cc
  LinearRandomWaveSimulation.cc
  LinearRegularWaveSimulation.cc
  TrochoidIrregularWaveSimulation.cc
  WaveSpectrum.cc
  WaveSpreadingFunction.cc
)

set(gtest_sources
  ${gtest_sources}
  Algorithm_TEST.cc
  Fft_TEST.cc
  Geometry_TEST.cc
  Grid_TEST.cc
  LinearRegularWaveSimulation_TEST.cc
  LinearRandomFFTWaveSimulation_TEST.cc
  MeshTools_TEST.cc
  Physics_TEST.cc
  TriangulatedGrid_TEST.cc
  Wavefield_TEST.cc
  WaveSpectrum_TEST.cc
  WaveSpreadingFunction_TEST.cc
)

# Create the library target
gz_create_core_library(SOURCES ${sources} CXX_STANDARD ${CMAKE_CXX_STANDARD})

target_link_libraries(${PROJECT_LIBRARY_TARGET_NAME}
  PUBLIC
    gz-math${GZ_MATH_VER}
    gz-msgs${GZ_MSGS_VER}::gz-msgs${GZ_MSGS_VER}
    gz-common${GZ_COMMON_VER}::gz-common${GZ_COMMON_VER}
    gz-common${GZ_COMMON_VER}::graphics
    gz-common${GZ_COMMON_VER}::profiler
    gz-transport${GZ_TRANSPORT_VER}::gz-transport${GZ_TRANSPORT_VER}
    sdformat${SDF_VER}::sdformat${SDF_VER}
  PRIVATE
    OpenMP::OpenMP_CXX
)
if(GZ_WAVES_WITH_EMBREE)
  target_link_libraries(${PROJECT_LIBRARY_TARGET_NAME} PRIVATE embree)
  target_compile_definitions(${PROJECT_LIBRARY_TARGET_NAME}
    PRIVATE GZ_WAVES_HAVE_EMBREE=1)
endif()
if (UNIX AND NOT APPLE)
  target_link_libraries(${PROJECT_LIBRARY_TARGET_NAME}
    PRIVATE stdc++fs)
endif()

target_include_directories(${PROJECT_LIBRARY_TARGET_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  ${EIGEN3_INCLUDE_DIR}
  PRIVATE
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/thirdparty/pocketfft
)

include_directories(${PROJECT_SOURCE_DIR}/test)

# Build the unit tests
gz_build_tests(TYPE UNIT
  SOURCES
    ${gtest_sources}
  LIB_DEPS
    ${PROJECT_LIBRARY_TARGET_NAME}
)
