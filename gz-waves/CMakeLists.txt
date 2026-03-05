cmake_minimum_required(VERSION 3.10.2)

#============================================================================
# Initialize the project
#============================================================================
project(gz-waves1 VERSION 1.0.0)

#============================================================================
# Find gz-cmake
#============================================================================
if("$ENV{GZ_VERSION}" STREQUAL "garden"
    OR "$ENV{GZ_VERSION}" STREQUAL "harmonic"
    OR NOT DEFINED "ENV{GZ_VERSION}")
  find_package(gz-cmake3 REQUIRED)
  set(GZ_CMAKE_VER ${gz-cmake3_VERSION_MAJOR})
elseif("$ENV{GZ_VERSION}" STREQUAL "ionic")
  find_package(gz-cmake4 REQUIRED)
  set(GZ_CMAKE_VER ${gz-cmake4_VERSION_MAJOR})
endif()

#============================================================================
# Configure the project
#============================================================================set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()
gz_configure_project()

#============================================================================
# Include cmake
#============================================================================

include(${PROJECT_SOURCE_DIR}/cmake/URL.conf.cmake)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#============================================================================
# Set project-specific options
#============================================================================

#============================================================================
# Search for project-specific dependencies
#============================================================================

#--------------------------------------
# Find gnuplot-iostream - the header is used by the examples in test/plots
include(Add_gnuplot-iostream)

# Garden (default)
if("$ENV{GZ_VERSION}" STREQUAL "garden" OR NOT DEFINED "ENV{GZ_VERSION}")

  #--------------------------------------
  # Find gz-math
  gz_find_package(gz-math7 REQUIRED COMPONENTS eigen3)
  set(GZ_MATH_VER ${gz-math7_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-plugin
  gz_find_package(gz-plugin2 REQUIRED COMPONENTS loader register)
  set(GZ_PLUGIN_VER ${gz-plugin2_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-common
  gz_find_package(gz-common5
    REQUIRED COMPONENTS graphics events
  )
  set(GZ_COMMON_VER ${gz-common5_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-msgs
  gz_find_package(gz-msgs9 REQUIRED)
  set(GZ_MSGS_VER ${gz-msgs9_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-transport
  gz_find_package(gz-transport12 REQUIRED)
  set(GZ_TRANSPORT_VER ${gz-transport12_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-sensors
  gz_find_package(gz-sensors7 REQUIRED)
  set(GZ_SENSORS_VER ${gz-sensors7_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-rendering
  gz_find_package(gz-rendering7 REQUIRED)
  set(GZ_RENDERING_VER ${gz-rendering7_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-sim
  gz_find_package(gz-sim7 REQUIRED)
  set(GZ_SIM_VER ${gz-sim7_VERSION_MAJOR})

  #--------------------------------------
  # Find SDFormat
  gz_find_package(sdformat13 REQUIRED)
  set(SDF_VER ${sdformat13_VERSION_MAJOR})

  message(STATUS "Compiling against Gazebo Garden")

# Harmonic
elseif("$ENV{GZ_VERSION}" STREQUAL "harmonic")

  #--------------------------------------
  # Find gz-math
  gz_find_package(gz-math7 REQUIRED COMPONENTS eigen3)
  set(GZ_MATH_VER ${gz-math7_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-plugin
  gz_find_package(gz-plugin2 REQUIRED COMPONENTS loader register)
  set(GZ_PLUGIN_VER ${gz-plugin2_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-common
  gz_find_package(gz-common5
    REQUIRED COMPONENTS graphics events
  )
  set(GZ_COMMON_VER ${gz-common5_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-msgs
  gz_find_package(gz-msgs10 REQUIRED)
  set(GZ_MSGS_VER ${gz-msgs10_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-transport
  gz_find_package(gz-transport13 REQUIRED)
  set(GZ_TRANSPORT_VER ${gz-transport13_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-sensors
  gz_find_package(gz-sensors8 REQUIRED)
  set(GZ_SENSORS_VER ${gz-sensors8_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-rendering
  gz_find_package(gz-rendering8 REQUIRED)
  set(GZ_RENDERING_VER ${gz-rendering8_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-sim
  gz_find_package(gz-sim8 REQUIRED)
  set(GZ_SIM_VER ${gz-sim8_VERSION_MAJOR})

  #--------------------------------------
  # Find SDFormat
  gz_find_package(sdformat14 REQUIRED)
  set(SDF_VER ${sdformat14_VERSION_MAJOR})

  message(STATUS "Compiling against Gazebo Harmonic")

# Ionic
elseif("$ENV{GZ_VERSION}" STREQUAL "ionic")

  #--------------------------------------
  # Find gz-math
  gz_find_package(gz-math8 REQUIRED COMPONENTS eigen3)
  set(GZ_MATH_VER ${gz-math8_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-plugin
  gz_find_package(gz-plugin3 REQUIRED COMPONENTS loader register)
  set(GZ_PLUGIN_VER ${gz-plugin3_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-common
  gz_find_package(gz-common6
    REQUIRED COMPONENTS graphics events
  )
  set(GZ_COMMON_VER ${gz-common6_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-msgs
  gz_find_package(gz-msgs11 REQUIRED)
  set(GZ_MSGS_VER ${gz-msgs11_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-transport
  gz_find_package(gz-transport14 REQUIRED)
  set(GZ_TRANSPORT_VER ${gz-transport14_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-sensors
  gz_find_package(gz-sensors9 REQUIRED)
  set(GZ_SENSORS_VER ${gz-sensors9_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-rendering
  gz_find_package(gz-rendering9 REQUIRED)
  set(GZ_RENDERING_VER ${gz-rendering9_VERSION_MAJOR})

  #--------------------------------------
  # Find gz-sim
  gz_find_package(gz-sim9 REQUIRED)
  set(GZ_SIM_VER ${gz-sim9_VERSION_MAJOR})

  #--------------------------------------
  # Find SDFormat
  gz_find_package(sdformat15 REQUIRED)
  set(SDF_VER ${sdformat15_VERSION_MAJOR})

  message(STATUS "Compiling against Gazebo Ionic")

else()
  message(FATAL_ERROR "Unsupported GZ_VERSION: $ENV{GZ_VERSION}")
endif()

#--------------------------------------
if("$ENV{ROS_DISTRO}" STREQUAL "jazzy")
  # Resolve https://github.com/srmainwaring/asv_wave_sim/issues/168 by applying the CMakeLists that provided by SeapathMasteron on Jul 4, 2024
  # Manually set OGRE2 paths
  set(OGRE2_INCLUDE_DIRS "/opt/ros/jazzy/opt/gz_ogre_next_vendor/include/OGRE-Next")
  set(OGRE2_LIBRARY_DIRS "/opt/ros/jazzy/opt/gz_ogre_next_vendor/lib")
  include_directories(${OGRE2_INCLUDE_DIRS})
  include_directories("${OGRE2_INCLUDE_DIRS}/Hlms/Common")
  link_directories(${OGRE2_LIBRARY_DIRS})

  # Manually create GzOGRE2 target
  add_library(GzOGRE2::GzOGRE2 INTERFACE IMPORTED)
  set_target_properties(GzOGRE2::GzOGRE2 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OGRE2_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${OGRE2_LIBRARY_DIRS}/libOgreNextMain.so"
  )

  # Manually create other OGRE2 components targets
  add_library(GzOGRE2::HlmsPbs INTERFACE IMPORTED)
  set_target_properties(GzOGRE2::HlmsPbs PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OGRE2_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${OGRE2_LIBRARY_DIRS}/libOgreNextHlmsPbs.so"
  )

  add_library(GzOGRE2::HlmsUnlit INTERFACE IMPORTED)
  set_target_properties(GzOGRE2::HlmsUnlit PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OGRE2_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${OGRE2_LIBRARY_DIRS}/libOgreNextHlmsUnlit.so"
  )

  add_library(GzOGRE2::Overlay INTERFACE IMPORTED)
  set_target_properties(GzOGRE2::Overlay PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OGRE2_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${OGRE2_LIBRARY_DIRS}/libOgreNextOverlay.so"
  )

  add_library(GzOGRE2::PlanarReflections INTERFACE IMPORTED)
  set_target_properties(GzOGRE2::PlanarReflections PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${OGRE2_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${OGRE2_LIBRARY_DIRS}/libOgreNextPlanarReflections.so"
  )
else()
  # Find OGRE
  list(APPEND gz_ogre_components "RTShaderSystem" "Terrain" "Overlay" "Paging")

  gz_find_package(GzOGRE VERSION 1.9.0
    COMPONENTS ${gz_ogre_components}
    REQUIRED_BY ogre
    PRIVATE_FOR ogre)
  #--------------------------------------
  # Find OGRE2: first try to find OGRE2 built with PlanarReflections support and
  # fallback to look for OGRE2 without it. Both seems to works for gz-rendering.
  # See https://github.com/gazebosim/gz-rendering/issues/597
  gz_find_package(GzOGRE2 VERSION 2.3
    COMPONENTS HlmsPbs HlmsUnlit Overlay PlanarReflections
    PRIVATE_FOR ogre2
    QUIET)

  if ("${OGRE2-PlanarReflections}" STREQUAL "OGRE2-PlanarReflections-NOTFOUND")
    message(STATUS "PlanarReflections component was not found. Try looking without it:")
    gz_find_package(GzOGRE2 VERSION 2.3
      COMPONENTS HlmsPbs HlmsUnlit Overlay
      REQUIRED_BY ogre2
      PRIVATE_FOR ogre2)
  endif()

  if (OGRE2_FOUND)
    set(HAVE_OGRE2 TRUE)
  endif()
endif()
#####################################

# Define compile-time default variables
if(MSVC)
  set(GZ_RENDERING_PLUGIN_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
else()
  set(GZ_RENDERING_PLUGIN_PATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()

#--------------------------------------
# Find eigen3

find_package(Eigen3 REQUIRED)

#--------------------------------------
# CGAL and its components...

find_package(CGAL REQUIRED COMPONENTS Core)

if(NOT CGAL_FOUND)
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

message("CGAL_VERSION_MAJOR: ${CGAL_VERSION_MAJOR}")

if (${CGAL_VERSION_MAJOR} LESS 6)
  # include helper file
  include(${CGAL_USE_FILE})
endif()

#--------------------------------------
# Find FFTW3 (double-precision) (GPL) as FFT library.

find_path(FFTW3-3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS $ENV{FFTW3_DIR}/include
  PATHS /usr/local/include
        /usr/include
)

find_library(FFTW3-3_LIBRARY
  NAMES fftw3 libfftw3
  HINTS $ENV{FFTW3_DIR}/lib
  PATHS /usr/local/lib
        /usr/lib
)

add_definitions(-DUSE_FFTW3)
set(FFT_INCLUDE_DIRS ${FFTW3-3_INCLUDE_DIR})
set(FFT_LIBRARIES ${FFTW3-3_LIBRARY})

#--------------------------------------
# Find OpenCL

# find_package(OpenCL REQUIRED)

#--------------------------------------
# Find clFFT

# find_package(clFFT REQUIRED)

#============================================================================

# Location of "fake install folder" used in tests
# Defined here at root scope so it is available for tests in src and test folders
set(FAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/fake/install")

#============================================================================
# Configure the build
#============================================================================

gz_configure_build(QUIT_IF_BUILD_ERRORS)

# add_subdirectory(examples)

#============================================================================
# Create package information
#============================================================================
gz_create_packages()

#============================================================================
# Configure documentation
#============================================================================

# TODO: add documentation api.md.in files

find_package(ament_cmake REQUIRED)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

# Following hooks are used to ensure that the correct environment variables
# will be set by exectuting 'source install/setup.bash' after compilation.
# When using this template for your project, change the filenames of the
# files in the 'hooks' folder, to correspond to your project name.
ament_environment_hooks("${CMAKE_CURRENT_SOURCE_DIR}/hooks/${PROJECT_NAME}.dsv.in")
ament_environment_hooks("${CMAKE_CURRENT_SOURCE_DIR}/hooks/${PROJECT_NAME}.sh.in")

ament_package()