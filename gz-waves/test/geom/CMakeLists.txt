# Golden-value tool for tests/geom/test_geom_baseline.py (CGAL removal).
# Built with the unit tests; it links only the core library.
add_executable(geom_baseline geom_baseline.cc)
target_link_libraries(geom_baseline
  PRIVATE
    ${PROJECT_LIBRARY_TARGET_NAME}
    gz-common${GZ_COMMON_VER}::graphics
)
set_target_properties(geom_baseline PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
