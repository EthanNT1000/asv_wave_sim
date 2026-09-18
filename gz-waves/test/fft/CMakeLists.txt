# Golden-value tool for tests/fft/test_fft_baseline.py (FFTW -> PocketFFT).
# Built with the unit tests; it links only the core library.
add_executable(fft_baseline fft_baseline.cc)
target_link_libraries(fft_baseline
  PRIVATE
    ${PROJECT_LIBRARY_TARGET_NAME}
)
target_include_directories(fft_baseline
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src
)
set_target_properties(fft_baseline PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
