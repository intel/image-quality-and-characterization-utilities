cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Teisko)

if (WIN32)
    set(ENABLE_WARNINGS_AS_ERRORS /W4 /WX)
    if (CMAKE_SIZEOF_VOID_P EQUAL 4)
        add_definitions(-DNO_64_BIT_REGS)
    endif()
elseif (UNIX)
    set(ENABLE_WARNINGS_AS_ERRORS -Wall -Wextra -Werror)

    if (CMAKE_SIZEOF_VOID_P EQUAL 4)
        set(GCC_ARCH_OPTION "-m32 -DNO_64_BIT_REGS")
    else()
        set(GCC_ARCH_OPTION "-m64")
    endif()

    set(CMAKE_CXX_FLAGS "-std=c++11 ${GCC_ARCH_OPTION} -msse4")
else()
    message(FATAL_ERROR "Unsupported system")
endif()

# TODO: Add inside macro when Eigen is taken in use
add_definitions(-DEIGEN_MPL2_ONLY)

set(HEADER_FILES
    # Algorithms
    include/Teisko/Algorithm/Bit.hpp
    include/Teisko/Algorithm/ConvexHull.hpp
    include/Teisko/Algorithm/DelaunayTriangulation.hpp
    include/Teisko/Algorithm/Functors.hpp
    include/Teisko/Algorithm/Histogram.hpp
    include/Teisko/Algorithm/Interpolate.hpp
    include/Teisko/Algorithm/Iterators.hpp
    include/Teisko/Algorithm/LinearSpace.hpp
    include/Teisko/Algorithm/NelderMead.hpp
    include/Teisko/Algorithm/PointXY.hpp
    include/Teisko/Algorithm/Pow2.hpp
    include/Teisko/Algorithm/ReduceTo.hpp
    include/Teisko/Algorithm/StandardDeviation.hpp
    include/Teisko/Algorithm/TrimmedMean.hpp
    include/Teisko/Algorithm/VectorMedian.hpp
    include/Teisko/Algorithm/VectorOperations.hpp
    # Data
    include/Teisko/Data/MacbethReflectance.hpp
    include/Teisko/Data/MunsellReflectance.hpp
    # Image related
    include/Teisko/Image/Algorithms.hpp
    include/Teisko/Image/API.hpp
    include/Teisko/Image/Conversion.hpp
    include/Teisko/Image/Point.hpp
    include/Teisko/Image/Polyscale.hpp
    include/Teisko/Image/RGB.hpp
    include/Teisko/Image/Support.hpp
    include/Teisko/Image/TIFF.hpp
    # Models and others
	include/Teisko/CalculateLux.hpp
    include/Teisko/BayerImage.hpp
    include/Teisko/BayerInfo.hpp
    include/Teisko/Chromaticity.hpp
    include/Teisko/Color.hpp
    include/Teisko/ColorCorrection.hpp
    include/Teisko/LateralChromaticAberration.hpp
    include/Teisko/LensShading.hpp
    include/Teisko/MacbethDetector.hpp
    include/Teisko/NoiseModel.hpp
    include/Teisko/Preprocessing.hpp
    include/Teisko/SpectralResponse.hpp
    include/Teisko/ValidImageArea.hpp
    )

set(TEST_FILES
    # Test runner and test data
    tests/catch.cpp
    tests/_data/icc_profiles.cpp
    tests/_data/icc_profiles.hpp
    # Algorithms
    tests/Algorithms/specs_bit.cpp
    tests/Algorithms/specs_convex_hull.cpp
    tests/Algorithms/specs_delaunay_triangulation.cpp
    tests/Algorithms/specs_interpolate.cpp
    tests/Algorithms/specs_linear_space.cpp
    tests/Algorithms/specs_nelder_mead.cpp
    tests/Algorithms/specs_trimmed_mean.cpp
    tests/Algorithms/specs_vector_median.cpp
    # Image related
    tests/Image/specs_algorithms.cpp
    tests/Image/specs_api.cpp
    tests/Image/specs_conversion.cpp
    tests/Image/specs_point.cpp
    tests/Image/specs_polyscale.cpp
    tests/Image/specs_recipes.cpp
    tests/Image/specs_rgb.cpp
    tests/Image/specs_support.cpp
    tests/Image/specs_tiff.cpp
    # Models and others
    tests/specs_calculate_lux.cpp
    tests/specs_bayer_info.cpp
    tests/specs_color.cpp
    tests/specs_color_correction.cpp
    tests/specs_lateral_chromatic_aberration.cpp
    tests/specs_lens_shading.cpp
    tests/specs_macbeth_detector.cpp
    tests/specs_noise_model.cpp
    tests/specs_preprocessing.cpp
    tests/specs_spectral_response.cpp
    tests/specs_valid_image_area.cpp
    tests/test_linker.cpp
    )

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/include/Teisko"
    PREFIX "Teisko"
    FILES  ${HEADER_FILES}
    )
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    PREFIX "Tests"
    FILES  ${TEST_FILES}
    )
add_executable(teisko_tester
    ${HEADER_FILES}
    ${TEST_FILES}
    )
set_target_properties(teisko_tester
    PROPERTIES DEBUG_POSTFIX d
    )
target_compile_options(teisko_tester
    PRIVATE ${ENABLE_WARNINGS_AS_ERRORS}
    )
target_include_directories(teisko_tester
    PRIVATE include
            tests
    )

# This property makes the test executable project the startup project in Visual
# Studio solutions. Requires CMake 3.6.3
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    PROPERTY VS_STARTUP_PROJECT teisko_tester
    )

# Code coverage is measured only when
# (a) opencppcoverage.exe is available, and
# (b) debug information (.pdb files) are to be generated.
#
# If opencppcoverage.exe could not be found or no .pdb files are about to be
# generated then we just run the tests as is.

enable_testing()
find_program(OPENCPPCOVERAGE_EXE opencppcoverage.exe)

if (OPENCPPCOVERAGE_EXE)
    message("-- Using code coverage tool: ${OPENCPPCOVERAGE_EXE}")
    string(REPLACE "/" "\\"
        PROJECT_SOURCES_WINDOWS_PATH
        "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    add_test(NAME ${CMAKE_PROJECT_NAME}-coverage-tester
        CONFIGURATIONS Debug | RelWithDebInfo
        COMMAND ${OPENCPPCOVERAGE_EXE}
            --export_type cobertura
            --sources ${PROJECT_SOURCES_WINDOWS_PATH}
            --excluded_sources ${PROJECT_SOURCES_WINDOWS_PATH}\\include\\Eigen
            --excluded_sources ${PROJECT_SOURCES_WINDOWS_PATH}\\tests
            --excluded_modules $ENV{SystemRoot}
            $<TARGET_FILE:teisko_tester>
        )
    add_test(NAME ${CMAKE_PROJECT_NAME}-tester
        CONFIGURATIONS Release | MinSizeRel
        COMMAND $<TARGET_FILE:teisko_tester>
        )
else()
    message("-- No code coverage tool found")
    add_test(NAME ${CMAKE_PROJECT_NAME}-tester
        COMMAND $<TARGET_FILE:teisko_tester>
        )
endif()

