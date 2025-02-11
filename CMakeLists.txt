cmake_minimum_required(VERSION 3.15)
project(AdmCoordConv VERSION 1.0.0 LANGUAGES CXX)

# Option to choose static or shared library
option(BUILD_SHARED_LIBS "Build shared (dynamic) libraries instead of static libraries" OFF)

# Library target
add_library(AdmCoordConv)
add_library(ADM::AdmCoordConv ALIAS AdmCoordConv)
if(BUILD_SHARED_LIBS)
    target_compile_definitions(AdmCoordConv PRIVATE ADM_COORD_CONV_EXPORT BUILD_SHARED_LIBS)
endif()

include(GenerateExportHeader)
generate_export_header(AdmCoordConv EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/adm_coord_conv_export.hpp)

# Source files
target_sources(AdmCoordConv
    PRIVATE
        src/adm_coord_conv.cpp
        include/adm_coord_conv/adm_coord_conv.hpp
)

# Include directories
target_include_directories(AdmCoordConv
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include> # Export header location
        $<INSTALL_INTERFACE:include>
)

set_target_properties(AdmCoordConv PROPERTIES EXPORT_NAME AdmCoordConv)

# C++ Standard
target_compile_features(AdmCoordConv PUBLIC cxx_std_17) 

# Install and export config
include(GNUInstallDirs)
install(TARGETS AdmCoordConv
    EXPORT AdmCoordConvTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT AdmCoordConvTargets
    FILE AdmCoordConvTargets.cmake
    NAMESPACE ADM::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AdmCoordConv
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/AdmCoordConvConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/AdmCoordConvConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/AdmCoordConvConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AdmCoordConv
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/AdmCoordConvConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/AdmCoordConvConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AdmCoordConv
)

# Option to build tests
option(BUILD_TESTS "Build the unit tests" OFF)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()