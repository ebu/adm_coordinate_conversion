# Ensure we static link Catch to keep things simple
set(BUILD_SHARED_LIBS_ ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
add_subdirectory("submodules/Catch2")
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_})

add_executable(AdmCoordConvTests tests.cpp)
target_compile_features(AdmCoordConvTests PUBLIC cxx_std_17) 

target_link_libraries(AdmCoordConvTests
    PRIVATE
        Catch2WithMain
)

add_custom_command(TARGET AdmCoordConvTests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/test_data-pos_extent_cart2polar.txt" $<TARGET_FILE_DIR:AdmCoordConvTests>
)
if(WIN32 AND ${BUILD_SHARED_LIBS})
    add_custom_command(TARGET AdmCoordConvTests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:AdmCoordConv> $<TARGET_FILE_DIR:AdmCoordConvTests>
    )
endif()

set_target_properties(AdmCoordConvTests PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY $<TARGET_FILE_DIR:AdmCoordConvTests>
)

target_link_libraries(AdmCoordConvTests PRIVATE AdmCoordConv)

include(Catch)
catch_discover_tests(AdmCoordConvTests WORKING_DIRECTORY $<TARGET_FILE_DIR:AdmCoordConvTests>)

set_property(DIRECTORY PROPERTY USE_FOLDERS ON)
set_property(TARGET AdmCoordConvTests PROPERTY FOLDER "Test")
set_property(TARGET Catch2 PROPERTY FOLDER "Test")
set_property(TARGET Catch2WithMain PROPERTY FOLDER "Test")