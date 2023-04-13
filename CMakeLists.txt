cmake_minimum_required(VERSION 3.16)
project(cg-shop-2023)

set(CMAKE_CXX_STANDARD 20)

include("${PROJECT_SOURCE_DIR}/conanbuildinfo.cmake")
conan_basic_setup()

file(GLOB files "${PROJECT_SOURCE_DIR}/src/*.cpp")
foreach (file ${files})
    get_filename_component(version "${file}" NAME_WE)
    add_executable("${version}" "src/${version}.cpp")
    target_link_libraries("${version}" ${CONAN_LIBS})
    target_include_directories("${version}" PRIVATE include)
endforeach ()
