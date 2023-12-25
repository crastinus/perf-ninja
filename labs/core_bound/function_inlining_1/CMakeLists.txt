cmake_minimum_required(VERSION 3.3)

project(lab)


add_compile_options(-fsanitize=address)
add_link_options(-fsanitize=address)



string(REGEX MATCH "^(.*)[\\/]labs[\\/].*$" repo "${CMAKE_CURRENT_SOURCE_DIR}")
include(${CMAKE_MATCH_1}/tools/labs.cmake)
