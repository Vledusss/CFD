cmake_minimum_required(VERSION 3.28)
project(CFD LANGUAGES CXX C)
message("Project name: " ${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(third_party)
add_subdirectory(src)

add_executable(main main.cpp)

target_link_libraries(main ${PROJECT_NAME})
