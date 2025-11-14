file(GLOB_RECURSE SOURCE *.cpp *.hpp)

add_library(${PROJECT_NAME} STATIC ${SOURCE})

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(${PROJECT_NAME} PUBLIC eigen)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})