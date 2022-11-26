cmake_minimum_required(VERSION 3.24.2)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Add your project files
file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)

igl_include(imgui)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::imgui)

igl_include(png)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::png)
