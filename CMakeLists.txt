cmake_minimum_required(VERSION 3.5)

project(open-pic LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


file(GLOB TARGET_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME} ${TARGET_SRC})

target_compile_options(${PROJECT_NAME} PRIVATE -O3 -w -ffast-math)

target_include_directories(${PROJECT_NAME} PRIVATE
                           ${CMAKE_SOURCE_DIR}/3rdparty/boost_1_52_0 
                           ${CMAKE_SOURCE_DIR}/3rdparty/lua-5.1/src
                           ${CMAKE_SOURCE_DIR}/3rdparty/lua-5.1/etc
                           ${CMAKE_SOURCE_DIR}/3rdparty/luabind-0.8.1/lib/include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_directories(${PROJECT_NAME} PRIVATE
                        ${CMAKE_SOURCE_DIR}/3rdparty/luabind-0.8.1/stage
                        /usr/lib)

target_link_libraries(${PROJECT_NAME} PRIVATE
                      lua5.1 luabind)
