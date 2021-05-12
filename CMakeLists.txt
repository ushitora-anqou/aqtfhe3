cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(aqtfhe3 LANGUAGES CXX)

set(AQTFHE3_CXXFLAGS -Wall -Wextra -pedantic -Wno-sign-compare)
set(AQTFHE3_CXXFLAGS_DEBUG   ${IYOKAN_CXXFLAGS} -O0 -g3)
set(AQTFHE3_CXXFLAGS_RELEASE ${IYOKAN_CXXFLAGS} -Ofast -march=native -DNDEBUG -g3)
set(AQTFHE3_LIBS tbb)

set(CMAKE_C_FLAGS "")
set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_C_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELEASE "")

add_executable(main main.cpp test.cpp benchmark.cpp)
target_compile_options(main PUBLIC "$<$<CONFIG:DEBUG>:${AQTFHE3_CXXFLAGS_DEBUG}>")
target_compile_options(main PUBLIC "$<$<CONFIG:RELEASE>:${AQTFHE3_CXXFLAGS_RELEASE}>")
target_link_libraries(main ${AQTFHE3_LIBS})
