cmake_minimum_required(VERSION 3.12)
project(SynDiv_c CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加所有的库
set(LIBRARIES strip_split_join save get_time multiinter coor no_syn cal window)

foreach(library IN LISTS LIBRARIES)
    add_library(${library} SHARED src/${library}.cpp)
    target_compile_features(${library} PRIVATE cxx_std_17)
    target_compile_options(${library} PRIVATE -Wall -Wextra -O3)
endforeach()

# 添加可执行文件
add_executable(SynDiv_c SynDiv_c.cpp)

# 添加头文件搜索路径
target_include_directories(SynDiv_c PRIVATE include)

# 查找threads库和zlib库
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# 设置公共链接库列表
set(COMMON_LIBRARIES 
    strip_split_join
    get_time
    save
    Threads::Threads
    ZLIB::ZLIB
    stdc++fs
)

# 添加链接库
target_link_libraries(SynDiv_c PRIVATE 
    multiinter
    coor
    no_syn
    cal
    window
    ${COMMON_LIBRARIES}	
)

# 添加编译选项
target_compile_options(SynDiv_c PRIVATE -O3)