cmake_minimum_required(VERSION 3.12)

project(SynDiv_c CXX)

# 设置 C++ 标准版本 (这里使用 C++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 添加源文件
set(SRC_FILES_SynDiv_c
    src/get_time.cpp
    src/save.cpp
    src/strip_split_join.cpp
    src/multiinter.cpp
    src/no_syn.cpp
    src/coor.cpp
    src/cal.cpp
    src/window.cpp
)

# 添加可执行文件
add_executable(SynDiv_c SynDiv_c.cpp)

# 添加源文件
target_sources(SynDiv_c PRIVATE ${SRC_FILES_SynDiv_c})

# 为目标添加编译选项
target_compile_options(SynDiv_c PRIVATE -O3)

# 为目标添加编译选项（包括 -O3、额外的警告和错误和CPU优化）
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # target_compile_options(SynDiv_c PRIVATE -O3 -Wall -Wextra -Werror -pedantic-errors -march=native)
    # target_compile_options(SynDiv_c PRIVATE -O3 -Wall -Wextra -march=native)
    target_compile_options(SynDiv_c PRIVATE -O3 -Wall -Wextra)
endif()

# 添加 include 目录
target_include_directories(SynDiv_c PRIVATE include)

# 添加链接库
find_package(Threads REQUIRED)
target_link_libraries(SynDiv_c Threads::Threads)

find_package(ZLIB REQUIRED)
target_link_libraries(SynDiv_c ZLIB::ZLIB)

target_link_libraries(SynDiv_c stdc++fs)

# install function
install(TARGETS SynDiv_c DESTINATION bin)

# 安装头文件
install(DIRECTORY include/ DESTINATION include)

install(FILES SynDiv.py SynDiv_p.py genome2SynDiv_config.py cal_Syn_Fst.py gene_Syn_Fst.py
        DESTINATION bin
)

# copy
install(CODE "execute_process(COMMAND cp -r include ${CMAKE_INSTALL_PREFIX}/bin/)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/SynDiv.py)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/SynDiv_p.py)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/genome2SynDiv_config.py)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/cal_Syn_Fst.py)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/gene_Syn_Fst.py)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/SynDiv.py ${CMAKE_INSTALL_PREFIX}/bin/SynDiv)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/SynDiv_p.py ${CMAKE_INSTALL_PREFIX}/bin/SynDiv_p)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/genome2SynDiv_config.py ${CMAKE_INSTALL_PREFIX}/bin/genome2SynDiv_config)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/cal_Syn_Fst.py ${CMAKE_INSTALL_PREFIX}/bin/cal_Syn_Fst)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/gene_Syn_Fst.py ${CMAKE_INSTALL_PREFIX}/bin/gene_Syn_Fst)")