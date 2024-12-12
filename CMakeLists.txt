cmake_minimum_required(VERSION 3.12)

project(SynDiv_c CXX)

# Set C++ standard version (using C++17 here)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add source files
set(SRC_FILES_SynDiv_c
    src/get_time.cpp
    src/save.cpp
    src/strip_split_join.cpp
    src/multiinter.cpp
    src/no_syn.cpp
    src/coor.cpp
    src/cal.cpp
    src/window.cpp
    src/Fst.cpp
)

# Add executable
add_executable(SynDiv_c SynDiv_c.cpp)

# Add source files
target_sources(SynDiv_c PRIVATE ${SRC_FILES_SynDiv_c})

# Add compile options for the target
target_compile_options(SynDiv_c PRIVATE -O3)

# Add compile options for the target (including -O3, additional warnings and errors, and CPU optimizations)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # target_compile_options(SynDiv_c PRIVATE -O3 -Wall -Wextra -Werror -pedantic-errors -march=native)
    # target_compile_options(SynDiv_c PRIVATE -O3 -Wall -Wextra -march=native)
    target_compile_options(SynDiv_c PRIVATE -O3 -Wall -Wextra)
endif()

# Add include directory
target_include_directories(SynDiv_c PRIVATE include)

# Add link libraries
find_package(Threads REQUIRED)
target_link_libraries(SynDiv_c Threads::Threads)

find_package(ZLIB REQUIRED)
target_link_libraries(SynDiv_c ZLIB::ZLIB)

target_link_libraries(SynDiv_c stdc++fs)

# install function
install(TARGETS SynDiv_c DESTINATION bin)

# Install header files
install(DIRECTORY include/ DESTINATION include)

install(FILES SynDiv.py SynDiv_p.py genome2SynDiv_config.py
        DESTINATION bin
)

# copy
install(CODE "execute_process(COMMAND cp -r include ${CMAKE_INSTALL_PREFIX}/bin/)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/SynDiv.py)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/SynDiv_p.py)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/genome2SynDiv_config.py)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/SynDiv.py ${CMAKE_INSTALL_PREFIX}/bin/SynDiv)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/SynDiv_p.py ${CMAKE_INSTALL_PREFIX}/bin/SynDiv_p)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/genome2SynDiv_config.py ${CMAKE_INSTALL_PREFIX}/bin/genome2SynDiv_config)")