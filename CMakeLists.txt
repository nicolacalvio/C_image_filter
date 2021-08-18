cmake_minimum_required(VERSION 3.15)
project(progetto_c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -Wextra -std=gnu89 ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -lm -g -O1 ")
add_executable(progetto_c main_iplib.c ip_lib.c bmp.c)
target_link_libraries(progetto_c m)