main_iplib: bmp.o ip_lib.o
	gcc -Wall main_iplib.c bmp.o ip_lib.o -o main_iplib -lm --ansi --pedantic -g3 -O3 -std=gnu89 -fsanitize=address -fsanitize=undefined -Wextra
bmp: bmp.c bmp.h
	gcc -c bmp.c -c -Wall -lm

ip_lib: ip_lib.c ip_lib.h
	gcc -c ip_lib.c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
