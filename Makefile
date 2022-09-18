all: main.c
	gcc -Wall -O3 main.c -o repseq -static
