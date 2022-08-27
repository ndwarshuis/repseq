all: main.c
	gcc -Wall -O2 main.c -o repseq -lm
