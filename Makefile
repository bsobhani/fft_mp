#CFLAGS =        -g -Wall -pedantic -Wstrict-prototypes
CFLAGS =       -O2 -Wall -ansi -pedantic -Wstrict-prototypes    

example: example.c fft.c
	gcc $(CFLAGS) example.c fft.c -o example -lm

clean:
	rm -f example
