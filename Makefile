brego_1d: brego_1d.c
	gcc brego_1d.c -o brego_1d -Wall -lm

all: brego_1d

clean:
	rm brego_1d
