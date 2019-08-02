CC=gcc
CFLAGS= -std=gnu11 -Wall -O3 -ggdb -lz -fopenmp
all:
	$(CC) $(CFLAGS)  *.c -o ./kssd -lm
cmpnt6:
	$(CC) $(CFLAGS) -DCOMPONENT_SZ=6 *.c -o ./kssd -lm
