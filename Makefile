CC=gcc
CFLAGS= -std=gnu11 -Wall -Wno-unused-result -O3 -ggdb -lz -fopenmp
all:
	$(CC) $(CFLAGS)  *.c -o ./kssd -lm
alert:
	$(CC) $(CFLAGS) -DCOMPONENT_SZ=6 *.c -o ./kssd -lm
strange:
	$(CC) $(CFLAGS) -DCTX_SPC_USE_L=10 *.c -o ./kssd -lm
