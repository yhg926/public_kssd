CC=gcc
CFLAGS= -std=gnu11 -Wall -Wno-unused-result -O3 -ggdb -lz -fopenmp
all:
	$(CC) $(CFLAGS)  *.c -o ./kssd -lm
alert:
	$(CC) $(CFLAGS) -DCOMPONENT_SZ=8 *.c -o ./kssd_CSZ8 -lm
strange:
	$(CC) $(CFLAGS) -DCTX_SPC_USE_L=10 *.c -o ./kssd_strange -lm
16S:
	$(CC) $(CFLAGS) -DMIN_KM_S=1 *.c -o ./kssd16S -lm

