.SUFFIXES: .c .u
CC= gcc
# CFLAGS_MAC = -g -Wall -O3 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -Winline -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -malign-natural -ffast-math -mdynamic-no-pic -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing -mcpu=7450 -faltivec
CFLAGS_MAC = -g -Wall -O3 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -Winline
CFLAGS_PTON = -g -Wall -O3 -DHAVE_INLINE=1 -DGSL_RANGE_CHECK_OFF=1
CFLAGS_DEBUG = -g -Wall
CFLAGS = -g -Wall -I/opt/local/include/gsl/ -I/usr/include/sys/ -I/usr/include/  -I/usr/local/include -L/usr/local/lib -I/usr/local/include -L/usr/local/lib -lm -lgsl -lgslcblas -I/opt/local/include/gsl

# MAC_LDFLAGS = -lgsl -latlas -lcblas -L/sw/li
MAC_LDFLAGS = -lgsl -lgslcblas -L/opt/local/lib  -I/opt/local/include/gsl/ -I/usr/include/sys/ -I/usr/include/  -I/usr/local/include -L/usr/local/lib -I/usr/local/include -L/usr/local/lib -lm   -lgsl -lgslcblas -I/opt/local/include/gsl
C2_LDFLAGS = -lgsl -lcblas -latlas
CYCLES_LDFLAGS = -lgsl -lgslcblas
LSOURCE = utils.c topic.c doc.c hyperparameter.c main.c gibbs.c
LOBJECTS = utils.o topic.o doc.o hyperparameter.o main.o gibbs.o

main:	$(LOBJECTS)
	$(CC) $(CFLAGS_MAC) $(LOBJECTS) -o hlda $(MAC_LDFLAGS)

c2:	$(LOBJECTS)
	$(CC) $(CFLAGS_PTON) $(LOBJECTS) -o hlda $(C2_LDFLAGS)

cycles:	$(LOBJECTS)
	$(CC) $(CFLAGS_PTON) $(LOBJECTS) -o hlda $(CYCLES_LDFLAGS)

debug:	$(LOBJECTS)
	$(CC) $(CFLAGS_DEBUG) $(LOBJECTS) -o hlda $(MAC_LDFLAGS)

clean:
	-rm -f *.o
