CC = cc
CFLAGS = -O3
FTN = gfortran
FFLAGS = -O2
DEFC =
DEFF =
PRECISION =

all: test1.Abs test2.Abs test3.Abs

test1.Abs:	coefficients.h f90_interfaces.inc test.inc zsetup.c timer.c kernel1.c test1.F90
	$(CC) $(CFLAGS) $(DEFC) $(PRECISION) -c zsetup.c timer.c kernel1.c
	$(FTN) $(FFLAGS) $(DEFF) $(PRECISION) test1.F90 zsetup.o timer.o kernel1.o -o test1.Abs

test2.Abs:	coefficients.h f90_interfaces.inc test.inc zsetup.c timer.c kernel1.c test1.F90
	$(CC) $(CFLAGS) $(DEFC) $(PRECISION) -c zsetup.c timer.c kernel2.c
	$(FTN) $(FFLAGS) $(DEFF) $(PRECISION) test2.F90 zsetup.o timer.o kernel2.o -o test2.Abs

test3.Abs:	coefficients.h f90_interfaces.inc test.inc zsetup.c timer.c kernel1.c test1.F90
	$(CC) $(CFLAGS) $(DEFC) $(PRECISION) -c zsetup.c timer.c kernel3.c
	$(FTN) $(FFLAGS) $(DEFF) $(PRECISION) test3.F90 zsetup.o timer.o kernel3.o -o test3.Abs

clean:
	rm -f *.Abs *.o
