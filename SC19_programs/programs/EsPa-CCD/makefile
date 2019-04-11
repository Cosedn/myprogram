CCD.exe:ccd.o main.o
	mpicc -o CCD.exe ccd.o main.o
ccd.o:ccd.c ccd.h
	mpicc -c ccd.c
main.o:main.c ccd.h
	mpicc -c main.c

clean:
	rm *.o -f
