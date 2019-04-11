DSGD.exe:dsgd.o main.o
	mpicc -o DSGD.exe dsgd.o main.o
dsgd.o:dsgd.c dsgd.h
	mpicc -c dsgd.c
main.o:main.c dsgd.h
	mpicc -c main.c

clean:
	rm *.o -f
