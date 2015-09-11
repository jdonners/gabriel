FC=mpiifort
FCFLAGS=-g -traceback

example: example.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o

halos.o: halos.f90
	${FC} ${FCFLAGS} -c $<

