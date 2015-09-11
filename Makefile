FC=mpiifort
FCFLAGS=-O0 -g -traceback -assume nobuffered_io

example: example.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o

halos.o: halos.f90
	${FC} ${FCFLAGS} -c $<

clean:
	rm example halos.o
