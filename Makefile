FC=mpiifort
FCFLAGS=-O0 -g -traceback -assume nobuffered_io

examples: example example1d

example: example.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o

example1d: example1d.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o

halos.o: halos.f90
	${FC} ${FCFLAGS} -c $<

clean:
	rm example halos.o
