FC=mpiifort
FCFLAGS=-O0 -g -traceback -stand f15
EXAMPLES=example example1d example_auto

examples: ${EXAMPLES}

example: example.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o

example1d: example1d.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o

example_auto: example_auto.f90 halos.o
	${FC} ${FCFLAGS} -o $@ $< halos.o


halos.o: halos.f90
	${FC} ${FCFLAGS} -c $<

clean:
	rm ${EXAMPLES} halos.o
