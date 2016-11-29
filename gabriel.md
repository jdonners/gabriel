project:   Gabriel
Summary: A Fortran library for convenient MPI
Author: John Donners
Date:    February 16, 2016
blank-value:
project_github: https://github.com/jdonners/gabriel
project_dir: .
src_dir: .
output_dir: ./doc
docmark: <
docmark_alt: !
exclude: bandwidth.f90

Gabriel is a Fortran library to simplify the use of MPI.
It offers verification, ease of use and a high performance, 
initially for models on regular grids.

The main innovations of gabriel are: 

- the user only needs to define the composition of the arrays from all ranks. Each rank provides the region of the array that is computed locally and  how that fits in the global domain. 
- the user then asks for a distribution for halo exchanges or a transformation from one composition to another. The library automatically deduces neighboring ranks and the subarrays to communicate.
- the user supplies array indices as defined in his application.

This is free and unencumbered software released into the public domain.
