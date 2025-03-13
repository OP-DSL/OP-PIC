# OP-PIC SimPIC

In simpic.h and funcs_from_simpic.cpp most of the original functions from https://bitbucket.org/lecad-peg/simpic/src/master/ remains. These will be used mainly for loading the mesh and particles.

The algorithm contains the original code used in simpic which has Matrix computations written explicitly using c++, which cannot be easily converted to kernel calls.

The content in the omp and seq folders and any other file with suffix _op (simpic_op.cpp) should be auto generated once the source to source translator is developed.
For now, only the auto generated code works (sequential)

Therefore, OP-PIC expects the user to write only simpic.cpp and kernels.h files.


