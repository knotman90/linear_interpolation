# linear_interpolation

To compile/run the tests:

1. `mkdir -p build  && cd build`
2. `cmake .. && make`
3. `./tests/bin/tests`

Two versions, non of them perform sanity check on input. Bad input results in undefined behaviour: 
include/linear_interpolation.cxx implements interpolation on cuboid 1/2/3D 
include/interpolate.cxx expose a single function that perform interpolation on N-Dimensional.

Test folder include tests for both version.
It tests all 1/2/3/4D cuboid interpolation on NPOLY random polynomials and  for each of them it tests NPOINTS points.
