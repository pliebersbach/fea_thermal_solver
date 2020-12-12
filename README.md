# fea_solver
Author: Piotr Liebersbach

C++ based finite element 2D thermal analysis code.  Analysis is done using 4 node quadrilateral isoparemetric elements.  All stiffness matrix and loads assembly are calculated using full Gauss Quadrature integration rules.  Linear system of equations is solved with Gauss elemination solver.

To build test cases move into test/ directory and build with:
  cmake .; make
  
Run test cases by running the executables, ie:
  ./plate

To plot the nodal solution run the provided Matlab/Octave script "plotter.m".  The script will generate a contour plot of the temperature profile over the mesh, and save the figure as "Temp.png".  All temperature contour plots for the test cases included in the test directory.
