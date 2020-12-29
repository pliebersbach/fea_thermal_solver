# fea_solver
Author: Piotr Liebersbach

C++ based finite element 2D thermal analysis code.  Analysis is done using 4 node quadrilateral isoparemetric elements.  All stiffness matrix and loads assembly are calculated using full Gauss Quadrature integration rules.  Linear system of equations is solved with Gauss elemination solver.

To build the test cases and custom analysis build the cmakelists file in the main directory and compile with make:

  cmake {path_to_project_dir}; make
  
Run test cases by running the executables, ie:
  
  ./plate

To plot the nodal solution run the provided Matlab/Octave script "plotter.m".  The script will generate a contour plot of the temperature profile over the mesh, and save the figure as "Temp.png".  All temperature contour plots for the test cases included in the test directory.

The read_mesh() function reads a custom mesh format.  Meshes in this format can be generated with the Commercial/Free meshing software Coreform Cubit.  To export the mesh from the program into the correct file format use the provided cubit_mesh_writer.py .