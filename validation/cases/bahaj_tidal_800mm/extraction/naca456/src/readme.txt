

NACA456 - COORDINATES OF NACA AIRFOILS          /naca456/readme.txt

The files for this program are in the directory naca456 
  naca456.f90  the main program source code
  nacax.f90    the source code to the module of auxiliary code
  splprocs.f90 the source code for the spline procedures
  epspsi.f90   the data module defining the epsilon and psi functions
  avd.f90      the source code for the program that makes web pages
                 in HTML with the data from Abbott and von Doenhoff
  input.txt    instructions for preparing input for naca456
  ver6.zip     the program,source code, instructions and examples as
                 distributed in versions 1-6 of the PDAS CD-ROM
  samples.zip    a selection of test cases (input and output)
  samplnx.zip    all of the above test cases with Unix end-of-line (ZIP)

The reference documents for this program may be accessed
from the web page https://www.pdas.com/naca456refs.html. 

To compile this program for your computer, use the command
   gfortran  -fdefault-real-8 naca456.f90 -o naca456.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.
All of the auxiliary modules will be included automatically.
It is best if you use the option that makes REAL variables 64-bit.
Using 32-bit arithmetic can give values that are off by one in
the sixth significant digit. 

This program is a complete revision of the NASA Langley programs for computing
the coordinates of NACA airfoils. The computational procedure is described in
AIAA 2001-5235 and also at the PDAS web site at 
https://www.pdas.com/naca456refs.html
The computational procedures are highly modularized and I hope this will prove
useful in other computing projects. All of the source code is public domain
and open source. Many advanced features of Fortran 95 are employed
to make the code compact and efficient. The source code is abundantly 
commented in the hope that an auxiliary document is not needed. Parametric
cubic splines are used extensively in the calculations and these procedures
are collected in a module called SplineProcdures. The epsilon and psi 
functions used in the definition of the 6-series airfoils are defined in
a module called EpsilonPsi. The remainder of the airfoil profile and mean
line calculations are in a module called NacaAuxiliary. The main program
in naca456.f90 is, therefore, quite compact.

This program simply asks for the name of the input file. This must be
a file written to conform to the format described in the file naca456.txt
which uses namelist input. The file input.txt describes the input
variables.

The program produces a file called naca.out that printed or scrolled
to your screen and a file called naca.gnu that may be plotted.
Using gnuplot, you plot the airfoil shape with
   gnuplot> plot 'naca.gnu' with lines

The default setting for gnuplot will expand the vertical scale to fit
the screen. This is useful for seeing the curvature, etc. of the airfoil
and I usually take a look at most airfoils with this expanded scale.
But you probably want to see the airfoil in true proportions. Say
   gnuplot> set size ratio -1
and you will get a short but wide window with the correct proportions.
Alternately, you may command
   gnuplot> set yrange [-0.5:0.5]
   gnuplot> set size square

A large number of additional sample cases is included.
