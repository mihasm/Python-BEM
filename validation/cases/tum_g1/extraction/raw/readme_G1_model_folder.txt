Structure of G1 WInd Turbine Model data folder
Contains the following 5 subfolders:
1) 3D CAD geometry: The file G1_Model_for_CFD_pitch_0d4.stp is a CAD model of the turbine including the tower and base
2) Airfoil profile: The file RG14 airfoil.txt contains the coordinates of the main airfoil profile   
3) Blade distributions: Contains the description of the blade planform geometry and the correspondence with the polar files.
The Excel file Blade_design_G1_RG14.xlsx contains (41 sections along the blade):
-	chord, twist and thickness distribution 
-	location of the pitch axis and the aerodynamic center  
-	geometry profiles of the blade sections 
The ascii file AeroDyn_blade.dat contains the chord / twist distribution and the corresponding name of the polar file (see "Polars" subfolder below) for a finer blade discretization (60 sections)
4) Polars: Contains the CL, CD, CM polars of the different airfoils for 4 Reynolds numbers (0.047, 0.057, 0.069, 0.084 million)  
section_1.dat  Polars for the cylinder profile
section_2.dat .... section_9.dat  Polars for the airfoils in the transition region (between the cylinder and the RG14 airfoil) 
section_10.dat  Polars for the RG14 airfoil
5) G1 OpenFast Model: Contains an OpenFast model of the turbine with all relevant data (polars, distributions etc).
   
