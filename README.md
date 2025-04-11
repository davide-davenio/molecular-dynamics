# molecular-dynamics
Molecular dynamics code in MatLab created for a course project. Course professors were Francesco Montalenti and Roberto Bergamaschini.

"MolecularDynamicsReport.pdf" contains a detailed description of the process
followed to create and test the code, plus a small example analysis of addatom
diffusion patterns on an infinite crystal slab.

Main simulation code is MolecularDynamics.m, while Minimization.m is a kinetik
steepest descent implementation to be used on initial configuations if 
periodic boundary conditions are not present.

Subroutines folder contains all secondary functions. Tests folder contains debugging scripts used to test the secondary functions.

Configurations folder contains premade crystal configurations used in testing.

Analysis folder contains the scripts used to visualize data along the building
process. Raw data folder contains this data in its original form. Images and movies
folders contains the visualizations. 

legacy_code folder contains older versions of the simulation code and secondary functions for backlogging purposes.
