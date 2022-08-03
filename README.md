# F90-MPAS
Augmentations of the occ package of M. Gorbunov to simulate RO soundings 
by wave optics integration using the output of the Model for Prediction 
Across Scales (MPAS). One new module is provided (MPAS) and two are 
modified to incorporate the additions (Wave, Wavelib). The last previous 
changes made to the base package were made in November, 2009. 

For this code to function correctly, the entire package must be built 
using a standard installation of NetCDF with Fortran 90 interfaces. The 
stripped down version of NetCDF 90 interfaces based on old Fortran 77 
subroutines provided by the old occ package will not function correctly 
with the MPAS augmentations. 

