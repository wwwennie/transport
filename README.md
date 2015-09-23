# Transport codes
## Description

This reposiroty contains transport codes based on the electron-LO phonon scattering.
Currently, these codes work for cubic crystals only. However, the generalization is straightforward.

There are several main programs:

* tauBZ_1.1.1_omp.f90 : Calculates electron-phonon interactions in the BZ
* transp.1.2_omp.f90  : Calculates transport coefficients as a function of Ef for constant scattering rate
* transp.2.1_omp.f90  : Calculates transport coefficients using the rates calculated from tauBZ_1.1.1_omp.f90
* transp_densities    : Calculates some densities that are useful for interpreting data
* veldis              : Calculates band velocity squared for as a function of Ef, useful for interpreting data

The required input files are:

* prefix.a2Fsave : for band energies and velocities
* fil_info       : at(3,3) and bg(3,3), i.e. direct and reciprocal lattice vectors
* Cnu.txt        : el-ph parameters
* wo.in          : LO frequencies

## Things to do

* Adding velocity dependencies to scattring rate calculations
* Unifying collinear and noncollinear calculations
* Unifying constant rate and full rate calculations
* MPI parallelization
