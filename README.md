# Transport codes
## Description

This reposiroty contains transport codes based on the electron-LO phonon scattering.
Currently, these codes work for cubic crystals only.
There are several main programs:

* tauBZ_1.1.1_omp.f90 : Calculates electron-phonon interactions in the BZ
* transp.1.2_omp.f90  : Calculates transport coefficients as a function of Ef for constant scattering rate
* transo.2.1_omp.f90  : Calculates transport coefficients using the rates calculated from tauBZ_1.1.1_omp.f90
* transp_densities    : Calculates some densities that are useful for interpreting data
* veldis              : Calculates band velocity squared for as a function of Ef, useful for interpreting data
