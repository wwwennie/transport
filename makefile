#
#
F90 = gfortran
FFLAGS = -fopenmp

tau_BZ.1.1.1_omp: tau_BZ.1.1.1_omp.f90
	${F90} ${FFLAGS} -o tau_BZ.1.1.1_omp.x tau_BZ.1.1.1_omp.f90 cryst_to_car.f90 lint.f90 find_kpq.f90

tau_BZ.1.1.1: tau_BZ.1.1.1.f90
	${F90} -o tau_BZ.1.1.1.x tau_BZ.1.1.1.f90 cryst_to_car.f90 lint.f90 find_kpq.f90

transp.1.1_omp: transp.1.1_omp.f90
	${F90} ${FFLAGS} -o transp.1.1_omp.x transp.1.1_omp.f90 cryst_to_car.f90 lint.f90 vband.f90 mband.f90

transp.1.1: transp.1.1.f90
	${F90} -o transp.1.1.x transp.1.1.f90 cryst_to_car.f90 lint.f90 vband.f90 mband.f90

transp.1.2_omp: transp.1.2_omp.f90
	${F90} ${FFLAGS} -o transp.1.2_omp.x transp.1.2_omp.f90 cryst_to_car.f90 lint.f90 vband.f90 mband.f90

transp.2.1_omp: transp.2.1_omp.f90
	${F90} ${FFLAGS} -o transp.2.1_omp.x transp.2.1_omp.f90 cryst_to_car.f90 lint.f90 vband.f90 mband.f90

transp.2.1: transp.2.1.f90
	${F90} -o transp.2.1.x transp.2.1.f90 cryst_to_car.f90 lint.f90 vband.f90 mband.f90

transp_map.0.1_omp: transp_map.0.1_omp.f90
	${F90} ${FFLAGS} -o transp_map.0.1_omp.x transp_map.0.1_omp.f90

transp_map.0.2_omp: transp_map.0.2_omp.f90
	${F90} ${FFLAGS} -o transp_map.0.2_omp.x transp_map.0.2_omp.f90

transp_densities: transp_densities.f90
	${F90} ${FFLAGS} -o transp_densities.x transp_densities.f90 cryst_to_car.f90 lint.f90 vband.f90

veldis: veldis.f90
	${F90} ${FFLAGS} -o veldis.x veldis.f90 cryst_to_car.f90 lint.f90 vband.f90
