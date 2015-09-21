      program veldis
      !---------------------------------------------------------------------------------
      ! Written by Burak H.
      ! This program computes the band velocities as a function of
      ! energy
      ! Files required:
      !  1. prefix.a2Fsave : for band energies and velocities
      !  2. fil_info       : at(3,3) and bg(3,3)
      !  
      !---------------------------------------------------------------------------------

!$    use omp_lib

      implicit none
      !
      ! Approximation to Dirac delta
      interface
        function w0gauss (x)
           double precision :: w0gauss
           double precision, intent(in) :: x
        end function w0gauss
      end interface
      !
      !
      integer :: i,j,k,iq,ik,ie,itemp,nu,nqtot,nsig,nat,nk1fit,nk2fit,  &
     &           nk3fit, nkfit, nksfit_real, nbnd, nksfit, npk, nsym,   &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f, ibnd_ph, &
     &           nphband, n, counter, Nen
      !
      double precision :: stemp, wk, at(3,3), bg(3,3), efermi, alat,    &
     &                    sig(3,3), vol, fac, tau_inv, temp_tau,        &
     &                    temperature, fd, dfd, emin, emax, En(1000),   &
     &                    degauss, v2e, nef
      ! 
      !double precision, allocatable :: wq(:),xq(:,:), w2(:)
      !
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                                 etfull(:,:),vk(:,:,:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      integer :: nthreads
      !
      logical :: lsoc
      !
      character*20 :: fil_kp, fil_a2F, fil_info, fil_tau, temp
      !
      double precision, PARAMETER :: Rytocm1 = 109737.57990d0,          &
     &                               RytoGHz = 3.289828D6,              &
     &                               RytoTHz = RytoGHz/1000.d0,         &
     &                               RytoeV = 13.60569253d0,            &
     &                               tpi = 6.283185307d0,               &
     &                               convsig = 2.89d5,                  &
     &                               KtoRy = 1.d0/38.681648/300/RytoeV
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 

      namelist /input/ fil_info, fil_a2F, alat,                         &
     &                 nk1fit, nk2fit, nk3fit, phband_i,                &
     &                 degauss, phband_f, efermi, nthreads, lsoc

      read(5,input)

      npk = 40000
      !
      !$ call omp_set_num_threads(nthreads)

      nkfit=nk1fit*nk2fit*nk3fit
      if ( lsoc .eqv. .true.) then
         wk = 1.d0 / nkfit
      else
         wk = 2.d0 / nkfit
      end if
      !
      ! Convert efermi to Ryd
      efermi = efermi / RytoeV
      !
      ! Convert to Ry
      !
      En = En / RytoeV
      !
      ! Convert degauss to Ry
      !
      degauss = degauss / RytoeV
      !
      nphband = phband_f - phband_i + 1 ! Total number of bands of interest (usually the number of relevant conduction bands)
      !
      ! Read a2F
      !
      open(11,file=fil_a2F,status='unknown')
      !
      read(11,*) nbnd, nksfit
      !
      allocate(etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
      !allocate(xq(3,nqtot))
      !allocate(wq(nqtot))
      !         etfit : band energies in IBZ
      !         xkfit : k-point coordinates
      !         wkfit :  weights
      ! Read data
      read(11,*) etfit
      read(11,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
      read(11,*) wkfit
      read(11,*) nk1fit, nk2fit, nk3fit
      read(11,* ) nsym
      do ns=1,nsym
         read(11,*)  ((s(i,j,ns),j=1,3),i=1,3)
      enddo
      !
      close(11)
      !
      ! Read info file on k-points (and lattice)
      !
      open(11,file=fil_info,status='unknown')
      !
      read(11,*)
      read(11,*) ((at(i,j),i=1,3),j=1,3)
      !
      read(11,*)
      read(11,*)
      !
      read(11,*) ((bg(i,j),i=1,3),j=1,3)
      ! 
      close(11)
      !
      allocate(etfull(nbnd,nkfit),vk(nphband,nkfit,3))
      !
      allocate (eqkfit(nkfit), sfit(nkfit))
      !         eqkfit : band energies in uniform grid
      !         sfit   : symmetries
      ! 
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      ! 
      ! eqkfit(nk) : maps IBZ to full grid. The full grid is in crystal coords
      !
      !$omp parallel do default(none) &
      !$omp collapse(2) &
      !$omp private(ibnd,ik) shared(etfull,etfit,eqkfit,nbnd,nkfit)
      do ibnd=1,nbnd
         do ik=1,nkfit
            etfull(ibnd,ik) = etfit(ibnd,eqkfit(ik))
         end do
      end do
      !$omp end parallel do
      ! 
      ! Deallocate unnecessary variables 
      !
      deallocate (sfit, xkfit, wkfit, etfit)
      !
      ! call the band velocities
      !
      call vband( nk1fit,nk2fit,nk3fit,nphband,etfull(phband_i:phband_f,:),at,vk )
      !
      vk = vk / tpi * alat
      !
      ! Compute the energy dependent velocity squared
      !
      nef = 0.d0
      v2e = 0.d0
      !
      do ibnd=phband_i, phband_f  
         !
         ibnd_ph = ibnd - phband_i + 1
         !
         !$omp parallel do default(shared) &
         !$omp private(ik) &
         !$omp reduction(+ : nef, v2e) 
         do ik=1,nkfit
            !
            ! DOS
            nef = nef + wk * w0gauss((efermi-etfull(ibnd,ik))/degauss)/degauss
            !
            ! v2e
            v2e = v2e + wk * vk(ibnd_ph,ik,1)**2 *                      &
     &               w0gauss((efermi-etfull(ibnd,ik))/degauss)/degauss
            !
         end do ! ik
         !$omp end parallel do
         !
      end do ! ibnd 
      !
      v2e = v2e / nef * RytoeV**2.0
      !
      ! Write to file
      !
      open(11,file='v2e.dat',status='unknown')
      write(11,'(2e14.6)') nef / RytoeV, v2e / alat**2 
      close(11)
      !
      ! Clean memory
      !
      deallocate(eqkfit,etfull)
      !
      end program veldis

      !
      FUNCTION w0gauss (x)

      implicit none

      double precision :: w0gauss

      double precision, intent(in) :: x

      double precision :: sqrtpm1, arg

      sqrtpm1 = 1.0d0/1.77245385090551602729d0

      !
      arg = min(200.d0,x**2.d0)
      !arg = (x/sg)**2.d0
      !
      w0gauss = sqrtpm1 * exp( -arg ) 
      !
      END FUNCTION w0gauss
    
