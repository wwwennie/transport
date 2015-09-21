      program transp_densities
      !---------------------------------------------------------------------------------
      ! Written by Burak H.
      ! Computes certain densities useful for transport
      ! Files required:
      !  1. prefix.a2Fsave : for band energies and velocities
      !  2. fil_info       : at(3,3) and bg(3,3)
      !  3. tauBZ.dat      : scattering rates in IBZ
      !  
      !---------------------------------------------------------------------------------

      !$ use omp_lib

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
      logical :: lsoc
      !
      integer :: i,j,k,iq,ik,imu,itemp,nu,nqtot,nsig,nat,nk1fit,nk2fit, &
     &           nk3fit, nkfit, nksfit_real, nbnd, nksfit, npk, nsym,   &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f, ibnd_ph, &
     &           nphband, n, counter
      !
      double precision :: stemp, wk, at(3,3), bg(3,3), efermi, alat,    &
     &                    sig(3,3), vol, fac, tau_inv, temp_tau,        &
     &                    temperature, fd, dfd, efmin, efmax,           &
     &                    degauss, wg0
      ! 
      !double precision, allocatable :: wq(:),xq(:,:), w2(:)
      !
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                                 etfull(:,:),vk(:,:,:),tauk(:,:), &
     &                                 tauk_full(:,:), DSig(:), Dv2(:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      integer :: nthreads
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

      namelist /input/ fil_info, fil_a2F, fil_tau, nat, alat, vol,      &
     &                 nk1fit, nk2fit, nk3fit, phband_i,                &
     &                 degauss, phband_f, efermi, lsoc, nthreads

      read(5,input)

      ! Set the number of threads
      !$ call omp_set_num_threads(nthreads)

      npk = 40000
      !

      nkfit=nk1fit*nk2fit*nk3fit
      !
      if ( lsoc .eqv. .true.) then
         wk = 1.d0 / nkfit
      else
         wk = 2.d0 / nkfit
      end if
      !
      ! Convert efermi to Ryd
      efermi = efermi / RytoeV
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
      ! Allocate and read tau-1 from file (IBZ)
      !
      allocate(tauk(nksfit,nphband))
      !
      tauk = 0.d0
      !
      open(11,file=fil_tau,status='unknown')
      !
      do ik=1,nksfit
         read(11,*) (tauk(ik,ibnd),ibnd=1,nphband)
      end do
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
      allocate( etfull(nbnd,nkfit) )
      allocate( vk(nphband,nkfit,3) )    ! Allocate only needed number of bands
      allocate(tauk_full(nkfit,nphband))
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
      do ibnd=1,nbnd
         do ik=1,nkfit
            etfull(ibnd,ik) = etfit(ibnd,eqkfit(ik))
         end do
      end do
      ! 
      do ibnd=1,nphband
         do ik=1,nkfit
            tauk_full(ik,ibnd) = tauk(eqkfit(ik),ibnd)
         end do
      end do
      ! 
      ! Deallocate unnecessary variables, allocate new ones 
      !
      deallocate (sfit, xkfit, wkfit, etfit, tauk, eqkfit)
      !
      allocate (DSig(nphband), Dv2(nphband))

      ! call the band velocities
      !
      call vband( nk1fit,nk2fit,nk3fit,nphband,etfull(phband_i:phband_f,:),at,vk )
      !
      vk = vk / tpi * alat
      !
      ! Compute the densities to be plotted 
      !
      DSig = 0.d0
      Dv2 = 0.d0
      !
      do ibnd=phband_i, phband_f
         !
         ibnd_ph = ibnd - phband_i + 1
         !
         !$omp parallel do default(shared) &
         !$omp private(ik,wg0) &
         !$omp reduction(+ : DSig, Dv2)
         do ik=1,nkfit
            !
            wg0 = w0gauss( (etfull(ibnd,ik)-efermi)/degauss ) / degauss
            !
            DSig(ibnd_ph) = DSig(ibnd_ph) + vk(ibnd_ph,ik,1)**2 * wg0      &
     &                      * 1.0 / tauk_full(ik,ibnd_ph) * wk 
            !
            Dv2(ibnd_ph) = Dv2(ibnd_ph) + vk(ibnd_ph,ik,1)**2 * wg0 * wk
            ! 
         end do
         !$omp end parallel do
         !
      end do
      !
      ! Print the densities
      open(11,file='DSig.dat',status='unknown')
      open(12,file='Dv2.dat',status='unknown')
      write(11,*) (DSig(ibnd),ibnd=1,nphband) 
      write(12,*) (Dv2(ibnd)*RytoeV,ibnd=1,nphband) 
      close(11)
      close(12)
      !
      ! Clean memory
      !
      deallocate(etfull,tauk_full,vk,DSig,Dv2)
      !
      end program transp_densities

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
    
