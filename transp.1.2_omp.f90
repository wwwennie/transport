      program transp
      !---------------------------------------------------------------------------------
      ! Written by Burak H. v1.2
      ! Computes conductivity, Seebeck coefficient, Hall coefficient
      ! (for cubic systems only)
      ! This program will be improved to compute other transport coefficients as well.
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
      integer :: i,j,iq,ik,imu,itemp,nu,nqtot,nsig,nat,nk1fit,nk2fit,   &
     &           nk3fit, nkfit, nksfit_real, nbnd, nksfit, npk, nsym,   &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f, ibnd_ph, &
     &           nphband, Nmu 
      !
      double precision :: stemp, wk, at(3,3), bg(3,3), Ef(1000), alat,  &
     &                    sig(1000), vol, fac, tau_inv, temp_tau,       &
     &                    temperature, fd, dfd, efmin, efmax, nelec1,   &
     &                    cbm, vbm, ef_mid, nelec2,                     &
     &                    Se(1000),l12(1000),nef(1000),sigRH(1000),     &
     &                    RH(1000)
      ! 
      !double precision, allocatable :: wq(:),xq(:,:), w2(:)
      !
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                               etfull(:,:),vk(:,:,:),mk(:,:,:,:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      ! OMP
      double precision :: t0
      !
      character*20 :: fil_kp, fil_a2F, fil_info, fil_sig,temp
      !
      double precision, PARAMETER :: Rytocm1 = 109737.57990d0,          &
     &                               RytoGHz = 3.289828D6,              &
     &                               RytoTHz = RytoGHz/1000.d0,         &
     &                               RytoeV = 13.60569253d0,            &
     &                               tpi = 6.283185307d0,               &
     &                               convsig = 2.89d5,                  &
     &                               KtoRy = 1.d0/38.681648/300/RytoeV, &
     &                               autocm = 5.2917721092d-9,          &
     &                               convRH = 9.2522431d-13
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 

      namelist /input/ fil_info, fil_a2F, nat, alat, vol,               &
     &                 nk1fit, nk2fit, nk3fit, tau_inv, phband_i,       &
     &                 temperature, phband_f, Nmu, efmin, efmax, vbm,   &
     &                 cbm

      read(5,input)

      npk = 40000
      !
      ! Chemical potentials
      !
      do imu=1,Nmu
        Ef(imu) = (efmax-efmin)/(Nmu-1) * (imu-1) + efmin
      end do
      ! 

      nkfit=nk1fit*nk2fit*nk3fit
      wk = 2.d0 / nkfit
      !
      ! Convert scattering rate from THz to Ryd  
      tau_inv = tau_inv / RytoTHz
      !
      ! Convert temperature from K to Ryd
      temperature = temperature * KtoRy 
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
      !
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
      allocate(etfull(nbnd,nkfit))
      allocate( vk(nphband,nkfit,3), mk(nphband,nkfit,3,3) ) ! Allocate only nphband number of bands 
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
      ! Check number of electrons, between IBZ and full-grid
      !
      ef_mid = (cbm + vbm)/2.0/RytoeV ! Fermi level at the middle of the gap
      !
      nelec1 = 0.d0
      nelec2 = 0.d0
      !  Irreducable BZ
      !$omp parallel do default(shared) &
      !$omp collapse(2) &
      !$omp private(ik,ibnd,fd) &
      !$omp reduction(+ : nelec1)
      do ik=1,nksfit
         do ibnd=1,nbnd  
            fd = 1.0/( exp( (etfit(ibnd,ik)- ef_mid)/temperature) +1.0)
            nelec1 = nelec1 + wkfit(ik)*fd
         end do
      end do
      !$omp end parallel do
      !
      !  Full BZ
      !$omp parallel do default(shared) &
      !$omp collapse(2) &
      !$omp private(ik,ibnd,fd) &
      !$omp reduction(+ : nelec2)
      do ik=1,nkfit
         do ibnd=1,nbnd  
            fd = 1.0/( exp( (etfull(ibnd,ik)- ef_mid)/temperature) +1.0)
            nelec2 = nelec2 + wk*fd
         end do
      end do
      !$omp end parallel do
      !
      if (abs(nelec1 - nelec2) .ge. 1.d-8) then
          write(6,*) 'Something wrong!'
          write(6,*) 'nelec1-nelec2= ', nelec1-nelec2
      end if
      ! 
      write(6,*) 'nelec=', nelec2
      !
      ! Deallocate unnecessary variables 
      !
      deallocate (sfit, xkfit, wkfit, etfit)
      !
      ! call the band velocities (only the required bands)
      !
      call vband( nk1fit,nk2fit,nk3fit,nphband,etfull(phband_i:phband_f,:),at,vk )
      !
      ! call the inverse mass tensor (only the required bands)
      !
      call mband ( nk1fit,nk2fit,nk3fit,nphband,vk,at,mk )
      !
      ! Divide by the reciprocal lattice size 2pi/a
      vk = vk / tpi * alat
      mk = mk / (tpi / alat)**2.0
      !
      ! Compute the conductivity tensor
      !
      sig = 0.d0
      l12 = 0.d0
      sigRH = 0.d0
      !
      ! Start timing here
      t0 = omp_get_wtime()
      !
      do imu=1,Nmu
         ! 
         do ibnd=phband_i, phband_f
            !
            ibnd_ph = ibnd - phband_i + 1
            !
            !$omp parallel do default(shared) &
            !$omp private(ik,fd,dfd,fac) &
            !$omp reduction(+ : nef,sig,l12,sigRH)
            do ik=1,nkfit
               !
               fd = 1.0/( exp( (etfull(ibnd,ik)- Ef(imu))/temperature) + 1.0)
               dfd = 1.0/temperature * fd * (1.0 - fd)
               fac = etfull(ibnd,ik)- Ef(imu)
               !
               ! N(ef)
               nef(imu) = nef(imu) + fd * wk
               !  
               !dfd = w0gauss((etfull(ibnd,ik)-Ef(imu))/temperature)/temperature
               ! 
               sig(imu) = sig(imu) + wk * vk(ibnd_ph,ik,1)**2.0*dfd/tau_inv
               !
               l12(imu) = l12(imu) + wk * vk(ibnd_ph,ik,1)**2.0*dfd/       &
     &                   tau_inv * fac 
               !
               sigRH(imu) = sigRH(imu) + wk * dfd *                     &
     &             ( vk(ibnd,ik,1)*vk(ibnd_ph,ik,1)*mk(ibnd_ph,ik,2,2) -      &
     &               vk(ibnd,ik,1)*vk(ibnd_ph,ik,2)*mk(ibnd_ph,ik,1,2) ) *    &
     &               (1.0/tau_inv)**2.d0
               ! 
            end do ! ik  
            !$omp end parallel do
         end do ! ibnd  
      end do ! imu
      !
      ! Measure OMP time here
      t0 = omp_get_wtime() - t0
      ! Seebeck 
      Se = l12/sig
      !
      ! Convert to V/K
      Se = - Rytoev/(temperature / KtoRy) * Se
      !
      ! Hall coefficient
      RH = -(sigRH/sig**2.0)*vol
      !
      ! Convert sig into (Ohm cm)-1
      sig = sig * convsig / vol
      !
      ! Convert RH into (m^3/Ohm)
      RH = RH * convRH
      !
      ! Convert n(ef) to cm-3
      nef = nef / (vol*autocm**3.0)  
     
      ! Print the conductivity tensor on file
      !
      open(11,file='sig.out',status='unknown')
      !
      do imu=1,Nmu
         write(11,"(5e14.6)") Ef(imu),sig(imu),Se(imu),RH(imu),nef(imu)
      end do
      !
      close(11)
      !
      write(6,*) 'K-point parallelized integration time= ', t0     
      ! Clean memory
      !
      deallocate(eqkfit,etfull,vk,mk)
      !
      end program transp

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
    
