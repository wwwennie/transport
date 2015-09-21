      program transp
      !---------------------------------------------------------------------------------
      ! Written by Burak H. (v2.1)
      ! Computes conductivity tensor,Seebeck coefficient and Hall coefficient.
      ! This program will be improved to compute other transport coefficients as well.
      ! Files required:
      !  1. prefix.a2Fsave : for band energies and velocities
      !  2. fil_info       : at(3,3) and bg(3,3)
      !  3. tauBZ.dat      : scattering rates in IBZ
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
     &           nphband
      !
      double precision :: stemp, wk, at(3,3), bg(3,3), efermi, alat,    &
     &                    sig,sigRH, vol, fac, tau_inv, temp_tau,       &
     &                    temperature, fd, dfd, efmin, efmax, nelec1,   &
     &                    cbm, vbm, ef_mid, nelec2, Se, l12, RH
      ! 
      !double precision, allocatable :: wq(:),xq(:,:), w2(:)
      !
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                                 etfull(:,:),vk(:,:,:),tauk(:,:), &
     &                                 tauk_full(:,:),mk(:,:,:,:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      double precision :: t0
      !
      character*20 :: fil_kp, fil_a2F, fil_info, fil_tau, temp
      !
      double precision, PARAMETER :: Rytocm1 = 109737.57990d0,          &
     &                               RytoGHz = 3.289828D6,              &
     &                               RytoTHz = RytoGHz/1000.d0,         &
     &                               RytoeV = 13.60569253d0,            &
     &                               tpi = 6.283185307d0,               &
     &                               convsig = 2.89d5,                  &
     &                               KtoRy = 1.d0/38.681648/300/RytoeV, &
                                     convRH = 0.092522431d-14
                                     ! convsig: conversion factor for 
                                     ! conductivity. From au to (Ohm cm)-1 
                                     ! convRH: From au to m^3/C

      namelist /input/ fil_info, fil_a2F, fil_tau, nat, alat, vol,      &
     &                 nk1fit, nk2fit, nk3fit, phband_i, tau_inv,       &
     &                 temperature, phband_f, efermi, vbm, cbm

      read(5,input)

      npk = 40000
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
      allocate(etfull(nbnd,nkfit),vk(nbnd,nkfit,3),mk(nbnd,nkfit,3,3))
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
      ! Maps
      !$omp parallel do default(shared) &
      !$omp collapse(2) &
      !$omp private(ik,ibnd) 
      do ibnd=1,nbnd
         do ik=1,nkfit
            etfull(ibnd,ik) = etfit(ibnd,eqkfit(ik))
         end do
      end do
      !$omp end parallel do
      ! 
      !
      !$omp parallel do default(shared) &
      !$omp collapse(2) &
      !$omp private(ik,ibnd)
      do ibnd=1,nphband
         do ik=1,nkfit
            tauk_full(ik,ibnd) = tauk(eqkfit(ik),ibnd)
         end do
      end do
      !$omp end parallel do 
      !
      ! Check number of electrons
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
      !
      if (abs(nelec1 - nelec2) .ge. 1.d-8) then
          write(6,*) 'Something wrong!'
          write(6,*) 'nelec1 - nelec2 = ', nelec1 - nelec2
      end if
      ! 
      write(6,*) 'nelec=', nelec2
      !
      ! Deallocate unnecessary variables 
      !
      deallocate (sfit, xkfit, wkfit, etfit, tauk)
      !
      ! call the band velocities
      !
      call vband( nk1fit,nk2fit,nk3fit,nbnd,etfull,at,vk )
      !
      ! call the inverse mass tensor
      !
      call mband ( nk1fit,nk2fit,nk3fit,nbnd,vk,at,mk )  
      !
      ! Divide by the reciprocal lattice size 2pi/a
      !
      vk = vk / tpi * alat
      mk = mk / (tpi * alat)**2.0
      !
      ! Compute the conductivity, Seebeck, and RH (cubic crystals only!)
      !
      sig   = 0.d0
      l12   = 0.d0
      sigRH = 0.d0
      !
      ! Start timing here
      t0 = omp_get_wtime()
      !
      do ibnd=phband_i, phband_f
         !
         ibnd_ph = ibnd - phband_i + 1
         !
         !$omp parallel do default(shared) &
         !$omp private(ik,fd,dfd,fac) &
         !$omp reduction(+ : sig,l12,sigRH)
         do ik=1,nkfit
            !
            fd = 1.0/( exp( (etfull(ibnd,ik)- efermi)/temperature) + 1.0)
            dfd = 1.0/temperature * fd * (1.0 - fd)
            fac = etfull(ibnd,ik)- efermi
            !
            !dfd = w0gauss((etfull(ibnd,ik)-Ef(imu))/temperature)/temperature
            !
            sig = sig + wk*vk(ibnd,ik,1)*vk(ibnd,ik,1) * dfd *          &
     &                  (1.0/tauk_full(ik,ibnd_ph))
            !
            l12 = l12 + wk*vk(ibnd,ik,1)*vk(ibnd,ik,1) * dfd *          &
     &                  (1.0/tauk_full(ik,ibnd_ph)) * fac
            !
            sigRH = sigRH + wk*dfd*                                     &
     &         ( vk(ibnd,ik,1)*vk(ibnd,ik,1)*mk(ibnd,ik,1,1) -          &
     &           vk(ibnd,ik,1)*vk(ibnd,ik,2)*mk(ibnd,ik,1,2) ) *        &
     &         (1.0/tauk_full(ik,ibnd_ph))**2.d0  
            !
         end do ! ik  
         !$omp end parallel do
      end do ! ibnd  
      !
      ! Measure OMP time here
      t0 = omp_get_wtime() - t0
      !
      ! Seebeck
      Se = l12/sig
      !
      ! Convert to V/K
      Se = -RytoeV/(temperature/KtoRy)*Se
      !
      ! Hall coefficient
      RH = -(sigRH/sig**2.0)*vol
      !
      ! Convert sig into (Ohm cm)-1
      sig = sig * convsig / vol 
      !
      ! Convert RH into (m^3/Ohm)
      sigRH = sigRH * convRH
      !
      ! Print the conductivity tensor on file
      !
      open(11,file='sig.out',status='unknown')
      !
      write(11,"(4e14.6)") efermi, sig, Se, RH
      !
      write(6,*) 'K-point parallelized integration time= ', t0
      ! 
      close(11)
      !     
      ! Clean memory
      !
      deallocate(eqkfit,etfull,tauk_full,vk,mk)
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
    
