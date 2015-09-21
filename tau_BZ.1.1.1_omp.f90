      program tau_BZ
      !---------------------------------------------------------------------------------
      ! Written by Burak H. v1.1.1_omp
      ! Computes el-ph scattering rates using the Frohlich model in the full BZ.
      ! We use only LO mode scattering, therefore g_{qv} is from
      ! Frohlich model. 
      ! 
      ! This version is to be used with transp.x.x.x.f90
      !
      ! Attempting to compute
      ! \tau_{nk}^{-1} = 2\pi/hbar \sum_{qv,m} |g_{q,v}|^2 {
      ! (n(w_v)+f_{m,k+q}) \delta(e_{m,k+q} - e_{nk} - \hbar w_v) +
      ! (1+n(w_v)-f_{m,k+q}) \delta(e_{m,k+q} - e_{nk} + \hbar w_v)}
      !
      !
      ! Files required: 
      !  1. prefix.a2Fsave : for band energies and velocities
      !  2. fil_info       : at(3,3) and bg(3,3)
      !  3. Cnu.txt        : el-ph parameters
      !  4. wo.in          : LO frequencies
      !  
      ! Parallelization:
      ! This version uses OpenMP parallelization. Set nthreads=# of cores you want to
      ! use. You have to use a single node! 
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
      ! Approximation to FD distribution function
      interface
        function w1gauss (x,ngauss)
           double precision :: w1gauss
           double precision, intent(in) :: x
           integer, intent(in) :: ngauss
        end function w1gauss
      end interface
      !
      !
      integer :: i,j,k,nu,itemp,nk1fit,nk2fit,nk3fit,nkfit,             &
     &           nksfit_real, nbnd, nksfit, npk, nsym,                  &
     &           s(3,3,48),ns,nrot,ibnd,io,phband_i, phband_f,          &
     &           nphband, n, nn, jbnd, nqtot, counter
      !
      integer :: ik, iq, ii, iik, iiq, jjq, ikq
      !
      double precision :: wk, at(3,3), bg(3,3), efermi,                 &
     &                    etk, etq, temp1, temp2, q2, sumk,             &
     &                    temperature, fd, be, wo(3), w0g1, w0g2,       &
     &                    degauss, al(3), k_plus_q(3), lenkpq, lenxq
      ! 
      logical :: soc
      !
      double precision, allocatable :: xkfit(:,:), etfit(:,:), wkfit(:),&
     &                                 tauk(:,:), xqfit(:,:), wqfit(:)
      !
      integer, allocatable :: eqkfit(:), sfit(:)
      !
      ! OMP variables
      double precision :: t0
      !
      integer :: nthreads
      !
      character*20 :: fil_a2F, fil_info
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

      namelist /input/ fil_info, fil_a2F,                               &
     &                 nk1fit, nk2fit, nk3fit, phband_i,                &
     &                 temperature, phband_f, efermi, degauss, nthreads,&
     &                 soc

      read(5,input)

      !$ call omp_set_num_threads(nthreads)

      npk = 40000
      !
      nkfit=nk1fit*nk2fit*nk3fit ! size of regular grid 
      !
      if (soc .eqv. .true.) then
         wk = 1.d0 / nkfit  ! weight for k-points in regular grid 
      else
         wk = 2.0 / nkfit
      end if
      !
      ! Convert temperature from K to Ryd
      temperature = temperature * KtoRy 
      !
      ! Convert efermi and degauss from eV to Ryd
      !
      efermi = efermi / RytoeV
      degauss = degauss / RytoeV
      !
      ! Read the POP frequencies in cm-1 (just Gamma should be)
      !
      open(11,file='wo.in',status='unknown')
          read(11,*) (wo(i),i=1,3)
      close(11)
      !
      ! Convert to Ryd
      wo = wo / Rytocm1 
      !
      ! LO-phonon couplings (Ryd*Ryd) 
      open(11,file='Cnu.txt',status='unknown')
      do nu=1,3
         read(11,*) al(nu)
      end do
      close(11)
      !
      !
      nphband = phband_f - phband_i + 1 ! Total number of bands of interest (usually the number of relevant conduction bands)
      !
      ! Read a2F
      !
      open(11,file=fil_a2F,status='unknown')
      !
      read(11,*) nbnd, nksfit
      !
      allocate(etfit(nbnd,nksfit), xkfit(3,nksfit), xqfit(3,nksfit))
      allocate(wkfit(nksfit), wqfit(nksfit))
      !         etfit : band energies in IBZ
      !         xkfit : k-point coordinates
      !         xqfit : q-point coordinates (equivalent to k)
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
      !
      xqfit = xkfit        ! Same grid
      if (soc .eqv. .true.) then
         wqfit = wkfit 
      else
         wqfit = wkfit / 2.0
      end if
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
      ! 
      allocate (eqkfit(nkfit), sfit(nkfit))
      !         eqkfit : band energies in uniform grid -- k
      !         sfit   : symmetries
      ! 
      call lint ( nsym, s, .true., at, bg, npk, 0,0,0,                  &
     &  nk1fit,nk2fit,nk3fit, nksfit, xkfit, 1, nkfit, eqkfit, sfit)
      ! 
      ! eqkfit(nk) : maps IBZ to full grid. The full grid is in crystal coords
      !
      deallocate (sfit)
      !
      ! Allocate tauk in IBZ 
      !
      allocate(tauk(nksfit,nphband))
      !
      ! Shift to all positive values
      !
      !$omp parallel do default(none) &
      !$omp collapse(2) &
      !$omp private(ik,i) shared(nksfit,xqfit,xkfit)
      do ik=1,nksfit
      do i=1,3
          !
          if (xkfit(i,ik) .le. -0.5 ) then
             xkfit(i,ik) = xkfit(i,ik) + 1.0
          end if
          !
          if (xqfit(i,ik) .le. -0.5 ) then
             xqfit(i,ik) = xqfit(i,ik) + 1.0
          end if
          !
      end do
      end do
      !$omp end parallel do
      !
      ! Compute the sum over q-points to obtain tauk
      !
      tauk = 0.d0
      !
      ! Start timimg here
      !$ t0 = omp_get_wtime()
      !
      do ik=1,nksfit
         !
         do ibnd=phband_i, phband_f 
            ii = ibnd - phband_i + 1
            !
            !$omp parallel do default(shared) &
            !$omp private(i,iq,ikq,jbnd,nu,q2,etq,be,temp1,temp2,w0g1,w0g2,k_plus_q) &
            !$omp reduction(+ : tauk) 
            do iq = 2, nksfit  ! Avoid the q=0 divergence
               !
               q2 = 0.d0
               do i=1,3
                  q2 = q2 + xqfit(i,iq)**2.0
               end do
               ! 
               ! Find the eigenvalue at k+q
               !
               k_plus_q(:) = xqfit(:,iq) + xkfit(:,ik)
               !
               call find_kpq (xqfit,k_plus_q,nksfit,ik,iq,ikq)
               !
               ! Sum over bands jbnd and modes nu
               !
               do jbnd=phband_i,phband_f
                  !
                  ! Band energies e_k and e_{q}
                  !
                  etk=etfit(ibnd,ik)   ! e_{ibnd,k}
                  etq=etfit(jbnd,ikq)   ! e_{jbnd,q}
                  ! 
                  do nu=1,3
                     !
                     ! BE distribution
                     be = 1.d0 / ( exp(wo(nu)/temperature) - 1.d0 )
                     !
                     ! n_{q\nu} + f_{m,q}
                     temp1 = w1gauss((efermi-etq)/temperature,1)+       &
     &                       be
                     ! 1 - f_{m,q} + n_{q\nu}
                     temp2 = 1.d0 - w1gauss((efermi-etq)/temperature,1)+&
     &                       be
                     !
                     ! Energy conserving delta-functions 
                     w0g1 = w0gauss((etq-etk-wo(nu))/degauss)/degauss
                     w0g2 = w0gauss((etq-etk+wo(nu))/degauss)/degauss
                     ! tau
                     !if ( abs(q2) .ge. 1.d-20 ) then ! Avoid the divergence at q-k=0
                     tauk(ik,ii) = tauk(ik,ii) + tpi * wqfit(iq) * al(nu) / q2 *  &
     &                      ( temp1 * w0g1 + temp2 * w0g2 )
                     !end if
                     !
                     !
                  end do ! nu
               end do ! jbnd
            end do ! iq
            !$omp end parallel do
            !
         end do ! ibnd
      end do ! ik   
      !
      ! Measure integration time
      !$ t0 = omp_get_wtime() - t0
      !
      open(11,file='tauBZ-tr.dat',status='unknown')
      !
      do ik=1,nksfit
         write(11,'(3f14.6)') (tauk(ik,ibnd),ibnd=1,nphband)
      end do
      !
      close(11)
      ! 
      write(6,*) 'Integration walltime= ', t0
      !     
      ! Clean memory
      !
      deallocate(eqkfit,wqfit,xqfit,tauk,wkfit,etfit,xkfit)
      !
      end program tau_BZ

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
      !
      !
      !    
      function w1gauss (x,ngauss)

      implicit none

      double precision :: w1gauss

      double precision, intent(in) :: x

      integer, intent(in) :: ngauss

      double precision :: sqrtpm1

      sqrtpm1 = 1.0d0/1.77245385090551602729d0

      ! ngauss=0 (gaussian smearing)
      ! ngauss=1 FD

      if ( ngauss .eq. 0 ) then
         w1gauss = 0.5d0 * ( 1.d0 + derf( x * sqrtpm1 ) )
      else if ( ngauss .eq. 1 ) then
         w1gauss = 1.d0 / ( exp(-x) + 1.d0)
      end if

      end function w1gauss

