      subroutine find_kpq(xqfit,k_plus_q,nksfit,ik,iq,ikq) 

      implicit none

      integer, intent(in) :: nksfit,ik,iq

      double precision, intent(in) :: xqfit(3,nksfit), k_plus_q(3)

      integer, intent(out) :: ikq

      integer :: i,j, iiq

      double precision :: temp, kpq(3)

      kpq = k_plus_q

      ! Fold back to 1st BZ
      !
      do i=1,3
         if (kpq(i) .gt. 0.5) then
             kpq(i) = 1.0 - kpq(i)
         else if (kpq(i) .ge. 1.0) then
             kpq(i) = kpq(i) - 1.0
         end if
      end do 
      !
      ! Swap to right order
      !
      do i=2,3
         j=i-1
         temp = kpq(i)
         do while ( (j >= 1) .and. (kpq(j) > temp)) 
            kpq(j+1) = kpq(j)
            j = j-1
         end do
         kpq(j+1) = temp
      end do
      !    
      ikq = 0
      !
      ! For cubic only!
      !
      do iiq=1,nksfit 
         !
         if (abs(kpq(1)-xqfit(1,iiq)) .lt. 1.d-8 .and.                  &
     &       abs(kpq(2)-xqfit(2,iiq)) .lt. 1.d-8 .and.                  &
     &       abs(kpq(3)-xqfit(3,iiq)) .lt. 1.d-8 ) then
             ikq = iiq
!             write(*,'(3i3,3f14.6)') ik, iq, ikq, (kpq(i),i=1,3)
          end if
          ! 
      end do


      end subroutine find_kpq
