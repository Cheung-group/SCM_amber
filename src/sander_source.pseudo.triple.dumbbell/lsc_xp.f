!#include "assert.h"
#include "dprec.h"
!
!------------------------------------------------------------------------------
subroutine lsc_xp(x,p)
!
  use lscivr_vars
  use constants, only : pi, hbar, kB
  implicit none
# include "extra.h"
# include "parallel.h"
!
  _REAL_, intent(inout):: x(ndof_lsc), p(ndof_lsc)
  integer:: i, ierr, j
  _REAL_ :: pxcom, pycom, pzcom, psqcom, psqnew, rat_psq, psum    
! an auxiliariary matrix for the Kubo-transformed corr. func.
  _REAL_, allocatable :: mat_kubo(:,:), p_kubo(:)
  integer:: info, nrand
  _REAL_:: xrand
!
   allocate(mat_kubo(ndof_lsc,ndof_lsc), stat = ierr)
   allocate(p_kubo(ndof_lsc), stat = ierr)
!   REQUIRE(ierr == 0)
!

! determine the randome number, this is necessary because each run has to be independent

    nrand =  0 
    !read(file_rhoa_lsc,*)nrand
    !if(nrand == 0)then
    !   call gauss(0.0d0, 1.0d0, xrand)  
    !   rewind(file_rhoa_lsc)
    !else      
    !   do i = 1, nrand
    !      call gauss(0.0d0, 1.0d0, xrand)
    !   enddo  
    !endif

! 2nd derivative of the mass-weighted potenital
   do ilsc = 1, ndof_lsc
      v2_lsc(ilsc,1:ndof_lsc) = v2_lsc(ilsc,1:ndof_lsc)/sqrt(mass_lsc(ilsc))
   end do
   do ilsc = 1, ndof_lsc
      v2_lsc(1:ndof_lsc,ilsc) = v2_lsc(1:ndof_lsc,ilsc)/sqrt(mass_lsc(ilsc))
   end do   


!   write(100*worldsize+worldrank, * ) 'dumping v2_lsc before dsyev'

!   do ilsc = 1, ndof_lsc
!       write(100*worldsize+worldrank, *) 'vs_lsc(:,', ilsc, ')'
!       write(100*worldsize+worldrank, '(8F10.4)' ) v2_lsc(:,ilsc)
!   end do


!   write(file_rhoa_lsc,*)natom_lsc, ndof_lsc
!   do ilsc = 1, ndof_lsc
!      do i = 1, ndof_lsc
!	write(file_rhoa_lsc,*)ilsc, i, v2_lsc(ilsc,i)
!      end do
!   end do

!
! diagonalize the mass-weighted Hesian matrix
   call dsyev( 'V', 'L', ndof_lsc, v2_lsc, ndof_lsc, v2_freq, & 
               work_lsc, 10*ndof_lsc, info )

!   write(100*worldsize+worldrank, * ) 'dumping v2_lsc(1:10,1) after dsyev'
!   write(100*worldsize+worldrank, '(5E16.4)' ) (v2_lsc(i,1), i=1,10)

!   write(100*worldsize+worldrank, * ) 'dumping v2_lsc(1:10,2) after dsyev'
!   write(100*worldsize+worldrank, '(5E16.4)' ) (v2_lsc(i,2), i=1,10)

!   write(100*worldsize+worldrank, * ) 'dumping v2_lsc(1:10,3) after dsyev'
!   write(100*worldsize+worldrank, '(5E16.4)' ) (v2_lsc(i,3), i=1,10)



! Now the array 'v2_ivr' stores the eigenvectors
   p_lsc = 0.0d0
   sigma_lsc = 0.0d0
   do ilsc = 1, ndof_lsc
   if ( v2_freq(ilsc) >= 1.0d-48 ) then
      omega_lsc = sqrt( v2_freq(ilsc) )
      alpha_lsc(ilsc) = omega_lsc/hbar*cosh(0.50d0*beta_lsc*hbar*omega_lsc) &
                    /sinh(0.50d0*beta_lsc*hbar*omega_lsc)
   elseif( abs(v2_freq(ilsc)) .lt. 1.0d-48 )then 
      omega_lsc = sqrt( abs(v2_freq(ilsc)) )
      alpha_lsc(ilsc) = cosh(0.50d0*beta_lsc*hbar*omega_lsc)/(0.5d0*beta_lsc*hbar**2)     
   else
      omega_lsc = sqrt( - v2_freq(ilsc) )
      if ( ( beta_lsc * hbar * omega_lsc ) < pi ) then
         alpha_lsc(ilsc) = omega_lsc/hbar*cos(0.5d0*beta_lsc*hbar*omega_lsc) &
                       /sin(0.5d0*beta_lsc*hbar*omega_lsc)
      else
         alpha_lsc(ilsc) = 0.0d0
      endif
   endif
   enddo
!      write(file_rhoa_lsc,*)v2_lsc
!
!!   do ilsc = 1, ndof_lsc
!!      write(file_rhoa_lsc,'(I6,2x,2x,e12.6,2x,e12.6)') ilsc, v2_freq(ilsc), alpha_lsc(ilsc)
!!   end do
!   write(file_rhoa_lsc,*)'-----------------------------------'
!   write(file_rhoa_lsc,*)'v2_lsc='
!   write(file_rhoa_lsc,*)v2_lsc
!   write(file_rhoa_lsc,*)'-----------------------------------'
!   stop
!

!   alpha_lsc = 1.0d0


   do ilsc = 1, ndof_lsc
      if(alpha_lsc(ilsc) .gt. 0.0d0)then
         sigma_lsc = sqrt( 0.5d0 * hbar**2 * alpha_lsc(ilsc) )
         call gauss( 0.0d0, sigma_lsc, p_lsc(ilsc) )
         nrand = nrand + 1
      else
         p_lsc(ilsc) = 0.0d0
      end if
   end do
!


!   write(100*worldsize+worldrank, * ) 'dumping p_lsc:'
!   write(100*worldsize+worldrank, '(5E16.8)') (p_lsc(i), i=1,ndof_lsc)

!   write(100*worldsize+worldrank, * ) 'dumping v2_lsc(1,:)'
!   write(100*worldsize+worldrank, '(5E16.4)' ) (v2_lsc(1,i), i=1,ndof_lsc)

   do i = 1, ndof_lsc
      if(i==1) then
         p(i) = 0.0
         do j=1, ndof_lsc
            p(i) = p(i) + v2_lsc(1,j)*p_lsc(j)
            write(100*worldsize+worldrank, '(4E16.4)') v2_lsc(1,j), p_lsc(j), v2_lsc(1,j)*p_lsc(j), p(i)
         end do
      end if
      p(i) = sum( v2_lsc(i,:) * p_lsc(:) )* sqrt(mass_lsc(i))
   enddo


!   write(100*worldsize+worldrank, * ) 'dumping velocity:'
!   write(100*worldsize+worldrank, '(5E16.8)') (p(i), i=1,ndof_lsc)

   psqcom = 0.0d0
   do i = 1, ndof_lsc
      psqcom = psqcom + 0.50d0*p(i)**2/mass_lsc(i)
   end do
!!   write(file_rhoa_lsc,*)'temp=',psqcom/kB/natom_lsc*2.0d0/3.0d0
!   stop

!
   ! remove the momentum of the center
   pxcom = 0.0d0
   pycom = 0.0d0
   pzcom = 0.0d0
   psqcom = 0.0d0
!
   do i = 1, ndof_lsc
      psqcom = psqcom + p(i)**2/mass_lsc(i)
   end do
!
   do i = 1, natom_lsc
      pxcom = pxcom + p(3*i-2)
      pycom = pycom + p(3*i-1)
      pzcom = pzcom + p(3*i)
   end do
!
     pxcom = pxcom / dble(natom_lsc)
     pycom = pycom / dble(natom_lsc)
     pzcom = pzcom / dble(natom_lsc)
!
!!   write(file_rhoa_lsc,*)pxcom, pycom, pzcom
!
   do i = 1, natom_lsc
      p(3*i-2) = p(3*i-2) - pxcom
      p(3*i-1) = p(3*i-1) - pycom
      p(3*i)   = p(3*i) -   pzcom
   end do
!
   psqnew = 0.0d0
   do i = 1, ndof_lsc
       psqnew = psqnew + p(i)**2/mass_lsc(i)
   end do
!
   rat_psq = psqcom/psqnew
!!   write(file_rhoa_lsc,*)'ratio=', rat_psq
   p = p * sqrt(rat_psq)

   psqnew = 0.0d0
   do i = 1, ndof_lsc
      psqnew = psqnew + 0.50d0*p(i)**2/mass_lsc(i)
   end do
!!   write(file_rhoa_lsc,*)'temp=',psqnew/kB/natom_lsc*2.0d0/3.0d0


!  At this point, transfer the momentum p(:) 
!   to the velocity v(:)
   do i = 1, ndof_lsc 
       p(i) = p(i) / mass_lsc(i)
!       write(file_rhoa_lsc,*)i,p(i)
   end do 
!
   p_kubo = 0.0d0
   do i = 1, ndof_lsc
      p_kubo(i) = sqrt(mass_lsc(i)) * p(i)
   end do
   p_kubo = matmul(transpose(v2_lsc),p_kubo) 

!   do i = 1, ndof_lsc
!      p_kubo(i) = p_kubo 
!   end do
!
   rewind(file_rhoa_lsc)
   write(file_rhoa_lsc,*)nrand
   write(file_rhoa_lsc,*)p
   write(file_rhoa_lsc,*)x
   write(file_rhoa_lsc,*)mass_lsc
!  stop

! the sign to determine which type of the correlation function
! icorf_lsc = 0,       scalar function
!             1,       position <x(0)x(t)>
!             2,       velocity <v(0)v(t)>
!             3,       nonlinear function <f(0)f(t)>
!             4,       Kubo-transformed <v(0)v(t)>
! Also see lsc_init.f lscivr_vars.f
   if(icorf_lsc == 1)then 
       rho_a = x
   elseif(icorf_lsc == 2)then
       rho_a = p
   elseif(icorf_lsc == 4)then
       mat_kubo = 0.0d0
       mat_kubo = transpose(v2_lsc)
       do i = 1, ndof_lsc
          mat_kubo(:,i) = mat_kubo(:,i) / sqrt(mass_lsc(i))            
       end do
!       write(file_rhoa_lsc,*) 1.0d0/(beta_lsc*kb)
!       write(file_rhoa_lsc,*)mass_lsc  
!       write(file_rhoa_lsc,*)mat_kubo

       do i = 1, ndof_lsc
          if(alpha_lsc(i) < 1.0d-100)then
             p_kubo(i) = 0.0d0 !0.50 * 1.0d0 * hbar**2 !*mass_lsc(i)
          else  
!             p_kubo(i) = p_lsc(i)**2 / alpha_lsc(i) !* mass_lsc(i)
!             p_kubo(i) = p_lsc(i)/alpha_lsc(i) * p_kubo(i)
              p_kubo(i) = p_kubo(i)/alpha_lsc(i)
          endif
       end do
! 
       p_kubo = matmul(v2_lsc, p_kubo) 
       p_kubo = sqrt(mass_lsc)*p_kubo
          
!       write(file_rhoa_lsc,*)'-----------------------------------'
!       write(file_rhoa_lsc,*) 'mat_kubo='
!       write(file_rhoa_lsc,*) mat_kubo
!       write(file_rhoa_lsc,*)'-----------------------------------'

!       mat_kubo = matmul(mat_kubo, transpose(v2_lsc))

!       do i = 1, ndof_lsc
!          mat_kubo(:,i) = mat_kubo(:,i) * sqrt(mass_lsc(i))
!          write(file_rhoa_lsc,*) i,'mat_kubo='
!          write(file_rhoa_lsc,*) mat_kubo(i,:)
!       end do

!       mat_kubo = transpose(mat_kubo)

!       write(file_rhoa_lsc,*)'-----------------------------------'
!       write(file_rhoa_lsc,*)'mat_kubo='
!       write(file_rhoa_lsc,*)mat_kubo
!       write(file_rhoa_lsc,*)'-----------------------------------' 

!       rho_a = 0.0d0
!       do i = 1, ndof_lsc
!         rho_a(i) = sum(mat_kubo(i,:)*p(:)) 
!       end do
       do i = 1, ndof_lsc
       rho_a(i) = p_kubo(i)  !/mass_lsc(i)
!       write(*,*)mass_lsc(i)
       end do

       rho_a = rho_a * 2.0d0/(hbar**2 * beta_lsc) /(2*kb)
   endif
!
   deallocate(mat_kubo, stat = ierr)
   deallocate(p_kubo, stat = ierr)

!  REQUIRE(ierr == 0)  
   if(icorf_lsc == 0)then
     write(file_rhoa_lsc, '(f21.9)')    &
          rho_a
   else
     psum = 0.0d0
     do i = 1, natom_lsc
!        psum = psum + rho_a(3*i-2) + rho_a(3*i-1) + rho_a(3*i)
        psum = psum + rho_a(3*i-2)*p(3*i-2) + rho_a(3*i-1)*p(3*i-1) + rho_a(3*i)*p(3*i)
!       write(file_rhoa_lsc,*)'vector ='
       write(file_rhoa_lsc, '(3f21.9)')  &
              rho_a(3*i-2:3*i)
     end do
     psum = psum / dble(natom_lsc)
     write(file_rhoa_lsc,*)'psum=',psum
   endif

!   stop

   return
end subroutine lsc_xp
