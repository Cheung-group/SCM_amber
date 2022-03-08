
! prologue: gather the data and put it in temporary arrays.

   n=ipairs(m)
!   write(72,*)"N=",n,"M=",m
   itran=ishft(n,-27)
!   write(72,*)itran,"N",n
!   call EXIT
   n = iand(n,mask27)
   j = bckptr(n)
!   write(72,*)"N",n,"I",i,"J",j
!   write(72,*)j
!   call EXIT 
   delx = imagcrds(1,n) + xktran(1,itran)
!    XKJ = imagcrds(1,n) - xk + tranvec(1,itran)
!write(72,*)n,imagcrds(1,n),xktran(1,itran),delx,i,j,m
!write(72,*)"prologue delx",delx
   dely = imagcrds(2,n) + xktran(2,itran)
   delz = imagcrds(3,n) + xktran(3,itran)
   delr2 = delx*delx + dely*dely+delz*delz

!write(74,*)delx,XKJ,n
   
!write(72,*)"icount",icount
   if ( delr2 < filter_cut2 )then
      icount = icount + 1
      cache_x(icount) = delx
      cache_y(icount) = dely
      cache_z(icount) = delz
      cache_r2(icount) = delr2
      cache_bckptr(icount) = j
   end if
