
! prologue: gather the data and put it in temporary arrays.

   n=ipairs(m)
   itran=ishft(n,-27)
   n = iand(n,mask27)
   j = bckptr(n)
   delx = imagcrds(1,n) + xktran(1,itran)
   dely = imagcrds(2,n) + xktran(2,itran)
   delz = imagcrds(3,n) + xktran(3,itran)
   delr2 = delx*delx + dely*dely+delz*delz
write(72,*)i,j
!   if(maphb(i,j).gt.0) then
!  write(72,*)i,j,maphb(i,j),NXI(maphb(i,j)),NXI1(maphb(i,j)),NXI2(maphb(i,j)),NXI3(maphb(i,j))
!   end if

   do counter=1,nhb_pair
!   write(72,*)NXI(counter),NXI1(counter)
!   if(((i == NXI(counter)).and.(j == NXI1(counter))).or.((j == NXI(counter)).and.(i == NXI1(counter)))) then
   if(i == NXI(counter))  then
!	write(72,*) NXI(counter),NXI1(counter),i,j 

   delx_ixi_ixj(maphb(i,j)) = imagcrds(1,n) + xktran(1,itran)
   dely_ixi_ixj(maphb(i,j)) = imagcrds(2,n) + xktran(2,itran)
   delz_ixi_ixj(maphb(i,j)) = imagcrds(3,n) + xktran(3,itran)


   end if
   enddo
   do counter=1,nhb_pair

   if(((i == NXI2(counter)) .and. (j == NXI3(counter))) .or. ((j == NXI2(counter)) .and. (i == NXI3(counter)))) then
!        write(72,*) NXI2(maphb(i,j)),NXI3(maphb(i,j)),i,j

   delx_ixk_ixl(maphb(i,j)) = imagcrds(1,n) + xktran(1,itran)
   dely_ixk_ixl(maphb(i,j)) = imagcrds(2,n) + xktran(2,itran)
   delz_ixk_ixl(maphb(i,j)) = imagcrds(3,n) + xktran(3,itran)
!write(72,*) delx,delx_ixk_ixl(maphb(i,j))

   end if
   enddo


   if ( delr2 < filter_cut2 )then
      icount = icount + 1
      cache_x(icount) = delx
      cache_y(icount) = dely
      cache_z(icount) = delz
      cache_r2(icount) = delr2
      cache_bckptr(icount) = j
   end if
