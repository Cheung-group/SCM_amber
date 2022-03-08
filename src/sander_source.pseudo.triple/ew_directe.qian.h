
! epilogue: 12-6 LF terms
do im_new = 1,icount
   j = cache_bckptr(im_new)
   dfee = cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2inv = cache_r2(im_new)

!  write(72,*)"epilogue delx",delx
   ic = ico(iaci+iac(j))
   r6 = delr2inv*delr2inv*delr2inv
#ifdef LES 
   lfac=lesfac(lestmp+lestyp(j))
   f6 = cn2(ic)*r6*lfac
   f12 = cn1(ic)*(r6*r6)*lfac
   if(ipimd>0.or.ineb>0) then
      if(cnum(i).eq.0.and.cnum(j).eq.0) then
         nrg_all(1:nbead)=nrg_all(1:nbead) + (f12-f6)/nbead
      else
         if(cnum(i).ne.0) then
            nrg_all(cnum(i)) = nrg_all(cnum(i)) + f12-f6
         else
            nrg_all(cnum(j)) = nrg_all(cnum(j)) + f12-f6
         endif
      endif
   endif
#else
#  ifdef TVDW
   !         "truncated" LJ repulsion:
   r6pinv = 1.d0/(1.d0/r6 + 0.01d0)
   f12 = cn1(ic)*r6pinv*r6pinv
   f6 = cn2(ic)*r6pinv
#  else
   f6 = cn2(ic)*r6
   f12 = cn1(ic)*(r6*r6)

!   write(72,*)"im_new",im_new
!Antonios added here 3/15/2010 HB
   if(maphb(iac(i),iac(j)).lt.1) then
	EHBV=EHBV+f12-f6
   ! -- ti decomp
!   write(80,*)idecomp
	   if(decpr .and. idecomp > 0) call decpair(3,i,j,(f12 - f6)/(nstlim/ntpr))
	evdw = evdw + f12 - f6
        df = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
	dfx = delx*df
	dfy = dely*df
	dfz = delz*df
#ifndef noVIRIAL
	vxx = vxx - dfx*delx
	vxy = vxy - dfx*dely
	vxz = vxz - dfx*delz
	vyy = vyy - dfy*dely
	vyz = vyz - dfy*delz
	vzz = vzz - dfz*delz
#endif
	dumx = dumx + dfx
	dumy = dumy + dfy
	dumz = dumz + dfz
	force(1,j) = force(1,j) + dfx
	force(2,j) = force(2,j) + dfy
	force(3,j) = force(3,j) + dfz

!       write(72,*)im_new,j
                  else

! Yes HB for non-periodic	
        num= maphb(i,j)
        U=f12-f6

!   write(72,*)f12,f6
        df =  dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
        dfx = delx*df
        dfy = dely*df
        dfz = delz*df
#ifndef noVIRIAL
!        vxx = vxx - dfx*delx
!        vxy = vxy - dfx*dely
!        vxz = vxz - dfx*delz
!        vyy = vyy - dfy*dely
!        vyz = vyz - dfy*delz
!        vzz = vzz - dfz*delz
#endif

              ixi=NXI(num)
              ixj=NXI1(num)
              ixk=NXI2(num)
              ixl=NXI3(num)

!         write(72,*)im_new,i,j,ixi,ixj,ixk,ixl 
!!         check if ixk=i,ixj=j
      

                      XIJ = delx_ixi_ixj(maphb(i,j)) 
                      YIJ = dely_ixi_ixj(maphb(i,j)) 
                      ZIJ = delz_ixi_ixj(maphb(i,j))
                      XKJ = delx_ixj_ixk(maphb(i,j))
                      YKJ = dely_ixj_ixk(maphb(i,j))
                      ZKJ = delz_ixj_ixk(maphb(i,j))
                      XLK = delx_ixk_ixl(maphb(i,j)) 
                      YLK = dely_ixk_ixl(maphb(i,j)) 
                      ZLK = delz_ixk_ixl(maphb(i,j)) 


!          write(72,*)XKJ,YKJ,ZKJ,i,j
!          write(72,*)XIJ,YIJ,ZIJ,i,j
!          write(72,*)XLK,YLK,ZLK,i,j

        TX=(YIJ*ZKJ-YKJ*ZIJ)
        TY=(XKJ*ZIJ-XIJ*ZKJ)
        TZ=(XIJ*YKJ-XKJ*YIJ)

!	write(72,*)TX,TY,TZ

        UX=(YLK*ZKJ-YKJ*ZLK)
        UY=(XKJ*ZLK-XLK*ZKJ)
        UZ=(XLK*YKJ-XKJ*YLK)

!	write(72,*)UX,UY,UZ,i,j

        temp=TX*UX+TY*UY+TZ*UZ
        RT=TX*TX+TY*TY+TZ*TZ
        RT=sqrt(RT)
        RU=UX*UX+UY*UY+UZ*UZ
        RU=sqrt(RU)
        cosphi=temp/(RT*RU)

!	write(72,*)TX,TY,TZ,RT,temp,cosphi
!! ----- calculate structural factor
!!        energy
!! the unperturbed cosphi
!!           sinphi=(1-cosphi**2)
!!           sinphi=sqrt(sinphi)


       tempA=1-cosphi
       tempB=cosphi+1
       tempC=cosphi/cosphia
       tempC=1-tempC

       dummy=tempA*tempB*tempC
       dummy=1+Amp*dummy*dummy
       Chi=1/dummy
!! Antonios (evdw does not seem to have a - sign anywhere so we do evdw+chi*U)

               evdw = evdw+Chi*U

!	write(72,*)evdw
!	write(72,*)Chi,U,i,j

              EHBA = EHBA+Chi*U

!	write(72,*)EHBA
!!  Forces

!! dChi/dcosphi
  
       temp=-tempB*tempC+tempA*tempC-tempA*tempB/cosphia
       temp=temp*tempA*tempB*tempC
       temp=temp*2*Amp
       dChi=-temp/(dummy**2)

!	write(72,*)temp,i,j
!! partial derivatives

        RT2=RT*RT
        RU2=RU*RU
        RTU=RT*RU
        DTXX=UX/RTU-cosphi*TX/RT2
        DTY=UY/RTU-cosphi*TY/RT2
        DTZ=UZ/RTU-cosphi*TZ/RT2

        DUX=TX/RTU-cosphi*UX/RU2
        DUY=TY/RTU-cosphi*UY/RU2
        DUZ=TZ/RTU-cosphi*UZ/RU2

!! force dcos/dxi 
        FXI=DTY*(-ZKJ)+DTZ*(YKJ)
!write(72,*)FXI,i,j
!! force dcos/dyi
        FYI=DTXX*(ZKJ)+DTZ*(-XKJ)

!! force dcos/dzi
        FZI=DTXX*(-YKJ)+DTY*(XKJ)


!! force dcos/dxj
        FXJ=DTY*(ZKJ-ZIJ)+DTZ*(-YKJ+YIJ)+DUY*(-ZLK) + DUZ*(YLK)

!! force dcos/dyj
        FYJ=DTXX*(ZIJ-ZKJ)+DTZ*(XKJ-XIJ)+DUX*(ZLK) + DUZ*(-XLK)

!! force dcos/dzj
        FZJ=DTXX*(YKJ-YIJ)+DTY*(XIJ-XKJ)+DUX*(-YLK) + DUY*(XLK)



!! force dcos/dxk
        FXK=DTY*(ZIJ)+DTZ*(-YIJ)+DUY*(ZKJ+ZLK) + DUZ*(-YKJ-YLK)

!! force dcos/dyk
        FYK=DTXX*(-ZIJ)+DTZ*(XIJ)+DUX*(-ZKJ-ZLK) + DUZ*(XLK+XKJ)

!! force dcos/dzk
       FZK=DTXX*(YIJ)+DTY*(-XIJ)+DUX*(YLK+YKJ) + DUY*(-XLK-XKJ)


!! force dcos/dxl
        FXl= DUY*(-ZKJ) + DUZ*(YKJ)

!! force dcos/dyl
       FYl= DUX*(ZKJ) + DUZ*(-XKJ)

!! force dcos/dzl
        FZl= DUX*(-YKJ) + DUY*(XKJ)

!! force -dP/dx

        FFxi=-U*dChi*FXI
        FFyi=-U*dChi*FYI
        FFzi=-U*dChi*FZI
        FFxj=-U*dChi*FXJ
        FFyj=-U*dChi*FYJ
        FFzj=-U*dChi*FZJ
        FFxk=-U*dChi*FXK
        FFyk=-U*dChi*FYK
        FFzk=-U*dChi*FZK
        FFxl=-U*dChi*FXL
        FFyl=-U*dChi*FYL
        FFzl=-U*dChi*FZL

!! -dp/dxi
                      N_HB=NXI(num)

                       force(1,N_HB) = force(1,N_HB) + FFxi
                       force(2,N_HB) = force(2,N_HB) + FFyi
                       force(3,N_HB) = force(3,N_HB) + FFzi

!! -dp/dxj
!                      N_HB=NXI1(num)
                      N_HB=j
                       force(1,N_HB) = force(1,N_HB) + FFxj + Chi*dfx
                       force(2,N_HB) = force(2,N_HB) + FFyj + Chi*dfy
                       force(3,N_HB) = force(3,N_HB) + FFzj + Chi*dfz

!              write(72,*)N_HB,i,j

!! -dp/dxk
!                      N_HB=NXI2(num)
                      N_HB=i
                       force(1,N_HB) = force(1,N_HB) + FFxk - Chi*dfx
                       force(2,N_HB) = force(2,N_HB) + FFyk - Chi*dfy
                       force(3,N_HB) = force(3,N_HB) + FFzk - Chi*dfz

!! dp/dxl
                      N_HB=NXI3(num)
                       force(1,N_HB) = force(1,N_HB) + FFxl
                       force(2,N_HB) = force(2,N_HB) + FFyl
                       force(3,N_HB) = force(3,N_HB) + FFzl

	
  	   endif
     
       
!Antonios end 

	if((im_new .eq. 1) .and. (i .eq. 1)) then
!	write(73,*)"haha"
	do i_chi=1,ICHI

        ixi=ICA(i_chi)
        ixi1=ICB(i_chi)
        ixi2=INC(i_chi)
        ixi3=ICC(i_chi)

        XAX = delx_xax(i_chi)
        XAY = dely_xay(i_chi)
        XAZ = delz_xaz(i_chi)

        XBX = delx_xbx(i_chi)
        XBY = dely_xby(i_chi)
        XBZ = delz_xbz(i_chi)

        XCX = delx_xcx(i_chi)
        XCY = dely_xcy(i_chi)
        XCZ = delz_xcz(i_chi)


        xtriple=  XAX*XBY*XCZ+ XAZ*XBX*XCY + XAY*XBZ*XCX- XAZ*XBY*XCX - XAY*XBX*XCZ- XAX*XBZ*XCY

! calculate energy

        temp_chi=xkchi*(xtriple-XCHI(i_chi))

!        write(85,*) 'temp_chi ', xkchi, xtriple, XCHI(i_chi) 


        ECHII=0.5*temp_chi*(xtriple-XCHI(i_chi))
        ECHI=ECHII+ECHI
!       evdw = evdw + ECHI
!       if(i_chi.eq.1) then
!       write(81,*)xchi(i_chi),xtriple,ECHII,temp_chi
!       write(72,*)ixi,ixi1,ixi2,ixi3
!       endif
! calculate force

        Fxi1_chi= -(XBY*XCZ-XBZ*XCY)*temp_chi
        Fxi2_chi= -(XAZ*XCY-XAY*XCZ)*temp_chi
        Fxi3_chi= -(XAY*XBZ-XAZ*XBY)*temp_chi
	Fxi_chi = -(Fxi1_chi+Fxi2_chi+Fxi3_chi)

        Fyi1_chi= -(XBZ*XCX-XBX*XCZ)*temp_chi
        Fyi2_chi= -(XAX*XCZ-XAZ*XCX)*temp_chi
        Fyi3_chi= -(XAZ*XBX-XAX*XBZ)*temp_chi
        Fyi_chi=   -(Fyi1_chi+Fyi2_chi+Fyi3_chi)

!       if(i_chi.eq.1) then
!       write(83,*)fyi1,fyi2,fyi3
!       endif

        Fzi1_chi= -(XBX*XCY-XBY*XCX)*temp_chi
        Fzi2_chi= -(XAY*XCX-XAX*XCY)*temp_chi
        Fzi3_chi= -(XAX*XBY-XAY*XBX)*temp_chi
        Fzi_chi = -(Fzi1_chi+Fzi2_chi+Fzi3_chi)

!!      if(i_chi.eq.1) then
!       write(79,*)fzi1,fzi2,fzi3
!       endif


        if(Fxi_chi.gt.MAXLIM) then
        write(84,*) 'ixi,x ', ixi,fxi_chi,temp_chi,xtriple
        Fxi_chi=MAXLIM
        endif
        if(Fyi_chi.gt.MAXLIM) then
        write(84,*) 'ixi,y ', ixi,fyi_chi,temp_chi,xtriple
        Fyi_chi=MAXLIM
        endif
        if(Fzi_chi.gt.MAXLIM) then
        Fzi_chi=MAXLIM
        write(84,*) 'ixi,z ', ixi,fzi_chi,temp_chi,xtriple
        endif

        if(Fxi1_chi.gt.MAXLIM) then
        write(84,*) 'ixi1,x ', ixi1,fxi1_chi,temp_chi,xtriple
        Fxi1_chi=MAXLIM
        endif
        if(Fyi1_chi.gt.MAXLIM) then
        write(84,*) 'ixi1,y ', ixi1,fyi1_chi,temp_chi,xtriple
        Fyi1_chi=MAXLIM
        endif
        if(Fzi1_chi.gt.MAXLIM) then
        write(84,*) 'ixi1,z ', ixi1,fzi1_chi,temp_chi,xtriple
        Fzi1_chi=MAXLIM
        endif
	if(Fxi2_chi.gt.MAXLIM) then
        write(84,*) 'ixi2,x ', ixi2,fxi2_chi,temp_chi,xtriple
        Fxi2_chi=MAXLIM
        endif
        if(Fyi2_chi.gt.MAXLIM) then
        write(84,*) 'ixi2,y ', ixi2,fyi2_chi,temp_chi,xtriple
        Fyi2_chi=MAXLIM
        endif
        if(Fzi2_chi.gt.MAXLIM) then
        write(84,*) 'ixi2,z ', ixi2,fzi2_chi,temp_chi,xtriple
        Fzi2_chi=MAXLIM
        endif

        if(Fxi3_chi.gt.MAXLIM) then
        write(84,*) 'ixi3,x ', ixi3,fxi3_chi,temp_chi,xtriple
        Fxi3_chi=MAXLIM
        endif
        if(Fyi3_chi.gt.MAXLIM) then
        write(84,*) 'ixi3,y ', ixi3,fyi3_chi,temp_chi,xtriple
        Fyi3_chi=MAXLIM
        endif
        if(Fzi3_chi.gt.MAXLIM) then
        write(84,*) 'ixi3,z ', ixi3,fzi3_chi,temp_chi,xtriple
        Fzi3_chi=MAXLIM
        endif

        force(1,ixi) = force(1,ixi)+Fxi_chi
        force(2,ixi) = force(2,ixi)+Fyi_chi
        force(3,ixi) = force(3,ixi)+Fzi_chi
        force(1,ixi1) = force(1,ixi1)+Fxi1_chi
        force(2,ixi1) = force(2,ixi1)+Fyi1_chi
        force(3,ixi1) = force(3,ixi1)+Fzi1_chi
        force(1,ixi2) = force(1,ixi2)+Fxi2_chi
        force(2,ixi2) = force(2,ixi2)+Fyi2_chi
        force(3,ixi2) = force(3,ixi2)+Fzi2_chi
        force(1,ixi3) = force(1,ixi3)+Fxi3_chi
        force(2,ixi3) = force(2,ixi3)+Fyi3_chi
        force(3,ixi3) = force(3,ixi3)+Fzi3_chi

        enddo
    evdw = evdw + echi
	endif
#  endif
#endif
end do  !  im_new = 1,icount
