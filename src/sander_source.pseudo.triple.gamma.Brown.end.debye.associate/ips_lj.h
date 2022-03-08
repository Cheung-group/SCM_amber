
! epilogue: 12-6 LF IPS terms  by Xiongwu Wu

do im_new = 1,icount
   j = cache_bckptr(im_new)

   dfee = cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2inv = cache_r2(im_new)

   ic = ico(iaci+iac(j))
   r6 = delr2inv*delr2inv*delr2inv
   f6 = cn2(ic)*r6
   f12 = cn1(ic)*(r6*r6)
!  L-J r6 term
        UIPS2=1.0D0/(DELR2INV*RIPS2)
        UIPS6=UIPS2*UIPS2*UIPS2
        TWOU2=2.0D0-UIPS2
        TWOU6=TWOU2*TWOU2*TWOU2
        TWOU12=TWOU6*TWOU6
        PVC=BIPSVC0+UIPS2*(BIPSVC1    &
             +UIPS2*(BIPSVC2+UIPS2*BIPSVC3))
            PVCU=2.0D0*BIPSVC1+       &
                  UIPS2*(4.0D0*BIPSVC2+6.0D0*BIPSVC3*UIPS2)
            DVCU=(PVCU+6.0D0*PVC/TWOU2)/TWOU6
!  L-J r12 term 
            PVA=BIPSVA0+UIPS2*(BIPSVA1    &
             +UIPS2*(BIPSVA2+UIPS2*BIPSVA3))
            PVAU=2.0D0*BIPSVA1+           &
             UIPS2*(4.0D0*BIPSVA2+6.0D0*BIPSVA3*UIPS2)
            DVAU=(PVAU+12.0D0*PVA/TWOU2)/TWOU12
   evdw = evdw + f12 - f6  &
         +(f12*(PVA/TWOU12-PIPSVAC)*UIPS6-f6*(PVC/TWOU6-PIPSVCC))*UIPS6
   df = dfee + (12.d0*f12 - 6.d0*f6  &
         -(f12*DVAU*uips6-f6*DVCU)*uips6*UIPS2)*delr2inv
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
end do  !  im_new = 1,icount
