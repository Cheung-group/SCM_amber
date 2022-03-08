! epilogue: softcore (modified 12-6) LJ terms for V0, equals vanishing atoms

do im_new = 1,icount
   j = cache_bckptr(im_new) ! atom# of atom j

   df =   cache_df(im_new)     ! electrostatic energy times 1/r^2 for the forces
   delx = cache_x(im_new)      ! delta x,y and z
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2 = cache_r2(im_new)    ! this contains r^2, note the difference from ew_directe2

   rfour=delr2*delr2
   r6=rfour*delr2

   ic = ico(iaci+iac(j))      ! locate the index in the vdW parameter array

   if ( nsc(i) == nsc(j) ) then
      !both atoms are softcore atoms and have normal 6-12 vdw
      denom = 1.0d0 / ( r6 * sigma6(ic) )
      denom2 = denom * denom

      sc_ener(7) = sc_ener(7) + foureps(ic) * ( denom2 - denom ) ! Potential goes into the softcore energy array
      df = df + oneweight * foureps(ic) * ( 12.0d0 * denom2 - 6.0d0 * denom ) / delr2 ! scaled up by oneweight
   else
      ! use the softcore potential
      denom = 1.0d0 / ( scalpha * clambda + r6 * sigma6(ic) ) !sigma6 is 1/(sigma^6)
      denom2 = denom * denom
      denom3 = denom2 * denom

      evdw = evdw + foureps(ic) * ( denom2 - denom ) ! softcore potential is part of van der Waals energy
      ! -- ti decomp
      if(decpr .and. idecomp > 0) call decpair(3,i,j,foureps(ic)*(denom2 - denom)/(nstlim/ntpr))

      sc_dvdl = sc_dvdl + foureps(ic) * ( -2.0d0 * scalpha * denom3 + scalpha * denom2 )
      ! -- ti decomp
      if(decpr .and. idecomp > 0) call decpair(3,i,j,weight0*foureps(ic)*(2.0d0*scalpha*denom3 - scalpha*denom2)/(nstlim/ntpr))

      df = df + foureps(ic) * ( 12.0d0 * rfour * sigma6(ic) * denom3 - 6.0d0 * rfour * sigma6(ic) * denom2 )
   end if

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
