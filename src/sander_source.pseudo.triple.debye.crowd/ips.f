! <compile=optimized>
#include "dprec.h"

!Performance updates by Ross Walker (TSRI, 2005)

      SUBROUTINE IPSSYS(NATOM,NTYPE,NTB,CG,CUT,CN1,CN2,IAC,X)
!-----------------------------------------------------------------------
!    Calculate system energy and force in 3D periodic systems
!
      use constants, only : FOURPI, third, sixth, INVSQRT2, SQRT2, eighth, half, one, zero
      implicit none

#include "extra.h"
#include "sgld.h"

      _REAL_ DIPSEC,DIPSVAC,DIPSVCC
      INTEGER NTYPE,I3,J3
      LOGICAL LVDW,LELEC
      INTEGER NATOM,NATC,NTB
      _REAL_ CG(*),CN1(*), CN2(*)
      _REAL_ X(*)
      INTEGER IAC(*)
      _REAL_  CUT
      _REAL_ ENBIJ,EELIJ,ECONTI
      INTEGER I,J,NSUMIT(500)
      INTEGER ITI,ITJ,ITMAX,IC,NITI,NITJ
      _REAL_  XI,YI,ZI,XIJ,YIJ,ZIJ,R2
      _REAL_  CGF,CGI,CGIJ,AIJ,CIJ,ANIJ
      _REAL_  SIG2,SIG6,SIG12,ROOT2
      _REAL_  PE0,PE1,PUE0,PUE1,PHE0,PHE1,DPUE,DPHE
      _REAL_  PVC0,PVC1,PUVC0,PUVC1,PHVC0,PHVC1,DPUVC,DPHVC
      _REAL_  PVA0,PVA1,PUVA0,PUVA1,PHVA0,PHVA1,DPUVA,DPHVA
      _REAL_  DUIJ,DHIJ,IPSX1,IPSY1,IPSZ1,BOXLRC(3)
      _REAL_  AIJSUM,CIJSUM,CGSUM,CGIJSUM, RIPS6

!  IPS Radius:

      RIPS2=cut
      RIPS6=RIPS2*RIPS2*RIPS2
      oneRIPS6=one/RIPS6
      oneRIPS12=oneRIPS6*oneRIPS6
      RIPSinv=one/SQRT(RIPS2)
      rips = rips2*ripsinv
      rips2inv = ripsinv*ripsinv

!  Ele IPS parameters:

      BIPSE1=1.109466D0
      BIPSE2=0.0708D0
      BIPSE3=0.0079175D0
      BIPSE0=1.0D0-3.0D0*BIPSE1-5.0D0*BIPSE2-7.0D0*BIPSE3
      !AIPSE=0.0D0*SQRT(2.0D0)-BIPSE0
      AIPSE=-BIPSE0

!  Dispersion IPS parameters:

      BIPSVC1=-1.70554D0
      BIPSVC2=-0.469571D0
      BIPSVC3=0.0000D0
      BIPSVC0=1.0D0-(4.0D0*BIPSVC1+5.0D0*BIPSVC2+6.0D0*BIPSVC3)*third
      AIPSVC=8.0D0*0.4376630D0-BIPSVC0

!  Repulsion IPS parameters:

      BIPSVA1=2.664085D0
      BIPSVA2=-0.611616D0
      BIPSVA3=-0.776091D0
      BIPSVA0=1.0D0-(7.0D0*BIPSVA1+8.0D0*BIPSVA2+9.0D0*BIPSVA3)*sixth
      AIPSVA=64.0D0*6.353604D-3-BIPSVA0

! Energy and force constants:

      PIPSEC=1.0D0+BIPSE0+BIPSE1+BIPSE2+BIPSE3
      DIPSEC=2.0D0*BIPSE1+4.0D0*BIPSE2+6.0D0*BIPSE3
      PIPSVAC=1.0D0+BIPSVA0+BIPSVA1+BIPSVA2+BIPSVA3
      DIPSVAC=2.0D0*BIPSVA1+4.0D0*BIPSVA2+6.0D0*BIPSVA3
      PIPSVCC=1.0D0+BIPSVC0+BIPSVC1+BIPSVC2+BIPSVC3
      DIPSVCC=2.0D0*BIPSVC1+4.0D0*BIPSVC2+6.0D0*BIPSVC3
      PIPSVA0=BIPSVA0/64.0D0-PIPSVAC
      PIPSVC0=BIPSVC0*eighth-PIPSVCC
      PIPSE0=BIPSE0*INVSQRT2-PIPSEC
!
      EIPSSNB=0.0D0
      EIPSSEL=0.0D0

!=======================================================================
!   Main loop begin
!=======================================================================

      DO ITI=1,NTYPE
        NSUMIT(ITI)=0
      ENDDO
      CGSUM=0.0D0
      ITMAX=0
      DO I=1,NATOM
        CGSUM=CGSUM+CG(I)
        ITI=IAC(I)
        IF(ITI.GT.NTYPE)STOP "problem in IAC!"
        NSUMIT(ITI)=NSUMIT(ITI)+1
      ENDDO
      CIJSUM=0.0D0
      AIJSUM=0.0D0
      NNBIPST=0

!  system energy is calculated based on all pairs:
      DO ITI=1,NTYPE
        IC=ITI*(ITI-1)/2+ITI
        AIJ=CN1(IC)
        CIJ=CN2(IC)
        NITI=NSUMIT(ITI)
        ANIJ=NITI*NITI*half
        CIJSUM=CIJSUM+CIJ*ANIJ
        AIJSUM=AIJSUM+AIJ*ANIJ
        NNBIPST=NNBIPST+NITI*NITI
        DO ITJ=ITI+1,NTYPE
          NITJ=NSUMIT(ITJ)
          IC=ITJ*(ITJ-1)/2+ITI
          AIJ=CN1(IC)
          CIJ=CN2(IC)
          ANIJ=NITI*NITJ
          CIJSUM=CIJSUM+CIJ*ANIJ
          AIJSUM=AIJSUM+AIJ*ANIJ
          NNBIPST=NNBIPST+2*NITI*NITJ
        ENDDO
      ENDDO
      CGIJSUM=CGSUM*CGSUM*half
      if( cgijsum .lt. 1.D-15 ) cgijsum = 0.d0
      EIPSSNB=AIJSUM*(PIPSVAC+AIPSVA/64.0D0)*oneRIPS12  &
                -CIJSUM*(PIPSVCC+AIPSVC*eighth)*oneRIPS6
      EIPSSEL=CGIJSUM*(PIPSEC+AIPSE*invsqrt2)*RIPSINV
      IF(.NOT.TEIPS)EIPSSEL=0.0D0
      IF(.NOT.TVIPS)EIPSSNB=0.0D0

!=======================================================================
!   Main loop end
!=======================================================================

! Calculate volume virials:

      VIRIPS=-(EIPSSNB+EIPSSEL)
      VBOXIPS=FOURPI*RIPS2*RIPS*third
      IF(NTB.EQ.0)THEN

        ! For non-periodic system, system energy is calculated based on cutoff
        NNBIPS=0
        DO I=1,NATOM
          ! Atom i long-range reference and self-interaction:
          I3=3*I-2
          XI=X(I3)
          YI=X(I3+1)
          ZI=X(I3+2)
          NNBIPS=NNBIPS+1
          DO J=I+1,NATOM
            J3=3*J-2
            XIJ=XI-X(J3)
            YIJ=YI-X(J3+1)
            ZIJ=ZI-X(J3+2)
            R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
            !  RIPS2 is infinite except for non periodic conditions
            IF(R2.LT.RIPS2)NNBIPS=NNBIPS+2
          ENDDO
        ENDDO
      ENDIF

      if( master ) then
         WRITE(6,'(" IPS Radius: ",F6.2," A")')RIPS
         IF(TEIPS)THEN
            WRITE(6,'(" IPS parameters for electrostatic energy:")')
            WRITE(6,'("  AIPSE,BIPSE0,BIPSE1,BIPSE2,BIPSE3="/5F10.6/)') &
              AIPSE,BIPSE0,BIPSE1,BIPSE2,BIPSE3
         ENDIF
         IF(TVIPS)THEN
            WRITE(6,'(" IPS parameters for L-J energy:")')
            WRITE(6,'("  AIPSVC,BIPSVC0,BIPSVC1,BIPSVC2,BIPSVC3="/5F10.6/)') &
	            AIPSVC,BIPSVC0,BIPSVC1,BIPSVC2,BIPSVC3
	        WRITE(6,'("  AIPSVA,BIPSVA0,BIPSVA1,BIPSVA2,BIPSVA3="/5F10.6/)') &
	         AIPSVA,BIPSVA0,BIPSVA1,BIPSVA2,BIPSVA3
	     ENDIF
         WRITE(6,'(" EIPSSNB, EIPSSEL=",2E20.7)') EIPSSNB,EIPSSEL
         IF(NTB.GT.0)THEN
            WRITE(6,'("  IPS region volume=",E20.7/)') VBOXIPS
         ELSE
           WRITE(6,'("  Enclosing atom pairs/total pairs=",I10, "/",I10/)') &
               NNBIPS,NNBIPST
         ENDIF
      else 
        EIPSSNB=0.0D0
        EIPSSEL=0.0D0
        VIRIPS=0.0D0
        NNBIPS=0
      end if
      RETURN 

end subroutine ipssys

SUBROUTINE IPSUPDATE(NTB)

!-----------------------------------------------------------------------
!     Update parameters once IPS radius or the box size changed
!-----------------------------------------------------------------------

      use nblist, only : volume
      implicit none
#ifdef MPI
#  include "parallel.h"
#  include "mpif.h"
   integer ierr
#endif
#include "sgld.h"
      INTEGER NTB
      INTEGER NNBTMP
      _REAL_ FIPS,CHANGE 

      IF(NTB.GT.0)THEN
        FIPS=VBOXIPS/volume
      ELSE
#ifdef MPI
        call mpi_allreduce(NNBIPS,NNBTMP,1,MPI_INTEGER, &
         mpi_sum,commsander,ierr)
        NNBIPS=NNBTMP 
#endif
        FIPS=NNBIPS*1.0D0/NNBIPST
      ENDIF
      CHANGE=FIPS-1.0D0
      CHANGE=CHANGE*CHANGE

      ! Update system energies and forces:
      IF(CHANGE.GT.1.0D-8)THEN
         VIRIPS=VIRIPS*FIPS
         EIPSSNB=EIPSSNB*FIPS
         EIPSSEL=EIPSSEL*FIPS
      ENDIF
      VBOXIPS=volume
      NNBIPST=NNBIPS
      RETURN
end subroutine ipsupdate

SUBROUTINE EEXIPS(ENB,EEL,IFRSTA,ILASTA,NTB,IAC,numex,natex,CG, &
                        CN1,CN2,DX,X,VIR)

!-----------------------------------------------------------------------
!   3D IPS interaction between excluded atom pairs
!   This routine must be called first to update IPS parameters when needed
!       and to initalize electrostatic and vdw energies
!
!   by Xiongwu Wu  - 9/10/2004
!
!-----------------------------------------------------------------------

      use constants, only : zero
      implicit none
#include "sgld.h"
      _REAL_ ENB, EEL
      INTEGER IFRSTA,ILASTA, NTB,IAC(*),numex(*),natex(*)
      _REAL_ CG(*),CN1(*), CN2(*),VIR(3,3)
      _REAL_ DX(*),X(*)
      ! local:
      INTEGER I,J,K,IC,ICX,ITI,ITJ
      INTEGER I3,I31,I32,J3,J31,J32
      INTEGER  ILAST, IFIRST,LJROW
      _REAL_ ESCALE,VSCALE
      _REAL_ DXI, DYI, DZI,DIJ,DIJX,DIJY,DIJZ
      _REAL_ ENBIJ,EELIJ,ECONTI
      _REAL_  XI,YI,ZI,XIJ,YIJ,ZIJ,x2,y2,z2
      _REAL_  CGF,CGI,CGIJ,AIJ,CIJ,AEXIJ,CEXIJ,SIG2,SIG6,SIG12
      _REAL_  R2,R6,R12,U2,TWOU1,TWOU2,oneTWOU2,TWOU6,TWOU12
      _REAL_  PE,PEU,DEU,PVA,PVAU,DVAU,PVC,PVCU,DVCU

      ! check to see if volume or atom pair changed 
      CALL IPSUPDATE(NTB)

      ! Setup constants for use in inner loops
      ENB=EIPSSNB
      EEL=EIPSSEL
      DEU=0.0D0
      DVAU=0.0D0
      DVCU=0.0D0
      CGIJ=0.0D0
      AIJ=0.0D0
      CIJ=0.0D0
      vir(1:3,1:3) = zero

!=======================================================================
!   Main loop begin
!=======================================================================

      ILAST = 0
      DO I=1,IFRSTA-1
         ILAST=ILAST+numex(I)
      ENDDO
      NNBIPS=0
      DO I=IFRSTA,ILASTA !1,natom
         NNBIPS=NNBIPS+1
         IFIRST = ILAST+1
         ILAST=ILAST+numex(I)
         IF(TEIPS)THEN
            CGI=CG(I)
            CGIJ=CGI*CGI
            EELIJ=0.5D0*CGIJ*PIPSE0*RIPSINV
            EEL=EEL+EELIJ
         ENDIF
         IF(TVIPS)THEN
            ITI=IAC(I)
            IC=ITI*(ITI-1)/2+ITI
            AIJ=CN1(IC)
            CIJ=CN2(IC)
            ! Atom i long-range reference and self-interaction
            ENBIJ=0.5D0*(AIJ*PIPSVA0*oneRIPS6-CIJ*PIPSVC0)*oneRIPS6
            ENB=ENB+ENBIJ
         ENDIF
         I32=3*I
         I31=I32-1
         I3=I32-2
         DXI=DX(I3)
         DYI=DX(I31)
         DZI=DX(I32)
         XI=X(I3)
         YI=X(I31)
         ZI=X(I32)
         DO K=IFIRST,ILAST !1,No interactions for this atom.
            J = natex(K) !negative if atom J should be excluded from interaction
                         ! with I.  This deals with 1-2,1-3 and 1-4 only. 
                         ! Does not deal with cutoff.
           IF(J.LE.0) cycle
           J32=3*J
           J31=J32-1
           J3=J32-2
           XIJ=X(J3)-XI
           YIJ=X(J31)-YI
           ZIJ=X(J32)-ZI
           R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
           NNBIPS=NNBIPS+2
           R6=R2*R2*R2
           R12=R6*R6
           U2=R2*RIPS2INV
           TWOU2=2.0D0-U2
           oneTWOU2=1.0d0/TWOU2
           IF(TEIPS)THEN
              CGIJ=CGI*CG(J)
              TWOU1=SQRT(oneTWOU2)
              PE=BIPSE0+U2*(BIPSE1+U2*(BIPSE2+U2*BIPSE3))
              PEU=2.0D0*(BIPSE1+U2*(2.0D0*BIPSE2+3.0D0*BIPSE3*U2))
              DEU=(PEU+PE*oneTWOU2)*TWOU1
              EELIJ=CGIJ*(PE*TWOU1-PIPSEC)*RIPSINV
              EEL=EEL+EELIJ
           ENDIF
           IF(TVIPS)THEN
              ITJ=IAC(J)
              IF(ITI.GT.ITJ)THEN
                 IC=ITI*(ITI-1)/2+ITJ
              ELSE
                 IC=ITJ*(ITJ-1)/2+ITI
              ENDIF
              AIJ=CN1(IC)
              CIJ=CN2(IC)
              TWOU6=oneTWOU2*oneTWOU2*oneTWOU2
              TWOU12=TWOU6*TWOU6

              !  L-J r6 term:

              PVC=BIPSVC0+U2*(BIPSVC1+U2*(BIPSVC2+U2*BIPSVC3))
              PVCU=2.0D0*(BIPSVC1+U2*(2.0D0*BIPSVC2+3.0D0*BIPSVC3*U2))
              DVCU=(PVCU+6.0D0*PVC*oneTWOU2)*TWOU6
 
              !  L-J r12 term:

              PVA=BIPSVA0+U2*(BIPSVA1+U2*(BIPSVA2+U2*BIPSVA3))
              PVAU=2.0D0*(BIPSVA1+U2*(2.0D0*BIPSVA2+3.0D0*BIPSVA3*U2))
              DVAU=(PVAU+12.0D0*PVA*oneTWOU2)*TWOU12
              ENBIJ=AIJ*(PVA*TWOU12-PIPSVAC)*oneRIPS12  &
                  -CIJ*(PVC*TWOU6-PIPSVCC)*oneRIPS6

              ENB=ENB+ENBIJ
              DIJ=-(CGIJ*DEU*RIPSINV +(AIJ*DVAU*oneRIPS6  &
                   -CIJ*DVCU)*oneRIPS6)*RIPS2INV
           ELSE
              DIJ=-CGIJ*DEU*RIPSINV*RIPS2INV
           ENDIF

           DIJX=DIJ*XIJ
           DIJY=DIJ*YIJ
           DIJZ=DIJ*ZIJ
           DXI=DXI-DIJX
           DYI=DYI-DIJY
           DZI=DZI-DIJZ
           DX(J3)=DX(J3)+DIJX
           DX(J31)=DX(J31)+DIJY
           DX(J32)=DX(J32)+DIJZ
           vir(1,1) = vir(1,1) - DIJX*XIJ
           vir(1,2) = vir(1,2) - DIJX*YIJ
           vir(1,3) = vir(1,3) - DIJX*ZIJ
           vir(2,1) = vir(2,1) - DIJY*XIJ
           vir(2,2) = vir(2,2) - DIJY*YIJ
           vir(2,3) = vir(2,3) - DIJY*ZIJ
           vir(3,1) = vir(3,1) - DIJZ*XIJ
           vir(3,2) = vir(3,2) - DIJZ*YIJ
           vir(3,3) = vir(3,3) - DIJZ*ZIJ
         end do
         DX(I3) =  DXI
         DX(I31) =  DYI
         DX(I32) =  DZI
      ENDDO

!=======================================================================
!   Main loop end
!=======================================================================

      RETURN
end subroutine eexips

