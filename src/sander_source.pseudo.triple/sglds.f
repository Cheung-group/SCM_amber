! <compile=optimized>
#include "dprec.h"

    SUBROUTINE PSGLD(NATOM,AMASS,V,VSG)
!-----------------------------------------------------------------------
!     This routine perform initiation for the self-guided        
!       Langevin Dynamcs (SGLD) simulaiton                 
!
      implicit none
#include "md.h"
#ifdef MPI
#  include "parallel.h"
#endif
#include "sgld.h"
      INTEGER NATOM,I,I3,M
      _REAL_ AMASS(*),V(*),VSG(*)
      _REAL_ AMASSI,VI3
!
!  Check for invalid sgld setting
      IF(ISGSTA < 1)ISGSTA=1
      IF(ISGEND > NATOM .OR. ISGEND < 1)ISGEND=NATOM
      IF(TSGAVG.LT.DT)TSGAVG=DT
      SGAVG1=DT/TSGAVG
      SGAVG0=1.0D0-SGAVG1
#ifdef MPI
      if(mytaskid.eq.0)THEN
#endif
        WRITE(6,910)ISGSTA,ISGEND,TSGAVG
        IF(TEMPSG > 0.0D0)THEN
          WRITE(6,920)TEMPSG
        ELSE
          WRITE(6,925)SGFT
        ENDIF
        IF(TLANGV)THEN
          WRITE(6,930)GAMMA_LN
        ELSE
          WRITE(6,935)
        ENDIF
#ifdef MPI
      ENDIF
#endif
      GAMMAS=GAMMA_LN/20.455d0
!    Initialize arrays
      AVGGG=0.0D0
      DO I=1,NATOM
        I3=3*I-3
        AMASSI = AMASS(I)
        IF((I-ISGSTA)*(ISGEND-I).LT.0)AMASSI=0.0D0
        DO M=1,3
          VI3=V(I3+M)
          VSG(I3+M)=VI3
          AVGGG=AVGGG+AMASSI*VI3*VI3
        END DO
      END DO
!   AVGGG=gammas*(tsg-dt)*<m*vsg*vsg>
      GAMMAT=TSGAVG-DT
      IF(GAMMAS>0.0D0)GAMMAT=GAMMA_LN*GAMMAT/(1.0D0+GAMMA_LN*DT)
      AVGGG=AVGGG*GAMMAT
910   FORMAT("************************************"/  &
      "Parameters for self-guided Langevin dynamics (SGLD) simulation"//  &
          "  Guiding range from ",I5,"  TO ",I5 /  &
          "  Averaging over ",F10.4," ps ")
920   FORMAT("  Guiding temperature is set at:",F8.2, " K" )
925   FORMAT("  Guiding factor is set at : ",F6.4 )
930   FORMAT("  friction constant:",F8.2," 1/ps" / &
          "*****************************************"/)
935   FORMAT("  SGLD method is applied to MD simulation" / &
          "*****************************************"/)
      RETURN
      END



      SUBROUTINE SGLDW(NATOM,ISTART,IEND, &
             DTX,TEMP0,RNDF,AMASS,WINV,F,V,VSG)
!-----------------------------------------------------------------------
!     This routine perform SGLD integration        
!
      implicit none
#include "sgld.h"
#ifdef MPI
# include "mpif.h"
      integer ierr
# include "parallel.h"
      _REAL_ temp1(3)
# ifndef USE_MPI_IN_PLACE
      _REAL_ :: temp2(3)
# endif
#endif
      INTEGER NATOM,ISTART,IEND
      _REAL_ DTX,TEMP0,RNDF
      _REAL_ AMASS(*),WINV(*),F(*),V(*),VSG(*)
!
      INTEGER I,M,I3
      _REAL_ BOLTZ,AMASSI,VSGI,ACMGG
      _REAL_ FACT,WFAC,VI3,FI3,DRAGI,PV1,PV2,RSD,FLN,SGSCAL
      PARAMETER (BOLTZ = 1.987192d-3)
!
!   Estimate guiding factor or guiding temperature
        IF(TEMPSG > 0.0D0)SGFT=RNDF*BOLTZ*TEMPSG/AVGGG
        TEMPSGI=SGFT*AVGGG/RNDF/BOLTZ
        PV1=0.0d0 
        PV2=0.0d0
        ACMGG=0.0d0
        DO  I = 1,NATOM 
          IF((I-ISTART)*(IEND-I).LT.0)THEN
!   Keep random number series the same as that in a single cpu simulation
            CALL GAUSS( 0.D0, 1.0D0, FLN )
            CALL GAUSS( 0.D0, 1.0D0, FLN )
            CALL GAUSS( 0.D0, 1.0D0, FLN )
            CYCLE
          ENDIF
          AMASSI = AMASS(I)
          WFAC =  DTX*0.5D0*WINV(I)
          RSD = SQRT(2.D0*GAMMAS*BOLTZ*TEMP0*AMASSI/DTX)
          IF((I-ISGSTA)*(ISGEND-I).LT.0)AMASSI=0.0D0
          DO  M = 1,3
            I3 = 3*(I-1)+M
            VI3=V(I3)
            VSGI=SGAVG0*VSG(I3)+SGAVG1*VI3
            ACMGG=ACMGG+AMASSI*VSG(I3)*VSGI
            VSG(I3)=VSGI
            DRAGI=SGFT*GAMMAS*AMASSI*VSGI
            CALL GAUSS( 0.D0, RSD, FLN )
            FI3=F(I3)+fln+DRAGI
            F(I3)=FI3
            VI3 = VI3 + FI3*wfac
            PV1=PV1+VI3*DRAGI
            PV2=PV2+AMASSI*VI3*VI3
          END DO
        END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=PV1
          TEMP1(2)=PV2
          TEMP1(3)=ACMGG
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,3,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          PV1=TEMP1(1)
          PV2=TEMP1(2)
          ACMGG=TEMP1(3)
#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,3, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          PV1=TEMP2(1)
          PV2=TEMP2(2)
          ACMGG=TEMP2(3)
# endif
        ENDIF
#endif
        AVGGG=SGAVG0*AVGGG+SGAVG1*GAMMAT*ACMGG
        SGSCAL=(1.0D0+0.5D0*GAMMAS*DTX)*PV1/(PV2-0.5D0*DTX*PV1)
        DO  I = ISTART,IEND 
          WFAC = WINV(I)*DTX
          IF((I-ISGSTA)*(ISGEND-I).GE.0)THEN
            FACT=0.5D0*DTX*(GAMMAS+SGSCAL)
          ELSE
            FACT=0.5D0*DTX*GAMMAS
          ENDIF
          DO  M = 1,3
            I3 = 3*(I-1)+M
            V(I3)=((1.0D0-FACT)*V(I3)+F(I3)*WFAC)/(1.0D0+FACT)
          END DO
        END DO
        RETURN
        END


        SUBROUTINE SGMDW(NATOM,ISTART,IEND, &
             DTX,RNDF,AMASS,WINV,F,V,VSG)
!-----------------------------------------------------------------------
!     This routine calculate guiding force using SGLD method 
!     for MD simulation        
!
      implicit none
#include "sgld.h"
#ifdef MPI
# include "mpif.h"
      integer ierr
# include "parallel.h"
      _REAL_ temp1(3)
# ifndef USE_MPI_IN_PLACE
      _REAL_ :: temp2(3)
# endif
#endif
      INTEGER NATOM,ISTART,IEND
      _REAL_ DTX,RNDF
      _REAL_ AMASS(*),WINV(*),F(*),V(*),VSG(*)
!
      INTEGER JSTA,JEND,I,M,I3
      _REAL_ BOLTZ,AMASSI,VSGI,ACMGG
      _REAL_ FACT,WFAC,VI3,FI3,DRAGI,PV1,PV2,SGSCAL
      PARAMETER (BOLTZ = 1.987192d-3)
!
        JSTA=ISTART
        JEND=IEND
        IF(JSTA.LT.ISGSTA)JSTA=ISGSTA
        IF(JEND.GT.ISGEND)JEND=ISGEND
!    Estimate guiding factor or guiding temperature
        IF(TEMPSG > 0.0D0)SGFT=RNDF*BOLTZ*TEMPSG/AVGGG
        TEMPSGI=SGFT*AVGGG/RNDF/BOLTZ
        PV1=0.0d0 
        PV2=0.0d0
        ACMGG=0.0d0
        DO  I = JSTA,JEND 
          AMASSI = AMASS(I)
          wfac =  DTX*0.5d0*WINV(I)
          DO  M = 1,3
            I3 = 3*(I-1)+M
            VI3=V(I3)
            VSGI=SGAVG0*VSG(I3)+SGAVG1*VI3
            ACMGG=ACMGG+AMASSI*VSG(I3)*VSGI
            DRAGI=SGFT*AMASSI*VSGI/20.455d0
            VSG(I3)=VSGI
            FI3=F(I3)+DRAGI
            F(I3)=FI3
            VI3 = VI3 + FI3*wfac
            PV1=PV1+VI3*DRAGI
            PV2=PV2+AMASSI*VI3*VI3
          END DO
        END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=PV1
          TEMP1(2)=PV2
          TEMP1(3)=ACMGG
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,3,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          PV1=TEMP1(1)
          PV2=TEMP1(2)
          ACMGG=TEMP1(3)
# else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,3, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          PV1=TEMP2(1)
          PV2=TEMP2(2)
          ACMGG=TEMP2(3)
# endif
        ENDIF
#endif
        AVGGG=SGAVG0*AVGGG+SGAVG1*GAMMAT*ACMGG
        SGSCAL=PV1/(PV2-0.5D0*DTX*PV1)
        DO  I = JSTA,JEND 
          WFAC = 0.5D0*WINV(I)*DTX
          FACT=AMASS(I)*SGSCAL/(1.0D0+0.5D0*SGSCAL*DTX)
          DO  M = 1,3
              I3 = 3*(I-1)+M
              VI3 = V(I3) + F(I3)*WFAC
              F(I3)=F(I3)-FACT*VI3
          END DO
        END DO
        RETURN
        END


