! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
subroutine qm2_rotate_qmqm(loop_count,IQM, JQM,NI,NJ,XI,XJ,W,KR,E1B,E2A, &
                       ENUC,qmitype,qmjtype)

!********************************************************
! Current routine maintained by: Ross Walker (TSRI, 2005)
! Inlining and Optimising by: Ross Walker (TSRI, 2005)
!********************************************************

!********************************************************
!   ROTATE CALCULATES THE TWO-PARTICLE INTERACTIONS.
!
!   ON INPUT
!             IQM    = qm atom I number (in 1 to nquant loop)
!             JQM    = qm atom J number in inner loop
!             NI     = ATOMIC NUMBER OF FIRST ATOM.
!             NJ     = ATOMIC NUMBER OF SECOND ATOM.
!             XI     = COORDINATE OF FIRST ATOM.
!             XJ     = COORDINATE OF SECOND ATOM.
!
! ON OUTPUT W      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS.
!           E1B,E2A= ARRAY OF ELECTRON-NUCLEAR ATTRACTION INTEGRALS,
!                    E1B = ELECTRON ON ATOM NI ATTRACTING NUCLEUS OF NJ.
!           ENUC   = NUCLEAR-NUCLEAR REPULSION TERM.
!
! *** THIS ROUTINE COMPUTES THE REPULSION AND NUCLEAR ATTRACTION
!     INTEGRALS OVER MOLECULAR-FRAME COORDINATES.  THE INTEGRALS OVER
!     LOCAL FRAME COORDINATES ARE EVALUATED BY SUBROUTINE REPP AND
!     STORED AS FOLLOWS (WHERE P-SIGMA = O,   AND P-PI = P AND P* )
!     IN RI
!     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!     (PP/P*P*)=21,   (P*P/P*P)=22.
!
!***********************************************************************
      use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, &
                              qm2_rij_eqns, EXPONENTIAL_CUTOFF
      use constants, only : one, A_TO_BOHRS, A2_TO_BOHRS2
      implicit none
!Passed in
      integer, intent(in) :: loop_count, iqm, jqm, ni, nj,qmitype,qmjtype
      integer, intent(inout) :: kr
      _REAL_, intent(in) :: xi(3), xj(3)
      _REAL_, intent(out) :: W(100), e1b(10), E2A(10), enuc

!Local
      _REAL_ :: X(3),Y(3),Z(3), RI(22)
      _REAL_ :: temp_real, temp_real2, anam1, C1, oneRIJ, RIJ
      _REAL_ :: rr,rr2, exp1i, exp1j, sqrtaee, BDD1i, BDD1j, bdd1ij
      _REAL_ :: a, xx11, xx21, xx22, xx31, xx32, xx33, yy11, yy21, yy22
      _REAL_ :: zz11, zz21, zz22, zz31, zz32, zz33, yyzz11, yyzz21, yyzz22
      _REAL_ :: xy11, xy21, xy22, xy31, xy32, xz11, xz21, xz22, xz31, xz32, xz33
      _REAL_ :: YZ11, yz21, yz22, yz31, yz32, css1, css2, csp1, cpps1, cppp1, csp2, cpps2, cppp2, scale
      _REAL_ :: PDDG_EXP1, PDDG_EXP2, PDDG_EXP3, PDDG_EXP4, PDDG_CORR
      integer :: i, ki
      logical :: SI,SJ

      X(1)=XI(1)-XJ(1)
      X(2)=XI(2)-XJ(2)
      X(3)=XI(3)-XJ(3)
#include "qm2_array_locations.h"
      if (qmmm_nml%qmqmrij_incore) then
        !We already have RIJ info in memory
        rr2=qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS2,loop_count)
        oneRIJ=qm2_rij_eqns%qmqmrijdata(QMQMONERIJ,loop_count)
        RIJ=qm2_rij_eqns%qmqmrijdata(QMQMRIJ,loop_count)
        rr=qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS,loop_count)
        exp1i=qm2_rij_eqns%qmqmrijdata(QMQMEXP1I,loop_count)
        exp1j=qm2_rij_eqns%qmqmrijdata(QMQMEXP1J,loop_count)
        sqrtaee=qm2_rij_eqns%qmqmrijdata(QMQMSQRTAEE,loop_count)
        if (qmmm_struct%PDDG_IN_USE) then
          PDDG_EXP1 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP1,loop_count)
          PDDG_EXP2 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP2,loop_count)
          PDDG_EXP3 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP3,loop_count)
          PDDG_EXP4 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP4,loop_count)
        end if
      else
        RIJ=X(1)*X(1)+X(2)*X(2)+X(3)*X(3) 
        rr2=RIJ*A2_TO_BOHRS2
        oneRIJ = one/SQRT(RIJ)
        RIJ=RIJ*oneRIJ!one/oneRIJ
        rr=RIJ*A_TO_BOHRS
        exp1i = EXP(-qm2_params%cc_exp_params(iqm)*RIJ)
        exp1j = EXP(-qm2_params%cc_exp_params(jqm)*RIJ)
        BDD1i=qm2_params%multip_2c_elec_params(3,iqm)
        BDD1j=qm2_params%multip_2c_elec_params(3,jqm)
        BDD1ij=bdd1i+bdd1j
        BDD1ij=bdd1ij*bdd1ij
        SQRTAEE=1.0d0/sqrt(RR2+bdd1ij)
        if (qmmm_struct%PDDG_IN_USE) then
!PDDG Specific terms
          PDDG_EXP1 = EXP(-10.0D0 * (RIJ - qm2_params%pddge1(iqm) - qm2_params%pddge1(jqm))**2)
          PDDG_EXP2 = EXP(-10.0D0 * (RIJ - qm2_params%pddge1(iqm) - qm2_params%pddge2(jqm))**2)
          PDDG_EXP3 = EXP(-10.0D0 * (RIJ - qm2_params%pddge2(iqm) - qm2_params%pddge1(jqm))**2)
          PDDG_EXP4 = EXP(-10.0D0 * (RIJ - qm2_params%pddge2(iqm) - qm2_params%pddge2(jqm))**2)
        end if
      end if
      X(1) = X(1)*oneRIJ
      X(2) = X(2)*oneRIJ
      X(3) = X(3)*oneRIJ

      CALL qm2_repp(iqm,jqm,rr,rr2,RI,SQRTAEE)
      IF(qmmm_nml%qmqm_erep_incore) then
! Ross Walker - We will need these repulsion integrals later to
!               calculate qm-qm analytical derivatives so store them
!               in our array.
        qm2_struct%qm_qm_e_repul(1:22,loop_count) = RI(1:22)
      end if

      IF (ABS(X(3)).GT.0.99999999D0) THEN
         X(3) = SIGN(1.D0,X(3))
         Y(1) = 0.D0
         Y(2) = 1.D0
         Y(3) = 0.D0
         Z(1) = 1.D0
         Z(2) = 0.D0
         Z(3) = 0.D0
      ELSE
         Z(3)=SQRT(1.D0-X(3)*X(3))
         A=1.D0/Z(3)
         Y(1)=-A*X(2)*SIGN(1.D0,X(1))
         Y(2)=ABS(A*X(1))
         Y(3)=0.D0
         Z(1)=-A*X(1)*X(3)
         Z(2)=-A*X(2)*X(3)
      ENDIF
      SI = (qm2_params%natomic_orbs(iqm) > 1)
      SJ = (qm2_params%natomic_orbs(jqm) > 1)
      IF ( SI .OR. SJ) THEN
         XX11 = X(1)*X(1)
         XX21 = X(2)*X(1)
         XX22 = X(2)*X(2)
         XX31 = X(3)*X(1)
         XX32 = X(3)*X(2)
         XX33 = X(3)*X(3)
         YY11 = Y(1)*Y(1)
         YY21 = Y(2)*Y(1)
         YY22 = Y(2)*Y(2)
         ZZ11 = Z(1)*Z(1)
         ZZ21 = Z(2)*Z(1)
         ZZ22 = Z(2)*Z(2)
         ZZ31 = Z(3)*Z(1)
         ZZ32 = Z(3)*Z(2)
         ZZ33 = Z(3)*Z(3)
         YYZZ11 = YY11+ZZ11
         YYZZ21 = YY21+ZZ21
         YYZZ22 = YY22+ZZ22
         XY11 = 2.D0*X(1)*Y(1)
         XY21 =      X(1)*Y(2)+X(2)*Y(1)
         XY22 = 2.D0*X(2)*Y(2)
         XY31 =      X(3)*Y(1)
         XY32 =      X(3)*Y(2)
         XZ11 = 2.D0*X(1)*Z(1)
         XZ21 =      X(1)*Z(2)+X(2)*Z(1)
         XZ22 = 2.D0*X(2)*Z(2)
         XZ31 =      X(1)*Z(3)+X(3)*Z(1)
         XZ32 =      X(2)*Z(3)+X(3)*Z(2)
         XZ33 = 2.D0*X(3)*Z(3)
         YZ11 = 2.D0*Y(1)*Z(1)
         YZ21 =      Y(1)*Z(2)+Y(2)*Z(1)
         YZ22 = 2.D0*Y(2)*Z(2)
         YZ31 =      Y(1)*Z(3)
         YZ32 =      Y(2)*Z(3)
      ENDIF
!     (S S/S S)
      W(1)=RI(1)
      KI = 1
      IF (SJ) THEN
!     (S S/PX S)
         W(2)=RI(5)*X(1)
!     (S S/PX PX)
         W(3)=RI(11)*XX11+RI(12)*YYZZ11
!     (S S/PY S)
         W(4)=RI(5)*X(2)
!     (S S/PY PX)
         W(5)=RI(11)*XX21+RI(12)*YYZZ21
!     (S S/PY PY)
         W(6)=RI(11)*XX22+RI(12)*YYZZ22
!     (S S/PZ S)
         W(7)=RI(5)*X(3)
!     (S S/PZ PX)
         W(8)=RI(11)*XX31+RI(12)*ZZ31
!     (S S/PZ PY)
         W(9)=RI(11)*XX32+RI(12)*ZZ32
!     (S S/PZ PZ)
         W(10)=RI(11)*XX33+RI(12)*ZZ33
         KI = 10
      ENDIF
      IF (SI) THEN
         IF (SJ) THEN
!     (PX S/S S)
            W(11)=RI(2)*X(1)
!     (PX S/PX S)
            W(12)=RI(6)*XX11+RI(7)*YYZZ11
!     (PX S/PX PX)
            W(13)=X(1)*(RI(13)*XX11+RI(14)*YYZZ11)+RI(15)*(Y(1)*XY11+Z(1)*XZ11)
!     (PX S/PY S)
            W(14)=RI(6)*XX21+RI(7)*YYZZ21
!     (PX S/PY PX)
            W(15)=X(1)*(RI(13)*XX21+RI(14)*YYZZ21)+RI(15)*(Y(1)*XY21+Z(1)*XZ21)
!     (PX S/PY PY)
            W(16)=X(1)*(RI(13)*XX22+RI(14)*YYZZ22)+RI(15)*(Y(1)*XY22+Z(1)*XZ22)
!     (PX S/PZ S)
            W(17)=RI(6)*XX31+RI(7)*ZZ31
!     (PX S/PZ PX)
            W(18)=X(1)*(RI(13)*XX31+RI(14)*ZZ31)+RI(15)*(Y(1)*XY31+Z(1)*XZ31)
!     (PX S/PZ PY)
            W(19)=X(1)*(RI(13)*XX32+RI(14)*ZZ32)+RI(15)*(Y(1)*XY32+Z(1)*XZ32)
!     (PX S/PZ PZ)
            W(20)=X(1)*(RI(13)*XX33+RI(14)*ZZ33)+RI(15)*(          Z(1)*XZ33)
!     (PX PX/S S)
            W(21)=RI(3)*XX11+RI(4)*YYZZ11
!     (PX PX/PX S)
            W(22)=X(1)*(RI(8)*XX11+RI(9)*YYZZ11)+RI(10)*(Y(1)*XY11+Z(1)*XZ11)
!     (PX PX/PX PX)
            W(23) =                                             &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX11+RI(18)*XX11*YYZZ11  &
           +RI(19)*(YY11*YY11+ZZ11*ZZ11)                        &
           +RI(20)*(XY11*XY11+XZ11*XZ11)+RI(21)*(YY11*ZZ11+ZZ11*YY11) &
           +RI(22)*YZ11*YZ11
!     (PX PX/PY S)
            W(24)=X(2)*(RI(8)*XX11+RI(9)*YYZZ11)+RI(10)*(Y(2)*XY11+Z(2)*XZ11)
!     (PX PX/PY PX)
            W(25) =                                                   &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX21+RI(18)*XX11*YYZZ21        &
           +RI(19)*(YY11*YY21+ZZ11*ZZ21)+RI(20)*(XY11*XY21+XZ11*XZ21) &
           +RI(21)*(YY11*ZZ21+ZZ11*YY21)+RI(22)*YZ11*YZ21
!     (PX PX/PY PY)
            W(26) =                                                   &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX22+RI(18)*XX11*YYZZ22        &
           +RI(19)*(YY11*YY22+ZZ11*ZZ22)+RI(20)*(XY11*XY22+XZ11*XZ22) &
           +RI(21)*(YY11*ZZ22+ZZ11*YY22)+RI(22)*YZ11*YZ22
!     (PX PX/PZ S)
            W(27)=X(3)*(RI(8)*XX11+RI(9)*YYZZ11)+RI(10)*(         +Z(3)*XZ11)
!     (PX PX/PZ PX)
            W(28) =                                         &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX31                &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ31      &
           +RI(20)*(XY11*XY31+XZ11*XZ31)+RI(22)*YZ11*YZ31
!     (PX PX/PZ PY)
            W(29) =                                         &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX32                &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ32      &
           +RI(20)*(XY11*XY32+XZ11*XZ32)+RI(22)*YZ11*YZ32
!     (PX PX/PZ PZ)
            W(30) =                                         &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX33                &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ33      &
           +RI(20)*XZ11*XZ33
!     (PY S/S S)
            W(31)=RI(2)*X(2)
!     (PY S/PX S)
            W(32)=RI(6)*XX21+RI(7)*YYZZ21
!     (PY S/PX PX)
            W(33)=X(2)*(RI(13)*XX11+RI(14)*YYZZ11)+RI(15)*(Y(2)*XY11+Z(2)*XZ11)
!     (PY S/PY S)
            W(34)=RI(6)*XX22+RI(7)*YYZZ22
!     (PY S/PY PX)
            W(35)=X(2)*(RI(13)*XX21+RI(14)*YYZZ21)+RI(15)*(Y(2)*XY21+Z(2)*XZ21)
!     (PY S/PY PY)
            W(36)=X(2)*(RI(13)*XX22+RI(14)*YYZZ22)+RI(15)*(Y(2)*XY22+Z(2)*XZ22)
!     (PY S/PZ S)
            W(37)=RI(6)*XX32+RI(7)*ZZ32
!     (PY S/PZ PX)
            W(38)=X(2)*(RI(13)*XX31+RI(14)*ZZ31)+RI(15)*(Y(2)*XY31+Z(2)*XZ31)
!     (PY S/PZ PY)
            W(39)=X(2)*(RI(13)*XX32+RI(14)*ZZ32)+RI(15)*(Y(2)*XY32+Z(2)*XZ32)
!     (PY S/PZ PZ)
            W(40)=X(2)*(RI(13)*XX33+RI(14)*ZZ33)+RI(15)*(         +Z(2)*XZ33)
!     (PY PX/S S)
            W(41)=RI(3)*XX21+RI(4)*YYZZ21
!     (PY PX/PX S)
            W(42)=X(1)*(RI(8)*XX21+RI(9)*YYZZ21)+RI(10)*(Y(1)*XY21+Z(1)*XZ21)
!     (PY PX/PX PX)
            W(43) =                                             &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX11+RI(18)*XX21*YYZZ11  &
           +RI(19)*(YY21*YY11+ZZ21*ZZ11)+RI(20)*(XY21*XY11+XZ21*XZ11) &
           +RI(21)*(YY21*ZZ11+ZZ21*YY11)+RI(22)*YZ21*YZ11
!     (PY PX/PY S)
            W(44)=X(2)*(RI(8)*XX21+RI(9)*YYZZ21)+RI(10)*(Y(2)*XY21+Z(2)*XZ21)
!     (PY PX/PY PX)
            W(45) =                                             &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX21+RI(18)*XX21*YYZZ21  &
           +RI(19)*(YY21*YY21+ZZ21*ZZ21)+RI(20)*(XY21*XY21+XZ21*XZ21) &
           +RI(21)*(YY21*ZZ21+ZZ21*YY21)+RI(22)*YZ21*YZ21
!     (PY PX/PY PY)
            W(46) =                                             &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX22+RI(18)*XX21*YYZZ22  &
           +RI(19)*(YY21*YY22+ZZ21*ZZ22)+RI(20)*(XY21*XY22+XZ21*XZ22) &
           +RI(21)*(YY21*ZZ22+ZZ21*YY22)+RI(22)*YZ21*YZ22
!     (PY PX/PZ S)
            W(47)=X(3)*(RI(8)*XX21+RI(9)*YYZZ21)+RI(10)*(         +Z(3)*XZ21)
!      (PY PX/PZ PX)
            W(48) =                                    &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX31            & 
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ31 &
           +RI(20)*(XY21*XY31+XZ21*XZ31)+RI(22)*YZ21*YZ31
!      (PY PX/PZ PY)
            W(49) =                                    &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX32            & 
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ32 &
           +RI(20)*(XY21*XY32+XZ21*XZ32)+RI(22)*YZ21*YZ32
!      (PY PX/PZ PZ)
            W(50) =                                    &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX33            &
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ33 &
           +RI(20)*XZ21*XZ33 
!     (PY PY/S S)
            W(51)=RI(3)*XX22+RI(4)*YYZZ22
!     (PY PY/PX S)
            W(52)=X(1)*(RI(8)*XX22+RI(9)*YYZZ22)+RI(10)*(Y(1)*XY22+Z(1)*XZ22)
!      (PY PY/PX PX)
            W(53) =                                             &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX11+RI(18)*XX22*YYZZ11  &
           +RI(19)*(YY22*YY11+ZZ22*ZZ11)+RI(20)*(XY22*XY11+XZ22*XZ11) &
           +RI(21)*(YY22*ZZ11+ZZ22*YY11)+RI(22)*YZ22*YZ11
!     (PY PY/PY S)
            W(54)=X(2)*(RI(8)*XX22+RI(9)*YYZZ22)+RI(10)*(Y(2)*XY22+Z(2)*XZ22)
!      (PY PY/PY PX)
            W(55) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX21+RI(18)*XX22*YYZZ21        &
           +RI(19)*(YY22*YY21+ZZ22*ZZ21)+RI(20)*(XY22*XY21+XZ22*XZ21) &
           +RI(21)*(YY22*ZZ21+ZZ22*YY21)+RI(22)*YZ22*YZ21
!      (PY PY/PY PY)
            W(56) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX22+RI(18)*XX22*YYZZ22        &
           +RI(19)*(YY22*YY22+ZZ22*ZZ22)+RI(20)*(XY22*XY22+XZ22*XZ22) &
           +RI(21)*(YY22*ZZ22+ZZ22*YY22)+RI(22)*YZ22*YZ22
!     (PY PY/PZ S)
            W(57)=X(3)*(RI(8)*XX22+RI(9)*YYZZ22)+RI(10)*(         +Z(3)*XZ22)
!      (PY PY/PZ PX)
            W(58) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX31                           &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ31                &
           +RI(20)*(XY22*XY31+XZ22*XZ31)+RI(22)*YZ22*YZ31
!      (PY PY/PZ PY)
            W(59) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX32                           &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ32                &
           +RI(20)*(XY22*XY32+XZ22*XZ32)+RI(22)*YZ22*YZ32
!      (PY PY/PZ PZ)
            W(60) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX33                           &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ33                &
           +RI(20)*XZ22*XZ33
!     (PZ S/SS)
            W(61)=RI(2)*X(3)
!     (PZ S/PX S)
            W(62)=RI(6)*XX31+RI(7)*ZZ31
!     (PZ S/PX PX)
            W(63)=X(3)*(RI(13)*XX11+RI(14)*YYZZ11)+RI(15)*(         +Z(3)*XZ11) 
!     (PZ S/PY S)
            W(64)=RI(6)*XX32+RI(7)*ZZ32
!     (PZ S/PY PX)
            W(65)=X(3)*(RI(13)*XX21+RI(14)*YYZZ21)+RI(15)*(         +Z(3)*XZ21)
!     (PZ S/PY PY)
            W(66)=X(3)*(RI(13)*XX22+RI(14)*YYZZ22)+RI(15)*(         +Z(3)*XZ22) 
!     (PZ S/PZ S)
            W(67)=RI(6)*XX33+RI(7)*ZZ33
!     (PZ S/PZ PX)
            W(68)=X(3)*(RI(13)*XX31+RI(14)*ZZ31)+RI(15)*(         +Z(3)*XZ31)
!     (PZ S/PZ PY)
            W(69)=X(3)*(RI(13)*XX32+RI(14)*ZZ32)+RI(15)*(         +Z(3)*XZ32)
!     (PZ S/PZ PZ)
            W(70)=X(3)*(RI(13)*XX33+RI(14)*ZZ33)+RI(15)*(         +Z(3)*XZ33)
!     (PZ PX/S S)
            W(71)=RI(3)*XX31+RI(4)*ZZ31
!     (PZ PX/PX S)
            W(72)=X(1)*(RI(8)*XX31+RI(9)*ZZ31)+RI(10)*(Y(1)*XY31+Z(1)*XZ31)
!      (PZ PX/PX PX)
            W(73) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX11+RI(18)*XX31*YYZZ11            &
           +RI(19)*ZZ31*ZZ11+RI(20)*(XY31*XY11+XZ31*XZ11)               &
           +RI(21)*ZZ31*YY11+RI(22)*YZ31*YZ11
!     (PZ PX/PY S)                                                              
            W(74)=X(2)*(RI(8)*XX31+RI(9)*ZZ31)+RI(10)*(Y(2)*XY31+Z(2)*XZ31)                                  
!      (PZ PX/PY PX)                                                            
            W(75) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX21+RI(18)*XX31*YYZZ21            &
           +RI(19)*ZZ31*ZZ21+RI(20)*(XY31*XY21+XZ31*XZ21)               &
           +RI(21)*ZZ31*YY21+RI(22)*YZ31*YZ21
!      (PZ PX/PY PY)                                                            
            W(76) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX22+RI(18)*XX31*YYZZ22            &
           +RI(19)*ZZ31*ZZ22+RI(20)*(XY31*XY22+XZ31*XZ22)               &
           +RI(21)*ZZ31*YY22+RI(22)*YZ31*YZ22
!     (PZ PX/PZ S)                                                              
            W(77)=X(3)*(RI(8)*XX31+RI(9)*ZZ31)+RI(10)*(         +Z(3)*XZ31)                                  
!     (PZ PX/PZ PX)                                                     
            W(78) =                                                     &
            (RI(16)*XX31+RI(17)*ZZ31)*XX31                              &        
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ31                              &        
           +RI(20)*(XY31*XY31+XZ31*XZ31)                                &        
           +RI(22)*YZ31*YZ31                                                    
!      (PZ PX/PZ PY)                                                            
            W(79) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX32                               &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ32                              &
           +RI(20)*(XY31*XY32+XZ31*XZ32)                                &
           +RI(22)*YZ31*YZ32
!      (PZ PX/PZ PZ)
            W(80) =                                                     &
            (RI(16)*XX31+RI(17)*ZZ31)*XX33                              &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ33                              &
           +RI(20)*XZ31*XZ33
!     (PZ PY/S S)
            W(81)=RI(3)*XX32+RI(4)*ZZ32
!     (PZ PY/PX S)                                                              
            W(82)=X(1)*(RI(8)*XX32+RI(9)*ZZ32)+RI(10)*(Y(1)*XY32+Z(1)*XZ32)                                 
!      (PZ PY/PX PX)                                                            
            W(83) =                                                      &
           (RI(16)*XX32+RI(17)*ZZ32)*XX11+RI(18)*XX32*YYZZ11             &
           +RI(19)*ZZ32*ZZ11+RI(20)*(XY32*XY11+XZ32*XZ11)                &
           +RI(21)*ZZ32*YY11+RI(22)*YZ32*YZ11                                                    
!     (PZ PY/PY S)                                                              
            W(84)=X(2)*(RI(8)*XX32+RI(9)*ZZ32)+RI(10)*(Y(2)*XY32+Z(2)*XZ32)                                  
!      (PZ PY/PY PX)
            W(85) =                                                      &
           (RI(16)*XX32+RI(17)*ZZ32)*XX21+RI(18)*XX32*YYZZ21             &
           +RI(19)*ZZ32*ZZ21+RI(20)*(XY32*XY21+XZ32*XZ21)                &
           +RI(21)*ZZ32*YY21+RI(22)*YZ32*YZ21
!      (PZ PY/PY PY)
            W(86) =                                                      &
           (RI(16)*XX32+RI(17)*ZZ32)*XX22+RI(18)*XX32*YYZZ22             &
           +RI(19)*ZZ32*ZZ22+RI(20)*(XY32*XY22+XZ32*XZ22)                &
           +RI(21)*ZZ32*YY22+RI(22)*YZ32*YZ22                                                    
!     (PZ PY/PZ S)                                                              
            W(87)=X(3)*(RI(8)*XX32+RI(9)*ZZ32)+RI(10)*(         +Z(3)*XZ32)                                  
!      (PZ PY/PZ PX)                                                            
            W(88) =                                                       &
            (RI(16)*XX32+RI(17)*ZZ32)*XX31+(RI(18)*XX32+RI(19)*ZZ32)*ZZ31 &
           +RI(20)*(XY32*XY31+XZ32*XZ31)+RI(22)*YZ32*YZ31                                                    
!      (PZ PY/PZ PY)
            W(89) =                                                       &
            (RI(16)*XX32+RI(17)*ZZ32)*XX32+(RI(18)*XX32+RI(19)*ZZ32)*ZZ32 &
           +RI(20)*(XY32*XY32+XZ32*XZ32)+RI(22)*YZ32*YZ32                                                    
!       (PZ PY/PZ PZ)
            W(90) =                                                       &
            (RI(16)*XX32+RI(17)*ZZ32)*XX33+(RI(18)*XX32+RI(19)*ZZ32)*ZZ33 &
           +RI(20)*XZ32*XZ33                                                    
!     (PZ PZ/S S)
            W(91)=RI(3)*XX33+RI(4)*ZZ33
!     (PZ PZ/PX S)
            W(92)=X(1)*(RI(8)*XX33+RI(9)*ZZ33)+RI(10)*(          Z(1)*XZ33)
!       (PZ PZ/PX PX)
            W(93) =                                                       &
           (RI(16)*XX33+RI(17)*ZZ33)*XX11+RI(18)*XX33*YYZZ11              &
           +RI(19)*ZZ33*ZZ11+RI(20)*XZ33*XZ11                             &
           +RI(21)*ZZ33*YY11                                                    
!     (PZ PZ/PY S)                                                              
            W(94)=X(2)*(RI(8)*XX33+RI(9)*ZZ33)+RI(10)*(         +Z(2)*XZ33)                                  
!       (PZ PZ/PY PX)                                                           
            W(95) =                                                       &
           (RI(16)*XX33+RI(17)*ZZ33)*XX21+RI(18)*XX33*YYZZ21              &
           +RI(19)*ZZ33*ZZ21+RI(20)*XZ33*XZ21                             &
           +RI(21)*ZZ33*YY21                                                    
!       (PZ PZ/PY PY)
            W(96) =                                                       &
           (RI(16)*XX33+RI(17)*ZZ33)*XX22+RI(18)*XX33*YYZZ22              &
           +RI(19)*ZZ33*ZZ22+RI(20)*XZ33*XZ22+RI(21)*ZZ33*YY22                                                    
!     (PZ PZ/PZ S)
            W(97)=X(3)*(RI(8)*XX33+RI(9)*ZZ33)+RI(10)*(         +Z(3)*XZ33)                                  
!       (PZ PZ/PZ PX)
            W(98) =                                                       &
            (RI(16)*XX33+RI(17)*ZZ33)*XX31+(RI(18)*XX33+RI(19)*ZZ33)*ZZ31 &
            +RI(20)*XZ33*XZ31                                                    
!       (PZ PZ/PZ PY)
            W(99) =                                                       &
            (RI(16)*XX33+RI(17)*ZZ33)*XX32+(RI(18)*XX33+RI(19)*ZZ33)*ZZ32 &
           +RI(20)*XZ33*XZ32                                                    
!       (PZ PZ/PZ PZ)
            W(100) =                                                      &
            (RI(16)*XX33+RI(17)*ZZ33)*XX33+(RI(18)*XX33+RI(19)*ZZ33)*ZZ33 &
           +RI(20)*XZ33*XZ33
            KI = 100                                                         
         ELSE
!     (PX S/S S)
            W(2)=RI(2)*X(1)
!     (PX PX/S S)
            W(3)=RI(3)*XX11+RI(4)*YYZZ11
!     (PY S/S S)
            W(4)=RI(2)*X(2)
!     (PY PX/S S)
            W(5)=RI(3)*XX21+RI(4)*YYZZ21
!     (PY PY/S S)
            W(6)=RI(3)*XX22+RI(4)*YYZZ22
!     (PZ S/SS)
            W(7)=RI(2)*X(3)
!     (PZ PX/S S)
            W(8)=RI(3)*XX31+RI(4)*ZZ31
!     (PZ PY/S S)
            W(9)=RI(3)*XX32+RI(4)*ZZ32
!     (PZ PZ/S S)
            W(10)=RI(3)*XX33+RI(4)*ZZ33
            KI = 10
         END IF
      END IF
! *** NOW ROTATE THE NUCLEAR ATTRACTION INTEGRALS.                              
! *** THE STORAGE OF THE NUCLEAR ATTRACTION INTEGRALS 

      CSS1=qm2_params%core_chg(jqm)*RI(1)
      CSS2=qm2_params%core_chg(iqm)*RI(1) 
      E1B(1)=-CSS1
      E2A(1)=-CSS2
      IF(qm2_params%natomic_orbs(iqm) == 4) THEN
         CSP1 = qm2_params%core_chg(jqm)*RI(2)
         CPPS1 = qm2_params%core_chg(jqm)*RI(3)
         CPPP1 = qm2_params%core_chg(jqm)*RI(4)
         E1B(2) = -CSP1 *X(1)
         E1B(3) = -CPPS1*XX11-CPPP1*YYZZ11
         E1B(4) = -CSP1 *X(2)
         E1B(5) = -CPPS1*XX21-CPPP1*YYZZ21
         E1B(6) = -CPPS1*XX22-CPPP1*YYZZ22
         E1B(7) = -CSP1 *X(3)
         E1B(8) = -CPPS1*XX31-CPPP1*ZZ31
         E1B(9) = -CPPS1*XX32-CPPP1*ZZ32
         E1B(10)= -CPPS1*XX33-CPPP1*ZZ33
      END IF
      IF(qm2_params%natomic_orbs(jqm) == 4) THEN
         CSP2 = qm2_params%core_chg(iqm)*RI(5)
         CPPS2 = qm2_params%core_chg(iqm)*RI(11)
         CPPP2 = qm2_params%core_chg(iqm)*RI(12)
         E2A(2) = -CSP2 *X(1)
         E2A(3) = -CPPS2*XX11-CPPP2*YYZZ11
         E2A(4) = -CSP2 *X(2)
         E2A(5) = -CPPS2*XX21-CPPP2*YYZZ21
         E2A(6) = -CPPS2*XX22-CPPP2*YYZZ22
         E2A(7) = -CSP2 *X(3)
         E2A(8) = -CPPS2*XX31-CPPP2*ZZ31
         E2A(9) = -CPPS2*XX32-CPPP2*ZZ32
         E2A(10)= -CPPS2*XX33-CPPP2*ZZ33
      END IF

!      SCALE = EXP(-ALP(NI)*RIJ)+EXP(-ALP(NJ)*RIJ)
      SCALE = EXP1i+EXP1j

      if(ni == 1 .AND. (nj == 7 .OR. nj == 8)) then
         SCALE = SCALE+(RIJ-1.0D0)*EXP1j
      elseif ((ni == 7 .OR. ni == 8) .AND. nj == 1) then
         SCALE = SCALE+(RIJ-1.0D0)*EXP1i
      end if
      C1 = qm2_params%core_chg(IQM)*qm2_params%core_chg(JQM)      
      ENUC = C1*RI(1)
      SCALE=ABS(SCALE*ENUC)
      IF(qmmm_struct%AM1_OR_PM3) THEN
!        Add gaussians.
         anam1=0.0d0
         do I=1,qm2_params%num_fn(qmitype)                                                      
            temp_real=RIJ-qm2_params%FN3(i,qmitype)
            temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
            if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
              anam1=anam1+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)
            end if
         end do
         do i=1,qm2_params%num_fn(qmjtype)
            temp_real=RIJ-qm2_params%FN3(i,qmjtype)
            temp_real2=qm2_params%FN2(i,qmjtype)*temp_real*temp_real
            if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
              anam1=anam1+qm2_params%FN1(i,qmjtype)*EXP(-temp_real2)
            end if
         end do
         anam1=anam1*c1*oneRIJ
         scale = scale + anam1 
      ENDIF

!PDDG Specific terms
      if (qmmm_struct%PDDG_IN_USE) then
         PDDG_CORR = qm2_params%PDDG_TERM1(qmitype,qmjtype)*PDDG_EXP1 + &
                     qm2_params%PDDG_TERM2(qmitype,qmjtype)*PDDG_EXP2 + &
                     qm2_params%PDDG_TERM3(qmitype,qmjtype)*PDDG_EXP3 + &
                     qm2_params%PDDG_TERM4(qmitype,qmjtype)*PDDG_EXP4
         SCALE = SCALE + PDDG_CORR
      end if

      ENUC=ENUC+SCALE
      KR=KR+KI
      RETURN
end subroutine qm2_rotate_qmqm

