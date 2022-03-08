! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
subroutine qm2_repp(IQM, JQM,R,RR2,RI,SQRTAEE)
!***********************************************************************
!  OPTIMISATION BY ROSS WALKER (TSRI,2005)
!
!  REPP CALCULATES THE TWO-ELECTRON REPULSION INTEGRALS
!
!     ON INPUT R,RR2     = INTERATOMIC DISTANCE, ^2
!              IQM, JQM = number of the QM pair we are dealing with from
!                         a 1 to nquant loop.
!
!    ON OUTPUT RI      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS
!
! *** THIS ROUTINE COMPUTES THE TWO-CENTRE REPULSION INTEGRALS
!
!***********************************************************************
      use qmmm_module, only : qm2_params
      use constants, only :  AU_TO_EV, HALF_AU_TO_EV, FOURTH_AU_TO_EV, EIGHTH_AU_TO_EV, SXNTH_AU_TO_EV

      implicit none

! Passed in
      integer, intent(in) :: iqm, jqm
      _REAL_, intent(in) :: R, RR2, sqrtaee
      _REAL_, intent(out) :: RI(22)

!Local
      LOGICAL SI,SJ
      _REAL_ ARG(71),SQR(71)
      _REAL_ DA,QA,ADE,AQE,XXX, DB, QB, AED, AEQ, AXX, ADQ, AQD, AQQ
      _REAL_ YYY,ZZZ,WWW,DZE,QZZE,QXXE,EDZ,EQZZ,EQXX,DXDX,DZDZ,DZQXX,QXXDZ
      _REAL_ DZQZZ,QZZDZ,QXXQXX,QXXQYY,QXXQZZ,QZZQXX,QZZQZZ,DXQXZ,QXZDX,QXZQXZ

!    ATOMIC UNITS ARE USED IN THE CALCULATION,
!    FINAL RESULTS ARE CONVERTED TO EV

      SI = (qm2_params%natomic_orbs(iqm) > 1)
      SJ = (qm2_params%natomic_orbs(jqm) > 1)

!ALL QM have S orbital

!     (SS/SS)
      RI(1) = AU_TO_EV*SQRTAEE

      IF (SI .AND. (.NOT.SJ)) THEN
!     HEAVY ATOM - HYDROGEN                                                     
         DA=qm2_params%multip_2c_elec_params(1,iqm)
         QA=qm2_params%multip_2c_elec_params(2,iqm) * 2.0D0
         ADE = qm2_params%multip_2c_elec_params(4,iqm) + qm2_params%multip_2c_elec_params(3,jqm)
         ADE = ADE * ADE
         AQE = qm2_params%multip_2c_elec_params(5,iqm) + qm2_params%multip_2c_elec_params(3,jqm)
         AQE = AQE * AQE
         XXX = R+DA
         ARG(1) = XXX*XXX + ADE
         XXX = R-DA
         ARG(2) = XXX*XXX + ADE
         XXX = R+QA
         ARG(3) = XXX*XXX + AQE
         XXX = R-QA
         ARG(4) = XXX*XXX + AQE
         ARG(5) = rr2 + AQE
         ARG(6) = ARG(5) + QA*QA
         ! do I = 1,6
            !Ross Walker - Inverted this for speed
            ! SQR(I) = 1.0D0/SQRT(ARG(I))
         ! end do
         call vdinvsqrt( 6, arg, sqr )
         RI(2) = HALF_AU_TO_EV*SQR(1) - HALF_AU_TO_EV*SQR(2)
         RI(3) = RI(1) + FOURTH_AU_TO_EV*SQR(3) + FOURTH_AU_TO_EV*SQR(4) - HALF_AU_TO_EV*SQR(5)
         RI(4) = RI(1) + HALF_AU_TO_EV*SQR(6) - HALF_AU_TO_EV*SQR(5)
      ELSE IF ((.NOT.SI).AND.SJ) THEN                                           
!     HYDROGEN - HEAVY ATOM                                                     
         DB=qm2_params%multip_2c_elec_params(1,jqm)
         QB=qm2_params%multip_2c_elec_params(2,jqm) * 2.0D0
         AED = qm2_params%multip_2c_elec_params(3,iqm) + qm2_params%multip_2c_elec_params(4,jqm)
         AED = AED * AED
         AEQ = qm2_params%multip_2c_elec_params(3,iqm) + qm2_params%multip_2c_elec_params(5,jqm)
         AEQ = AEQ * AEQ
         XXX = R-DB
         ARG(1) = XXX*XXX + AED
         XXX = R+DB
         ARG(2) = XXX*XXX + AED
         XXX = R-QB
         ARG(3) = XXX*XXX + AEQ 
         XXX = R+QB
         ARG(4) = XXX*XXX + AEQ
         ARG(5) = rr2 + AEQ
         ARG(6) = ARG(5) + QB*QB
         ! do I = 1,6
            !Ross Walker - Inverted this for speed
            ! SQR(I) = 1.0D0/SQRT(ARG(I))
         ! end do
         call vdinvsqrt( 6, arg, sqr )
         RI(5) = HALF_AU_TO_EV*SQR(1) - HALF_AU_TO_EV*SQR(2)
         RI(11) = RI(1) + FOURTH_AU_TO_EV*SQR(3) + FOURTH_AU_TO_EV*SQR(4) &
                  - HALF_AU_TO_EV*SQR(5)
         RI(12) = RI(1) + HALF_AU_TO_EV*SQR(6) - HALF_AU_TO_EV*SQR(5)
      ELSE IF (SI .AND. SJ) then
!     HEAVY ATOM - HEAVY ATOM                                                   
!     DEFINE CHARGE SEPARATIONS.                                                
         DA=qm2_params%multip_2c_elec_params(1,iqm)
         DB=qm2_params%multip_2c_elec_params(1,jqm)
         QA=qm2_params%multip_2c_elec_params(2,iqm) * 2.0D0
         QB=qm2_params%multip_2c_elec_params(2,jqm) * 2.0D0
         ADE = qm2_params%multip_2c_elec_params(4,iqm) + qm2_params%multip_2c_elec_params(3,jqm)
         ADE = ADE * ADE
         AQE = qm2_params%multip_2c_elec_params(5,iqm) + qm2_params%multip_2c_elec_params(3,jqm)
         AQE = AQE * AQE
         AED = qm2_params%multip_2c_elec_params(3,iqm) + qm2_params%multip_2c_elec_params(4,jqm)
         AED = AED * AED
         AEQ = qm2_params%multip_2c_elec_params(3,iqm) + qm2_params%multip_2c_elec_params(5,jqm)
         AEQ = AEQ * AEQ
         AXX = qm2_params%multip_2c_elec_params(4,iqm) + qm2_params%multip_2c_elec_params(4,jqm)
         AXX = AXX * AXX
         ADQ = qm2_params%multip_2c_elec_params(4,iqm) + qm2_params%multip_2c_elec_params(5,jqm)
         ADQ = ADQ * ADQ
         AQD = qm2_params%multip_2c_elec_params(5,iqm) + qm2_params%multip_2c_elec_params(4,jqm)
         AQD = AQD * AQD
         AQQ = qm2_params%multip_2c_elec_params(5,iqm) + qm2_params%multip_2c_elec_params(5,jqm)
         AQQ = AQQ * AQQ
         XXX = R + DA
         ARG(1) = XXX * XXX + ADE
         XXX = R - DA
         ARG(2) = XXX*XXX + ADE
         XXX = R - QA
         ARG(3) = XXX*XXX + AQE
         XXX = R + QA
         ARG(4) = XXX*XXX + AQE
         ARG(5) = rr2 + AQE
         ARG(6) = ARG(5) + QA*QA
         XXX = R-DB
         ARG(7) = XXX*XXX + AED
         XXX = R+DB
         ARG(8) = XXX*XXX + AED
         XXX = R - QB
         ARG(9) = XXX*XXX + AEQ
         XXX = R + QB
         ARG(10) = XXX*XXX + AEQ
         ARG(11) = rr2 + AEQ
         ARG(12) = ARG(11) + QB*QB
         XXX = DA-DB
         ARG(13) = rr2 + AXX + XXX*XXX
         XXX = DA+DB
         ARG(14) = rr2 + AXX + XXX*XXX
         XXX = R + DA - DB
         ARG(15) = XXX*XXX + AXX
         XXX = R - DA + DB
         ARG(16) = XXX*XXX + AXX
         XXX = R - DA - DB
         ARG(17) = XXX*XXX + AXX
         XXX = R + DA + DB
         ARG(18) = XXX*XXX + AXX
         XXX = R + DA
         ARG(19) = XXX*XXX + ADQ
         ARG(20) = ARG(19) + QB*QB
         XXX = R - DA
         ARG(21) = XXX*XXX + ADQ
         ARG(22) = ARG(21) + QB*QB
         XXX = R - DB
         ARG(23) = XXX*XXX + AQD
         ARG(24) = ARG(23) + QA*QA
         XXX = R + DB
         ARG(25) = XXX*XXX + AQD
         ARG(26) = ARG(25) + QA*QA
         XXX = R + DA - QB
         ARG(27) = XXX*XXX + ADQ
         XXX = R - DA - QB
         ARG(28) = XXX*XXX + ADQ
         XXX = R + DA + QB
         ARG(29) = XXX*XXX + ADQ
         XXX = R - DA + QB
         ARG(30) = XXX*XXX + ADQ
         XXX = R + QA - DB
         ARG(31) = XXX*XXX + AQD
         XXX = R + QA + DB
         ARG(32) = XXX*XXX + AQD
         XXX = R - QA - DB
         ARG(33) = XXX*XXX + AQD
         XXX = R - QA + DB
         ARG(34) = XXX*XXX + AQD
         ARG(35) = rr2 + AQQ
         XXX = QA - QB
         ARG(36) = ARG(35) + XXX*XXX
         XXX = QA + QB
         ARG(37) = ARG(35) + XXX*XXX
         ARG(38) = ARG(35) + QA*QA
         ARG(39) = ARG(35) + QB*QB
         ARG(40) = ARG(38) + QB*QB
         XXX = R - QB
         ARG(41) = XXX*XXX + AQQ
         ARG(42) = ARG(41) + QA*QA
         XXX = R + QB
         ARG(43) = XXX*XXX + AQQ
         ARG(44) = ARG(43) + QA*QA
         XXX = R + QA
         ARG(45) = XXX*XXX + AQQ
         ARG(46) = ARG(45) + QB*QB
         XXX = R - QA
         ARG(47) = XXX*XXX + AQQ
         ARG(48) = ARG(47) + QB*QB
         XXX = R + QA - QB
         ARG(49) = XXX*XXX + AQQ
         XXX = R + QA + QB
         ARG(50) = XXX*XXX + AQQ
         XXX = R - QA - QB
         ARG(51) = XXX*XXX + AQQ
         XXX = R - QA + QB
         ARG(52) = XXX*XXX + AQQ
         QA=qm2_params%multip_2c_elec_params(2,iqm)
         QB=qm2_params%multip_2c_elec_params(2,jqm)
         XXX = DA - QB
         XXX = XXX*XXX
         YYY = R - QB
         YYY = YYY*YYY
         ZZZ = DA + QB
         ZZZ = ZZZ*ZZZ
         WWW = R + QB
         WWW = WWW*WWW
         ARG(53) = XXX + YYY + ADQ
         ARG(54) = XXX + WWW + ADQ
         ARG(55) = ZZZ + YYY + ADQ
         ARG(56) = ZZZ + WWW + ADQ
         XXX = QA - DB
         XXX = XXX*XXX
         YYY = QA + DB
         YYY = YYY*YYY
         ZZZ = R + QA
         ZZZ = ZZZ*ZZZ
         WWW = R - QA
         WWW = WWW*WWW
         ARG(57) = ZZZ + XXX + AQD
         ARG(58) = WWW + XXX + AQD
         ARG(59) = ZZZ + YYY + AQD
         ARG(60) = WWW + YYY + AQD
         XXX = QA - QB
         XXX = XXX*XXX
         ARG(61) = ARG(35) + 2.0D0*XXX
         YYY = QA + QB
         YYY = YYY*YYY
         ARG(62) = ARG(35) + 2.0D0*YYY
         ARG(63) = ARG(35) + 2.0D0*(QA*QA+QB*QB)
         ZZZ = R + QA - QB
         ZZZ = ZZZ*ZZZ
         ARG(64) = ZZZ + XXX + AQQ
         ARG(65) = ZZZ + YYY + AQQ
         ZZZ = R + QA + QB
         ZZZ = ZZZ*ZZZ
         ARG(66) = ZZZ + XXX + AQQ
         ARG(67) = ZZZ + YYY + AQQ
         ZZZ = R - QA - QB
         ZZZ = ZZZ*ZZZ
         ARG(68) = ZZZ + XXX + AQQ
         ARG(69) = ZZZ + YYY + AQQ
         ZZZ = R - QA + QB
         ZZZ = ZZZ*ZZZ
         ARG(70) = ZZZ + XXX + AQQ
         ARG(71) = ZZZ + YYY + AQQ
         ! do I = 1,71
            !Ross Walker - Inverted this for speed
            ! SQR(I) = 1.0d0/SQRT(ARG(I))
         ! end do
         call vdinvsqrt( 71, arg, sqr )
         DZE = -HALF_AU_TO_EV*SQR(1) + HALF_AU_TO_EV*SQR(2)
         QZZE = FOURTH_AU_TO_EV*SQR(3) + FOURTH_AU_TO_EV*SQR(4) - HALF_AU_TO_EV*SQR(5)
         QXXE = HALF_AU_TO_EV*SQR(6) - HALF_AU_TO_EV*SQR(5)
         EDZ = - HALF_AU_TO_EV*SQR(7) + HALF_AU_TO_EV*SQR(8)
         EQZZ  = FOURTH_AU_TO_EV*SQR(9) + FOURTH_AU_TO_EV*SQR(10) - HALF_AU_TO_EV*SQR(11)
         EQXX  = HALF_AU_TO_EV*SQR(12) - HALF_AU_TO_EV*SQR(11)
         DXDX  = HALF_AU_TO_EV*SQR(13) - HALF_AU_TO_EV*SQR(14)
         DZDZ  = FOURTH_AU_TO_EV*SQR(15) + FOURTH_AU_TO_EV*SQR(16) &
                 - FOURTH_AU_TO_EV*SQR(17) - FOURTH_AU_TO_EV*SQR(18)
         DZQXX =  FOURTH_AU_TO_EV*SQR(19) - FOURTH_AU_TO_EV*SQR(20) &
                  - FOURTH_AU_TO_EV*SQR(21) + FOURTH_AU_TO_EV*SQR(22)
         QXXDZ =  FOURTH_AU_TO_EV*SQR(23) - FOURTH_AU_TO_EV*SQR(24) &
                  - FOURTH_AU_TO_EV*SQR(25) + FOURTH_AU_TO_EV*SQR(26)
         DZQZZ = -EIGHTH_AU_TO_EV*SQR(27) + EIGHTH_AU_TO_EV*SQR(28) &
                 - EIGHTH_AU_TO_EV*SQR(29) + EIGHTH_AU_TO_EV*SQR(30) &
                 - FOURTH_AU_TO_EV*SQR(21) + FOURTH_AU_TO_EV*SQR(19)
         QZZDZ = -EIGHTH_AU_TO_EV*SQR(31) + EIGHTH_AU_TO_EV*SQR(32) &
                 - EIGHTH_AU_TO_EV*SQR(33) + EIGHTH_AU_TO_EV*SQR(34) &
                 + FOURTH_AU_TO_EV*SQR(23) - FOURTH_AU_TO_EV*SQR(25)
         QXXQXX = EIGHTH_AU_TO_EV*SQR(36) + EIGHTH_AU_TO_EV*SQR(37) &
                  - FOURTH_AU_TO_EV*SQR(38) - FOURTH_AU_TO_EV*SQR(39) &
                  + FOURTH_AU_TO_EV*SQR(35)
         QXXQYY = FOURTH_AU_TO_EV*SQR(40) - FOURTH_AU_TO_EV*SQR(38) &
                  - FOURTH_AU_TO_EV*SQR(39) + FOURTH_AU_TO_EV*SQR(35)
         QXXQZZ = EIGHTH_AU_TO_EV*SQR(42) + EIGHTH_AU_TO_EV*SQR(44) &
                  - EIGHTH_AU_TO_EV*SQR(41) - EIGHTH_AU_TO_EV*SQR(43) &
                  - FOURTH_AU_TO_EV*SQR(38) + FOURTH_AU_TO_EV*SQR(35)
         QZZQXX = EIGHTH_AU_TO_EV*SQR(46) + EIGHTH_AU_TO_EV*SQR(48) &
                  - EIGHTH_AU_TO_EV*SQR(45) - EIGHTH_AU_TO_EV*SQR(47) &
                  - FOURTH_AU_TO_EV*SQR(39) + FOURTH_AU_TO_EV*SQR(35)
         QZZQZZ = SXNTH_AU_TO_EV*SQR(49) + SXNTH_AU_TO_EV*SQR(50)  &
                  + SXNTH_AU_TO_EV*SQR(51) + SXNTH_AU_TO_EV*SQR(52) &
                  - EIGHTH_AU_TO_EV*SQR(47) - EIGHTH_AU_TO_EV*SQR(45) &
                  - EIGHTH_AU_TO_EV*SQR(41) - EIGHTH_AU_TO_EV*SQR(43) &
                  + FOURTH_AU_TO_EV*SQR(35)
         DXQXZ = -FOURTH_AU_TO_EV*SQR(53) + FOURTH_AU_TO_EV*SQR(54) &
                 + FOURTH_AU_TO_EV*SQR(55) - FOURTH_AU_TO_EV*SQR(56)
         QXZDX = -FOURTH_AU_TO_EV*SQR(57) + FOURTH_AU_TO_EV*SQR(58) &
                 + FOURTH_AU_TO_EV*SQR(59) - FOURTH_AU_TO_EV*SQR(60)
         QXZQXZ = EIGHTH_AU_TO_EV*SQR(64) - EIGHTH_AU_TO_EV*SQR(66) &
                  - EIGHTH_AU_TO_EV*SQR(68) + EIGHTH_AU_TO_EV*SQR(70) &
                  - EIGHTH_AU_TO_EV*SQR(65) + EIGHTH_AU_TO_EV*SQR(67) &
                  + EIGHTH_AU_TO_EV*SQR(69) - EIGHTH_AU_TO_EV*SQR(71)
         RI(2) = -DZE
         RI(3) = RI(1) + QZZE
         RI(4) = RI(1) + QXXE
         RI(5) = -EDZ
         RI(6) = DZDZ
         RI(7) = DXDX
         RI(8) = -EDZ -QZZDZ
         RI(9) = -EDZ -QXXDZ
         RI(10) = -QXZDX
         RI(11) =  RI(1) + EQZZ
         RI(12) =  RI(1) + EQXX
         RI(13) = -DZE -DZQZZ
         RI(14) = -DZE -DZQXX
         RI(15) = -DXQXZ
         RI(16) = RI(1) +EQZZ +QZZE +QZZQZZ
         RI(17) = RI(1) +EQZZ +QXXE +QXXQZZ
         RI(18) = RI(1) +EQXX +QZZE +QZZQXX
         RI(19) = RI(1) +EQXX +QXXE +QXXQXX
         RI(20) = QXZQXZ
         RI(21) = RI(1) +EQXX +QXXE +QXXQYY
         RI(22) = 0.5D0 * (QXXQXX -QXXQYY)
      END IF
      RETURN
end subroutine qm2_repp
