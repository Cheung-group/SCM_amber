#include "dprec.h"
!!9.2.2010
!!CHI.h modified by Antonios for amber10
!!Chiral energy term header file

_REAL_ XCHI,xxkchi,xxchi,MAXLIM,temp_chi
integer ICHI,ICA,ICB,INC,ICC,nbeta
integer iica,iicb,iinc,iicc,i_chi
_REAL_ XAX,XAY,XAZ,XBX,XBY,XBZ,XCX,XCY,XCZ,xtriple,echii
_REAL_ delx_xax,dely_xay,delz_xaz,delx_xbx,dely_xby,delz_xbz,delx_xcx,dely_xcy,delz_xcz
! Qian change
!data MAXLIM /1.0d4/
data MAXLIM /1.0d2/
! Qian change end
_REAL_ Fxi1_chi,Fxi2_chi,Fxi3_chi,Fxi_chi
_REAL_ Fyi1_chi,Fyi2_chi,Fyi3_chi,Fyi_chi
_REAL_ Fzi1_chi,Fzi2_chi,Fzi3_chi,Fzi_chi


COMMON/XCHIR/XCHI(10000)
COMMON/NCHI/ICHI
COMMON/IBETA/ICA(10000),ICB(10000),INC(10000),ICC(10000)
!!AS   delx,dely,delz arrays for the Chiral interaction
common/DELMAPCHI1/delx_xax(10000),dely_xay(10000),delz_xaz(10000)
common/DELMAPCHI2/delx_xbx(10000),dely_xby(10000),delz_xbz(10000)
common/DELMAPCHI3/delx_xcx(10000),dely_xcy(10000),delz_xcz(10000)

