#include "dprec.h"
!!3.10.2010
!!HB.h modified by Antonios for amber10
!! input GO hydrogen bonds
!! store atom number of theta,phi angles 
!! H-B between JTHEI-JPHII, JTHEI<JPHII

_REAL_ amp,EHBA,EHBV,ECHI,XKCHI,THE0,PHI0,U,ENBT
integer nhb_pair,IXI,IXI1,IXI2,IXI3,num,ixj,ixk,ixl,N_HB,maphb
integer NXI,NXI1,NXI2,NXI3,counter,dummyi,dummyj
_REAL_ cosphi,cosphip,acosphi
_REAL_ XIJ,YIJ,ZIJ,XKJ,YKJ,ZKJ
_REAL_ XLK,YLK,ZLK,TX,TY,TZ,UX,UY,UZ
_REAL_ DTXX, DTY,DTZ,DUX,DUY,DUZ
_REAL_ RT,RU,dummy,tempA,tempB,tempC,Chi,dChi,cosphia 
_REAL_ RT2,RTU,RU2,FXI,FYI,FZI,FXJ,FYJ,FZJ,FXK,FYK,FZK 
_REAL_ FXl,FYl,FZl 
_REAL_ FFxi,FFyi,FFzi,FFxj,FFyj,FFzj,FFxk,FFyk,FFzk,FFxl,FFyl,FFzl 
_REAL_ delx_ixi_ixj,dely_ixi_ixj,delz_ixi_ixj,delx_ixk_ixl,dely_ixk_ixl,delz_ixk_ixl
_REAL_ delx_ixj_ixk,dely_ixj_ixk,delz_ixj_ixk

data cosphia /0.893373d0/

common/HYDROB/NXI(92000),NXI1(92000),NXI2(92000),NXI3(92000)
!! store the native Go phi,psi angles 
common/HBANG/THE0(10000),PHI0(10000)
!!     check contact HB map
common/MAP/maphb(10000,10000)
!!     output ener 
common/ENERTW/EHBA,EHBV,ECHI,XKCHI
!!     input Amp 
common/DISPLACE/amp,nhb_pair
!!AS   delx,dely,delz arrays for the HB many body interactions
common/DELMAP1/delx_ixi_ixj(92000),dely_ixi_ixj(92000),delz_ixi_ixj(92000)
common/DELMAP2/delx_ixk_ixl(92000),dely_ixk_ixl(92000),delz_ixk_ixl(92000)
common/DELMAP3/delx_ixj_ixk(92000),dely_ixj_ixk(92000),delz_ixj_ixk(92000)
