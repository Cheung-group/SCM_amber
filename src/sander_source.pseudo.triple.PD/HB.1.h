#include "dprec.h"
!!HB.1.h modified by Antonios for amber10
!! margaret added. 9.18.01 Hydrogen Bond. real.
!! input GO hydrogen bonds
!! store atom number of theta,phi angles 
!! H-B between JTHEI-JPHII, JTHEI<JPHII
!! To be included in ene.f in amber 10. If HB.h is used there are conflicts for some variable names




_REAL_ EHBA,EHBV,maphb

common/MAP/maphb(10000,10000)
common/ENERTW/EHBA,EHBV
