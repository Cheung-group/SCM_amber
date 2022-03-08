#include "dprec.h"
!-------------BEGIN    md.h  ------------------------------------------------
integer BC_MDI  ! size in integers of common block mdi
integer BC_MDR  ! size in Reals of common block mdr

! ... integer variables:

integer nrp,nspm,ig,ntx,ntcx,           &!5
      ntxo,ntt,ntp,ntr,init,             &!10
      ntcm,nscm,isolvp,nsolut,klambda,   &!15
      ntc,ntcc,ntf,ntid,ntn,             &!20
      ntnb,nsnb,ndfmin,nstlim,nrc,       &!25
      ntrx,npscal,imin,maxcyc,ncyc,      &!30
      ntmin,irest,jfastw,                &!33
      ibgwat,ienwat,iorwat,              &!36
      iwatpr,nsolw,igb,alpb,iyammp,           &!41
      gbsa,vrand,iwrap,nrespa,irespa,nrespai,icfe,  &!48
      rbornstat,ivcap,iconstreff,        &!51
      neb,vv,tmode,ipol,iesp,ievb,nodeid,num_noshake,    &!59
      idecomp,icnstph,ntcnstph,maxdup,numexchg,repcrd,numwatkeep,hybridgb,  &!67
      ibgion,ienion,profile_mpi                      !70
parameter (BC_MDI=70)

common/mdi/nrp,nspm,ig, &
      ntx,ntcx,ntxo,ntt,ntp,ntr,init,ntcm,nscm, &
      isolvp,nsolut,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,ndfmin, &
      nstlim,nrc,ntrx,npscal,imin,maxcyc,ncyc,ntmin, &
      irest,jfastw,ibgwat,ienwat,iorwat, &
      iwatpr,nsolw,igb,alpb,iyammp,gbsa,vrand,numexchg,repcrd,numwatkeep,hybridgb, &
      iwrap,nrespa,irespa,nrespai,icfe,rbornstat, &
      ivcap,iconstreff,idecomp,klambda,icnstph,ntcnstph,maxdup,neb,vv, &
  tmode,ipol,iesp,ievb,nodeid,num_noshake,ibgion,ienion, profile_mpi

! ... floats:

_REAL_ t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, & !9
      tol,taur,dx0,drms,vlimit,rbtarg(9),tmass,tmassinv,  & !25
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt,  & !32
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon,  & !38
      solvph,rgbmax,fsmax,restraint_wt, &  !42
      skmin,skmax,vfac,gbneckscale,v11,v12,v22,kevb,evbt,Arad   !52
parameter (BC_MDR=52)
common/mdr/t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, &
      tol,taur,dx0,drms,vlimit,rbtarg,tmass,tmassinv, &
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt, &
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon, &
      solvph,rgbmax,fsmax,restraint_wt,skmin,skmax,vfac,gbneckscale, &
      v11,v12,v22,kevb,evbt,Arad

! ... strings:

character(len=4) iwtnm,iowtnm,ihwtnm
  character(len=256) restraintmask,bellymask,tgtfitmask,tgtrmsmask,noshakemask,crgmask
common/mds/ restraintmask,bellymask,tgtfitmask,tgtrmsmask,noshakemask,crgmask,  &
            iwtnm,iowtnm,ihwtnm(2)

!-------------END    md.h  ------------------------------------------------

