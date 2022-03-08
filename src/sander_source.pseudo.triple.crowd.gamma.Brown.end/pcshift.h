parameter (mshfd=12000, maxfe=5, &
      mpar=5, mp=maxfe*mpar+1, np=maxfe*mpar)
common /simplex/ &
      coox(mshfd),cooy(mshfd),cooz(mshfd), &
      cmx(maxfe),cmy(maxfe),cmz(maxfe), &
      tolpro(mshfd),wt(mshfd),shift(mshfd), &
      obs(mshfd),mltpro(mshfd),iprot(mshfd), &
      obsp(mshfd),shavp(mshfd), &
      phi(maxfe),teta(maxfe),omg(maxfe), &
      a1dip(maxfe),a2dip(maxfe), &
      optphi(maxfe),opttet(maxfe),optomg(maxfe), &
      opta1(maxfe),opta2(maxfe), &
      toldip,resid,oldresid,optkon, &
      ippmc(maxfe),ioldvio,iviolation, &
      nhp,nfe,nprot,ipear,nstampa,nstampaten
common/derive/ d(3,mshfd)

