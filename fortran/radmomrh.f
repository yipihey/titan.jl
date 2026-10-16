      subroutine radmomrh
c
c***********************************************************************
c     radiation momentum equation
c.......................................................................
c     calling sequence: matgen > radmomrh > advecti
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq2'
c
c-----------------------------------------------------------------------
c     signum (x) 
c-----------------------------------------------------------------------
c
      sgn (x) = cvmgp ( +1.D0, -1.D0, x )
c
c***********************************************************************
c     radiation diffusion; equation rd1 - rd4
c***********************************************************************
c
c++++++++ ltran .gt. 1 +++++++++++++++++++++++++++++++++++++++++++ begin
      if (ltran .gt. 1) then
c
      xlam = (2.0D0 / 3.0D0) * float(lam)
c
      do 1 k = ngrs, ngre
c
c     equations rd6 - rd13
c
      chdvol   = 0.5D0 *  avchin(k) *  ( dn(k-1)  * dvoln   (k-1)
     .                                 + dn(k  )  * dvoln   (k  ))
c
      chdvolrp = 0.5D0 *  avchin(k) *  ( dn(k  )  * rmup1n  (k+1))
c
      chdvolr0 = 0.5D0 *  avchin(k) *  ( dn(k-1)  * rmup1n  (k  ) 
     .                                 - dn(k  )  * rmup1n  (k  ))
c
      chdvolrm = 0.5D0 *  avchin(k) *  (-dn(k-1)  * rmup1n  (k-1))
c
      chdvold0 = 0.5D0 *  avchin(k) *  ( dn(k  )  * dvoln   (k  ))
     .         + 0.5D0 * (avchin(k) / chifn(k  )) * dlcfdldn(k  )
     .                                            * chdvol
c
      chdvoldm = 0.5D0 *  avchin(k) *   (dn(k-1)  * dvoln   (k-1))
     .         + 0.5D0 * (avchin(k) / chifn(k-1)) * dlcfdldn(k-1)
     .                                            * chdvol
c
      chdvolt0 = 0.5D0 * (avchin(k) / chifn(k  )) * dlcfdltn(k  )
     .                                            * chdvol
c
      chdvoltm = 0.5D0 * (avchin(k) / chifn(k-1)) * dlcfdltn(k-1) 
     .                                            * chdvol
c
c     equations rd14 - rd17
c
      elim   =   xlam * rmun(k) * abs( ern(k-1) - ern(k) )
     .                                        /    ( ern(k-1) + ern(k) )
c
      elimr0 =   xmu  * elim
c
      elime0 = 
     .     - (elim * ern(k  ) + xlam * rmun(k) * sgn(ern(k-1) - ern(k)))
     .                                         /    (ern(k-1) + ern(k))
c
      elimem = 
     .     - (elim * ern(k-1) - xlam * rmun(k) * sgn(ern(k-1) - ern(k)))
     .                                         /    (ern(k-1) + ern(k))
c
c     equations rd18 - rd21
c
       fd    =   cc * rmun(k) * (ern(k-1) - ern(k)) / 3.0D0
      dfdde0 = - cc * rmun(k) *             ern(k)  / 3.0D0
      dfddem =   cc * rmun(k) *  ern(k-1)           / 3.0D0
      dfddr0 =  xmu * fd
c
c     equation rd22
c
      denom  (k) =  1.0D0 / ( chdvol + elim + 1.0D-30 )
c
c     equations rd23 - rd32
c
       dif   (k) =   denom(k) * fd
      ddifdrp(k) =                       - dif(k) * denom(k) * chdvolrp
      ddifdr0(k) =   denom(k) * dfddr0   - dif(k) * denom(k) * chdvolr0
     .                                   - dif(k) * denom(k) * elimr0
      ddifdrm(k) =                       - dif(k) * denom(k) * chdvolrm
      ddifdd0(k) =                       - dif(k) * denom(k) * chdvold0
      ddifddm(k) =                       - dif(k) * denom(k) * chdvoldm
      ddifdt0(k) =                       - dif(k) * denom(k) * chdvolt0
      ddifdtm(k) =                       - dif(k) * denom(k) * chdvoltm
      ddifde0(k) =   denom(k) * dfdde0   - dif(k) * denom(k) * elime0
      ddifdem(k) =   denom(k) * dfddem   - dif(k) * denom(k) * elimem
    1 continue
c
      do 2 k = ngrs + 1, ngre
c
c     equations rd33 - rd42
c
      ep1(if, jr, k) = ep1(if, jr, k)              - ddifdrp(k)
      e00(if, jr, k) = e00(if, jr, k)              - ddifdr0(k)
      em1(if, jr, k) = em1(if, jr, k)              - ddifdrm(k)
c
      e00(if, jd, k) = e00(if, jd, k)              - ddifdd0(k)
      em1(if, jd, k) = em1(if, jd, k)              - ddifddm(k)
c
      e00(if, jt, k) = e00(if, jt, k)              - ddifdt0(k)
      em1(if, jt, k) = em1(if, jt, k)              - ddifdtm(k)
c
      e00(if, je, k) = e00(if, je, k)              - ddifde0(k)
      em1(if, je, k) = em1(if, je, k)              - ddifdem(k)
c
      e00(if, jf, k) = e00(if, jf, k) + frnom(k)
c
c     equation rd43
c
      rhs(if, k) =                      frn  (k)   - dif(k)
    2 continue
c     
      end if
c
c++++++++ ltran .gt. 1 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c***********************************************************************
c     radiation transport; equation rm2
c***********************************************************************
c
c++++++++ ltran .le. 1 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (ltran .le. 1) then
c
c=======================================================================
c     set cell-center fluxes for advection
c=======================================================================
c
      do 4 l = 1, madv
      do 3 k = 1, mgr
      adv(k, l) = 0.0D0
    3 continue
    4 continue
c
      do 5   k  = ngrs - 1, ngre + 2
      qn    (k) = frn (k) * 2.0D0 / ( dn(k) + dn(k-1) )
      qo    (k) = fro (k) * 2.0D0 / ( do(k) + do(k-1) )
    5 continue
      do 55  k  = ngrs, ngre + 1
      qso   (k) = frso(k)
   55 continue
      do 555 k  = ngrs, ngre
      flow  (k) = - ( dmdt(k) + dmdt(k+1) )
  555 continue
c
      call advecti
c
      do 6  k = ngrs, ngre
      frb  (k) = qb (k)
    6 continue
      do 66 k  = ngrs, ngre + 1
      frsn (k) = qsn(k)
   66 continue
c
c=======================================================================
c     set auxiliary derivatives
c=======================================================================
c
c     equations rm12 - rm15
c
      do 7 k = ngrs, ngre
      dqbdlfp2(k) = + dqbdqp2(k) * frnom(k+2) 
     .                           / ( dn(k+2) + dn(k+1) ) * 2.D0 
      dqbdlfp1(k) = + dqbdqp1(k) * frnom(k+1) 
     .                           / ( dn(k+1) + dn(k  ) ) * 2.D0
      dqbdlf00(k) = + dqbdq00(k) * frnom(k  ) 
     .                           / ( dn(k  ) + dn(k-1) ) * 2.D0
      dqbdlfm1(k) = + dqbdqm1(k) * frnom(k-1) 
     .                           / ( dn(k-1) + dn(k-2) ) * 2.D0
c
c     equations rm16 - rm19
c
      dqbdldp2(k) = - dqbdqp2(k) * qn(k+2) * dn(k+2)
     .                         / ( dn(k+2) + dn(k+1) )
c
      dqbdldp1(k) = - dqbdqp1(k) * qn(k+1) * dn(k+1)
     .                         / ( dn(k+1) + dn(k  ) )
     .              - dqbdqp2(k) * qn(k+2) * dn(k+1)
     .                         / ( dn(k+2) + dn(k+1) )
c
      dqbdld00(k) = - dqbdq00(k) * qn(k  ) * dn(k  )
     .                         / ( dn(k  ) + dn(k-1) )
     .              - dqbdqp1(k) * qn(k+1) * dn(k  )
     .                         / ( dn(k+1) + dn(k  ) )
c
      dqbdldm1(k) = - dqbdqm1(k) * qn(k-1) * dn(k-1)
     .                         / ( dn(k-1) + dn(k-2) )
     .              - dqbdq00(k) * qn(k  ) * dn(k-1)
     .                         / ( dn(k  ) + dn(k-1) )
c
      dqbdldm2(k) = - dqbdqm1(k) * qn(k-1) * dn(k-2)
     .                         / ( dn(k-1) + dn(k-2) )
    7 continue
c
c=======================================================================
c     interior zones
c=======================================================================
c
      do 8 k = ngrs + 1, ngre
c
c-----------------------------------------------------------------------
c     time derivative
c-----------------------------------------------------------------------
c
c     equations rm4 - rm6
c
      ep1(if, jr, k) = ep1(if, jr, k) + qn(k) *   dn(k  ) * rmup1n(k+1) 
     .                                         / (2.0D0 * cc**2 * dtime)
c
      e00(if, jr, k) = e00(if, jr, k) + qn(k) * ( dn(k-1) * rmup1n(k  ) 
     .                                          - dn(k  ) * rmup1n(k  ))
     .                                         / (2.0D0 * cc**2 * dtime)
c
      em1(if, jr, k) = em1(if, jr, k) - qn(k) *   dn(k-1) * rmup1n(k-1) 
     .                                         / (2.0D0 * cc**2 * dtime)
c
c     equations rm7 & rm8
c
      e00(if, jd, k) = e00(if, jd, k)   + qn(k) * (dn(k  ) * dvoln(k  ))
     .                                         / (2.0D0 * cc**2 * dtime)
     .  - (dn(k  ) / (dn(k-1) + dn(k))) * qn(k) * (dn(k-1) * dvoln(k-1) 
     .                                           + dn(k  ) * dvoln(k  ))
     .                                         / (2.0D0 * cc**2 * dtime)
c
      em1(if, jd, k) = em1(if, jd, k)   + qn(k) * (dn(k-1) * dvoln(k-1))
     .                                         / (2.0D0 * cc**2 * dtime)
     .  - (dn(k-1) / (dn(k-1) + dn(k))) * qn(k) * (dn(k-1) * dvoln(k-1) 
     .                                           + dn(k  ) * dvoln(k  ))
     .                                         / (2.0D0 * cc**2 * dtime)
c
c     equation rm9
c
      e00(if, jf, k) = e00(if, jf, k)   + ( frnom(k) / (cc**2 * dtime) )
     .                                          * (dn(k-1) * dvoln(k-1)
     .                                          +  dn(k  ) * dvoln(k  ))
     .                                          / (dn(k-1) + dn   (k  ))
c
c-----------------------------------------------------------------------
c     advection
c-----------------------------------------------------------------------
c
c     equations rm20 - rm21
c
      ep1(if, jm, k) = ep1(if, jm, k) 
     .                 + 0.5D0 * xmen(k+1) * frb(k  ) / (cc**2 * dtime)
c
      e00(if, jm, k) = e00(if, jm, k) 
     .                 + 0.5D0 * xmen(k  ) * frb(k  ) / (cc**2 * dtime)
     .                 - 0.5D0 * xmen(k  ) * frb(k-1) / (cc**2 * dtime)
c
      em1(if, jm, k) = em1(if, jm, k) 
     .                 - 0.5D0 * xmen(k-1) * frb(k-1) / (cc**2 * dtime)
c
c     equations rm23 - rm27
c
      ep2(if, jd, k) = ep2(if, jd, k) 
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdldp2(k  ) / cc**2
c
      ep1(if, jd, k) = ep1(if, jd, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdldp2(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdldp1(k  ) / cc**2
c
      e00(if, jd, k) = e00(if, jd, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdldp1(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdld00(k  ) / cc**2
c
      em1(if, jd, k) = em1(if, jd, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdld00(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdldm1(k  ) / cc**2
c
      em2(if, jd, k) = em2(if, jd, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdldm1(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdldm2(k  ) / cc**2
c
c-----em3(if, jd, k) = em3(if, jd, k) 
c----.         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdldm2(k-1) / cc**2
c
c     equations rm28 - rm32
c
      ep2(if, jf, k) = ep2(if, jf, k) 
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdlfp2(k  ) / cc**2
c
      ep1(if, jf, k) = ep1(if, jf, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdlfp2(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdlfp1(k  ) / cc**2
c
      e00(if, jf, k) = e00(if, jf, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdlfp1(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdlf00(k  ) / cc**2
c
      em1(if, jf, k) = em1(if, jf, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdlf00(k-1) / cc**2
     .         - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdlfm1(k  ) / cc**2
c
      em2(if, jf, k) = em2(if, jf, k) 
     .         + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdlfm1(k-1) / cc**2
c
c-----------------------------------------------------------------------
c     pressure gradient
c-----------------------------------------------------------------------
c
c     equation rm33
c
      e00(if, jr, k) = e00(if, jr, k)  
     .         + thet * xmu * rn(k) * rmum1(k) * ( fedd(k  ) * er(k  ) 
     .                                           - fedd(k-1) * er(k-1) )
     .
     .         + thet * xmu * rn(k) * rmum1(k) * ( fr  (k  ) / cc**2 ) 
     .                                              * 0.5D0 * ( u(k+1)
     .                                                        - u(k-1) )
c
c     equations rm34 & rm35
c
      ep1(if, ju, k) = ep1(if, ju, k) 
     .             + thet * rmu(k) * (fr(k) / cc**2) * 0.5D0 * unom(k+1)
c
      em1(if, ju, k) = em1(if, ju, k) 
     .             - thet * rmu(k) * (fr(k) / cc**2) * 0.5D0 * unom(k-1)
c
c     equations rm36 & rm37
c
      e00(if, je, k) = e00(if, je, k) 
     .                            + thet * rmu(k) * fedd(k  ) * ern(k  ) 
c
      em1(if, je, k) = em1(if, je, k) 
     .                            - thet * rmu(k) * fedd(k-1) * ern(k-1) 
c
c     equation rm38
c
      e00(if, jf, k) = e00(if, jf, k) 
     .         + thet * rmu(k) * (frnom(k) / cc**2) * 0.5D0 * ( u(k+1)
     .                                                        - u(k-1) )
c
c-----------------------------------------------------------------------
c     isotropy
c-----------------------------------------------------------------------
c
c     equations rm39 - rm41
c
      ep1(if, jr, k) = ep1(if, jr, k)
     .  + 0.25D0 * thet * xmu * rmu  (k+1) * (rn(k+1) / r(k))
     .                          * (3.D0 * fedd(k  ) - 1.D0) * er (k  ) 
c
      e00(if, jr, k) = e00(if, jr, k) 
     .  + 0.25D0 * thet * xmu * rmum1(k  ) *  rn(k  )
     .                        * ( (3.D0 * fedd(k-1) - 1.D0) * er (k-1)
     .                          - (3.D0 * fedd(k  ) - 1.D0) * er (k  ) ) 
c
      em1(if, jr, k) = em1(if, jr, k) 
     .  - 0.25D0 * thet * xmu * rmu  (k-1) * (rn(k-1) / r(k))
     .                          * (3.D0 * fedd(k-1) - 1.D0) * er (k-1)
c
      e00(if, jr, k) = e00(if, jr, k) 
     .  - 0.25D0 * thet * xmu * (rn(k  ) / r(k)**2)
     .            * ( dvol(k-1) * (3.D0 * fedd(k-1) - 1.D0) * er (k-1)
     .              + dvol(k  ) * (3.D0 * fedd(k  ) - 1.D0) * er (k  ) ) 
c
c     equations rm42 & rm43
c
      e00(if, je, k) = e00(if, je, k)
     .  + 0.25D0 * thet * xmu * (dvol(k  ) / r(k))
     .                          * (3.D0 * fedd(k  ) - 1.D0) * ern(k  )
c
      em1(if, je, k) = em1(if, je, k) 
     .  + 0.25D0 * thet * xmu * (dvol(k-1) / r(k))
     .                          * (3.D0 * fedd(k-1) - 1.D0) * ern(k-1)
c
c-----------------------------------------------------------------------
c     radiation force
c-----------------------------------------------------------------------
c
c     equations rm45 - rm47
c
      ep1(if, jr, k) = ep1(if, jr, k)        + avchi(k) * (fr  (k) / cc)
     .                 * 0.5D0 * thet * (+ d(k  ) * rmu(k+1) ) * rn(k+1) 
c
c
      e00(if, jr, k) = e00(if, jr, k)        + avchi(k) * (fr  (k) / cc)
     .                 * 0.5D0 * thet * (- d(k  ) * rmu(k  ) 
     .                                   + d(k-1) * rmu(k  ) ) * rn(k  )
c
      em1(if, jr, k) = em1(if, jr, k)        + avchi(k) * (fr  (k) / cc)
     .                 * 0.5D0 * thet * (- d(k-1) * rmu(k-1) ) * rn(k-1) 
c
c     equations rm48 & rm49
c
      e00(if, jd, k) = e00(if, jd, k)        + avchi(k) * (fr  (k) / cc)
     .        * 0.5D0 * thet                         * dn(k) * dvol(k)
     .
     .        +               (avchi(k) / chif(k  ))**2 * (fr  (k) / cc)
     .        * 0.5D0 * thet * ( d (k-1) * dvol(k-1) + d (k) * dvol(k) )
     .        * 0.5D0 * chifn(k  ) * dlcfdldn(k  ) 
c
      em1(if, jd, k) = em1(if, jd, k)        + avchi(k) * (fr  (k) / cc)
     .        * 0.5D0 * thet *   dn(k-1) * dvol(k-1)   
     .
     .        +               (avchi(k) / chif(k-1))**2 * (fr  (k) / cc)
     .        * 0.5D0 * thet * ( d (k-1) * dvol(k-1) + d (k) * dvol(k) )
     .        * 0.5D0 * chifn(k-1) * dlcfdldn(k-1) 
c
c     equations rm50 & rm51
c
      e00(if, jt, k) = e00(if, jt, k) 
     .        +               (avchi(k) / chif(k  ))**2 * (fr  (k) / cc)
     .        * 0.5D0 * thet * ( d (k-1) * dvol(k-1) + d (k) * dvol(k) )
     .        * 0.5D0 * chifn(k  ) * dlcfdltn(k  ) 
c
      em1(if, jt, k) = em1(if, jt, k) 
     .        +               (avchi(k) / chif(k-1))**2 * (fr  (k) / cc)
     .        * 0.5D0 * thet * ( d (k-1) * dvol(k-1) + d (k) * dvol(k) )
     .        * 0.5D0 * chifn(k-1) * dlcfdltn(k-1) 
c
c     equation rm52
c
      e00(if, jf, k) = e00(if, jf, k)       + avchi(k) * (frnom(k) / cc)
     .        * 0.5D0 * thet * ( d (k-1) * dvol(k-1) + d (k) * dvol(k) )
c
c-----------------------------------------------------------------------
c     right hand side
c-----------------------------------------------------------------------
c
c     equation rm53
c
      rhs(if, k) = ( qn(k) * (dn(k-1) * dvoln(k-1) + dn(k) * dvoln(k))
     .             - qo(k) * (do(k-1) * dvolo(k-1) + do(k) * dvolo(k)) )
     .                                         / (2.0D0 * cc**2 * dtime)
     .
     .                        - ( ( dmdt(k  ) + dmdt(k+1) ) * frb(k  )
     .                          - ( dmdt(k-1) + dmdt(k  ) ) * frb(k-1) )
     .                                         / (2.0D0 * cc**2        )
     .
     .                              + rmu(k) * ( fedd(k  ) *  er(k  ) 
     .                                         - fedd(k-1) *  er(k-1) )
     .
     .                   + 0.5D0 * rmu(k) * (fr(k) / cc**2) * ( u(k+1)
     .                                                      -   u(k-1) )
     .
     .           + ( 0.25D0 * xmu / r(k) )
     .           * ( (3.0D0 * fedd(k-1) - 1.0D0) * er(k-1) * dvol(k-1)
     .             + (3.0D0 * fedd(k  ) - 1.0D0) * er(k  ) * dvol(k  ) )
     .
     .               + 0.5D0 * ( d (k) * dvol(k) + d (k-1) * dvol(k-1) )
     .               *  avchi(k) * fr(k) / cc 
    8 continue
c
      end if
c
c++++++++ ltran .le. 1 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c=======================================================================
c     inner boundary condition
c=======================================================================
c
      k = ngrs
c
c     optically transmitting boundary; equation bc22
c
      if (lribc .eq. 1) then
            rhs(if,     k) = - cpi * ( 2.0D0 * geddl + 1.0D0 ) * xipl 
     .                     + frn  (k)   + cc * geddl * ern(k)
            e00(if, je, k) =            + cc * geddl * ern(k)
            e00(if, jf, k) = frnom(k)
      end if
c
c     optically reflecting boundary; equation bc23
c
      if (lribc .eq. 2) then
            rhs(if,     k) = 0.0D0
            e00(if, jf, k) = 1.0D0
      end if
c
c     imposed net flux; equations bc24 & bc25
c
      if (lribc .eq. 3) then
            rhs(if,     k) = 0.0D0
            e00(if, jf, k) = 1.0D0
      end if
c
c=======================================================================
c     outer boundary condition
c=======================================================================
c
      k = ngre + 1
c
c     optically transmitting boundary; equation bc35
c
      if (lrobc .eq. 1) then
            rhs(if,     k) = + cpi * ( 2.0D0 * geddr + 1.0D0 ) * ximr 
     .                     + frn  (k)   - cc * geddr * ern(k-1) 
            e00(if, jf, k) = frnom(k)
            em1(if, je, k) =            - cc * geddr * ern(k-1)
      end if
c
c     optically reflecting boundary; equation bc36
c
      if (lrobc .eq. 2) then
            rhs(if,     k) = 0.0D0
            e00(if, jf, k) = 1.0D0
      end if
c
c     imposed net flux; equation bc37
c
      if (lrobc .eq. 3) then
            rhs(if,     k) = rmun(k-1) * frn  (k-1) - rmun(k) * frn  (k)
            e00(if, jr, k) =                  - xmu * rmun(k) * frn  (k)
            em1(if, jr, k) = rmun(k-1) * frn  (k-1) * xmu
c
            e00(if, jf, k) =                        - rmun(k) * frnom(k)
            em1(if, jf, k) = rmun(k-1) * frnom(k-1)
      end if
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 9 k = 1, ngrs - 2
      rhs(if,     k) = 0.0D0
      e00(if, jf, k) = 1.0D0
    9 continue
c
      do 10 k = ngre + 3, mgr
      rhs(if,     k) = 0.0D0
      e00(if, jf, k) = 1.0D0
   10 continue
c
c-----------------------------------------------------------------------
c     inner boundary
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
c     optically transmitting boundary; equations pz33 & pz34
c
      if (lribc .eq. 1) then
            rhs(if,     k) = frn  (k+1) - frn  (k)
            ep1(if, jf, k) = frnom(k+1)
            e00(if, jf, k) =            - frnom(k)
      end if
c
c     optically reflecting boundary; equations pz35 & pz36
c
      if (lribc .eq. 2) then
            rhs(if,     k) = frn  (k+2) + frn  (k)
            ep2(if, jf, k) = frnom(k+2)
            e00(if, jf, k) =            + frnom(k)
      end if
c
c     imposed net flux; equations pz37 & pz38
c
      if (lribc .eq. 3) then
            rhs(if,     k) = rmun(k+1) * frn  (k+1) - rmun(k) * frn  (k)
            ep1(if, jr, k) = rmun(k+1) * frn  (k+1) * xmu
            e00(if, jr, k) =                  - xmu * rmun(k) * frn  (k)
            ep1(if, jf, k) = rmun(k+1) * frnom(k+1)
            e00(if, jf, k) =                        - rmun(k) * frnom(k)
      end if
c
c-----------------------------------------------------------------------
c     outer boundary 
c-----------------------------------------------------------------------
c
      k = ngre + 2
c
c     optically transmitting boundary; equations pz81 & pz82
c
      if (lrobc .eq. 1) then
            rhs(if,     k) = - frn  (k-1) + frn  (k)
            e00(if, jf, k) =              + frnom(k)
            em1(if, jf, k) = - frnom(k-1)
      end if
c
c     optically reflecting boundary; equations pz85 & pz86
c
      if (lrobc .eq. 2) then
            rhs(if,     k) = + frn  (k-2) + frn  (k)
            e00(if, jf, k) =              + frnom(k)
            em2(if, jf, k) = + frnom(k-2)
      end if
c
c     imposed net flux; equations pz89 & pz90
c
      if (lrobc .eq. 3) then
            rhs(if,     k) = rmun(k-1) * frn  (k-1) - rmun(k) * frn  (k)
            e00(if, jr, k) =                  - xmu * rmun(k) * frn  (k)
            em1(if, jr, k) = rmun(k-1) * frn  (k-1) * xmu
            e00(if, jf, k) =                        - rmun(k) * frnom(k)
            em1(if, jf, k) = rmun(k-1) * frnom(k-1)
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
