      subroutine gasnrgrh
c
c***********************************************************************
c     radiating fluid energy equation fe3
c.......................................................................
c     calling sequence: matgen > gasnrgh > advectc, diffuse
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq1'
      include'titan.eq3'
c
c=======================================================================
c     define various quantities and their derivatives
c=======================================================================
c
c-----------------------------------------------------------------------
c     f1 = <u/r>; equations re16 - re18
c-----------------------------------------------------------------------
c
      f1  (xk,xp,uk,up) =   0.5D0*(  up/xp +   uk/xk   )
      f1rk(xk,xp,uk,up) = - 0.5D0*(        +   uk/xk**2)
      f1rp(xk,xp,uk,up) = - 0.5D0*(  up/xp**2          )
      f1uk(xk,xp,uk,up) = + 0.5D0*(        + 1.D0/xk   )
      f1up(xk,xp,uk,up) = + 0.5D0*(1.D0/xp             )
c
c=======================================================================
c     set auxiliary quantities for radiation anisotropy, 
c     equations re19 -  re23
c=======================================================================
c
      do 1 k = ngrs, ngre
c
       aur    (k) = f1  (r(k),r(k+1),u(k),u(k+1))
      dardlrp1(k) = f1rp(r(k),r(k+1),u(k),u(k+1)) * thet * rn  (k+1)
      dardlr00(k) = f1rk(r(k),r(k+1),u(k),u(k+1)) * thet * rn  (k  )
      dardlup1(k) = f1up(r(k),r(k+1),u(k),u(k+1)) * thet * unom(k+1)
      dardlu00(k) = f1uk(r(k),r(k+1),u(k),u(k+1)) * thet * unom(k+1)
    1 continue
c
c=======================================================================
c     set interface energy densities for advection
c=======================================================================
c
      do 3 l = 1, madv
      do 2 k = 1, mgr
      adv(k, l) = 0.0D0
    2 continue
    3 continue
c
      do 4   k  = ngrs - 1, ngre + 2
      qn    (k) = egrn (k)
      qo    (k) = egro (k)
    4 continue
      do 44  k  = ngrs    , ngre + 1
      qso   (k) = egrso(k)
   44 continue
      do 444 k  = ngrs + 1, ngre + 1
      flow  (k) = - dmdt(k)
  444 continue
c
      call advectc
c
      do 5  k  = ngrs + 1, ngre + 1
      egrb (k) = qb (k)
    5 continue
      do 55 k  = ngrs    , ngre + 1
      egrsn(k) = qsn(k)
   55 continue
c
c-----------------------------------------------------------------------
c     set auxiliary derivatives for advection
c-----------------------------------------------------------------------
c
      do 6 k = ngrs - 1, ngre + 2
c
c     equations fe9 - fe11
c
      dlqdlt(k) = ( egn(k) * dlegdltn(k)                  ) / egrn(k)
      dlqdld(k) = ( egn(k) * dlegdldn(k) - ern(k) / dn(k) ) / egrn(k)
      dlqdle(k) = (                      + ern(k) / dn(k) ) / egrn(k)
    6 continue
c
      do 7 k = ngrs + 1, ngre + 1
c
c     equations fe12 - fe23
c
      dqbdlep1(k) = dqbdlqp1(k) * dlqdle(k+1)
      dqbdle00(k) = dqbdlq00(k) * dlqdle(k  )
      dqbdlem1(k) = dqbdlqm1(k) * dlqdle(k-1)
      dqbdlem2(k) = dqbdlqm2(k) * dlqdle(k-2)
c
      dqbdltp1(k) = dqbdlqp1(k) * dlqdlt(k+1)
      dqbdlt00(k) = dqbdlq00(k) * dlqdlt(k  )
      dqbdltm1(k) = dqbdlqm1(k) * dlqdlt(k-1)
      dqbdltm2(k) = dqbdlqm2(k) * dlqdlt(k-2)
c
      dqbdldp1(k) = dqbdlqp1(k) * dlqdld(k+1)
      dqbdld00(k) = dqbdlq00(k) * dlqdld(k  )
      dqbdldm1(k) = dqbdlqm1(k) * dlqdld(k-1)
      dqbdldm2(k) = dqbdlqm2(k) * dlqdld(k-2)
    7 continue
c
c=======================================================================
c     set gas energy diffusion terms
c=======================================================================
c
      do 9 l = 1, mdfz
      do 8 k = 1, mgr
      dfz(k, l) = 0.0D0
    8 continue
    9 continue
c
      sig = sige
      zet = 1.0D0
c
      call diffuse
c
c=======================================================================
c     all terms except advection and rhs, note that this loop covers the
c     entire domain, including boundaries
c=======================================================================
c
      do 10 k = ngrs, ngre
c
c-----------------------------------------------------------------------
c     time derivative
c-----------------------------------------------------------------------
c
c     equations fe4 - fe8
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .                           + dn(k) * egrn(k) * rmup1n(k+1) / dtime
      e00(it, jr, k) = e00(it, jr, k)
     .                           - dn(k) * egrn(k) * rmup1n(k  ) / dtime
c
      e00(it, jd, k) = e00(it, jd, k)
     .       + (1.0D0 + dlegdldn(k)) * dn(k) * egn(k) * dvoln(k) / dtime
c
      e00(it, jt, k) = e00(it, jt, k)
     .                + dlegdltn(k)  * dn(k) * egn(k) * dvoln(k) / dtime
c
      e00(it, je, k) = e00(it, je, k)        + ern(k) * dvoln(k) / dtime
c
c-----------------------------------------------------------------------
c     flux divergence
c-----------------------------------------------------------------------
c
c     equations re6 - re7
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .                  + thet * xmu * rn(k+1) * rmum1(k+1) * fr   (k+1)
      e00(it, jr, k) = e00(it, jr, k) 
     .                  - thet * xmu * rn(k  ) * rmum1(k  ) * fr   (k  )
c
c     equations re8 - re9
c
      ep1(it, jf, k) = ep1(it, jf, k)   + thet * rmu  (k+1) * frnom(k+1)
      e00(it, jf, k) = e00(it, jf, k)   - thet * rmu  (k  ) * frnom(k  )
c
c-----------------------------------------------------------------------
c     work
c-----------------------------------------------------------------------
c
c     equations fe41 - fe47
c
      ep1(it, jr, k) = ep1(it, jr, k)  + thet * xmu * rn(k+1)
     .             * ( pg(k) + fedd(k) * er(k) ) * rmum1(k+1) * u(k+1)
      e00(it, jr, k) = e00(it, jr, k)  - thet * xmu * rn(k  )
     .             * ( pg(k) + fedd(k) * er(k) ) * rmum1(k  ) * u(k  )
c
      e00(it, jd, k) = e00(it, jd, k) 
     .             + thet * pgn(k) * dlpgdldn(k) * ( rmu(k+1) * u(k+1) 
     .                                             - rmu(k  ) * u(k  ))
c
      e00(it, jt, k) = e00(it, jt, k) 
     .             + thet * pgn(k) * dlpgdltn(k) * ( rmu(k+1) * u(k+1) 
     .                                             - rmu(k  ) * u(k  ))
c
      e00(it, je, k) = e00(it, je, k)
     .             + thet   * fedd(k) * ern(k)   * ( rmu(k+1) * u(k+1) 
     .                                             - rmu(k  ) * u(k  ))
c
      ep1(it, ju, k) = ep1(it, ju, k)
     .     + thet * ( pg(k) + fedd(k) * er (k) ) * rmu(k+1) * unom(k+1)
      e00(it, ju, k) = e00(it, ju, k) 
     .     - thet * ( pg(k) + fedd(k) * er (k) ) * rmu(k  ) * unom(k  )
c
c-----------------------------------------------------------------------
c     radiation anisotropy
c-----------------------------------------------------------------------
c
c     equations re24 - re28
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     . *( + thet * rn(k+1) * rmu(k+1) * aur(k) +  dardlrp1(k) * dvol(k))
      e00(it, jr, k) = e00(it, jr, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     . *( - thet * rn(k  ) * rmu(k  ) * aur(k) +  dardlr00(k) * dvol(k))
c
      ep1(it, ju, k) = ep1(it, ju, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     .                                         *  dardlup1(k) * dvol(k) 
      e00(it, ju, k) = e00(it, ju, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     .                                         *  dardlu00(k) * dvol(k) 
c
      e00(it, je, k) = e00(it, je, k) 
     .           + thet * 0.5D0 * xmu * ern(k) * (1.0D0 - 3.0D0*fedd(k))
     .                                * aur(k) *                dvol(k)
c
c-----------------------------------------------------------------------
c     (artificial) energy diffusion
c-----------------------------------------------------------------------
c
c     equations fe48 - fe57
c
      ep2(it, jr, k) = ep2(it, jr, k)                    - ddfdlrp1(k+1)
      ep1(it, jr, k) = ep1(it, jr, k)    + ddfdlrp1(k  ) - ddfdlr00(k+1)
      e00(it, jr, k) = e00(it, jr, k)    + ddfdlr00(k  ) - ddfdlrm1(k+1)
      em1(it, jr, k) = em1(it, jr, k)    + ddfdlrm1(k  )
c
      ep1(it, jd, k) = ep1(it, jd, k)                    - ddfdld00(k+1)
      e00(it, jd, k) = e00(it, jd, k)    + ddfdld00(k  ) - ddfdldm1(k+1)
      em1(it, jd, k) = em1(it, jd, k)    + ddfdldm1(k  ) 
c
      ep1(it, jt, k) = ep1(it, jt, k)                    - ddfdlt00(k+1)
      e00(it, jt, k) = e00(it, jt, k)    + ddfdlt00(k  ) - ddfdltm1(k+1)
      em1(it, jt, k) = em1(it, jt, k)    + ddfdltm1(k  )
c
c-----------------------------------------------------------------------
c     (artificial) energy dissipation
c-----------------------------------------------------------------------
c
c     equations fe91 - fe96
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .                     + qf(k) * dudr(k) * thet * rn(k+1) * rmu(k+1)
     .                             + dqfdlrp1(k) *  dudr   (k) * dvol(k)
     .                             +  qf     (k) * durdlrp1(k) * dvol(k)
      e00(it, jr, k) = e00(it, jr, k) 
     .                     - qf(k) * dudr(k) * thet * rn(k  ) * rmu(k  )
     .                             + dqfdlr00(k) *  dudr   (k) * dvol(k)
     .                             +  qf     (k) * durdlr00(k) * dvol(k)
c
      ep1(it, ju, k) = ep1(it, ju, k) 
     .                             + dqfdlup1(k) *  dudr   (k) * dvol(k)
     .                             +  qf     (k) * durdlup1(k) * dvol(k)
      e00(it, ju, k) = e00(it, ju, k) 
     .                             + dqfdlu00(k) *  dudr   (k) * dvol(k)
     .                             +  qf     (k) * durdlu00(k) * dvol(k)
c
      e00(it, jt, k) = e00(it, jt, k) 
     .                             + dqfdlt00(k) *  dudr   (k) * dvol(k)
c
      e00(it, jd, k) = e00(it, jd, k) 
     .                             + dqfdld00(k) *  dudr   (k) * dvol(k)
   10 continue
c
c=======================================================================
c     advection and right-hand side
c=======================================================================
c
c-----------------------------------------------------------------------
c     interior zones
c-----------------------------------------------------------------------
c
      do 11 k = ngrs + 1, ngre - 1
c
c.......................................................................
c     advection
c.......................................................................
c
c     equations fe24 & fe25
c
      ep1(it, jm, k) = ep1(it, jm, k)   + egrb(k+1) * xmen(k+1) / dtime
      e00(it, jm, k) = e00(it, jm, k)   - egrb(k  ) * xmen(k  ) / dtime
c
c     equations fe26 - fe40
c
      ep2(it, jd, k) = ep2(it, jd, k)        - dmdt(k+1) * dqbdldp1(k+1)
      ep1(it, jd, k) = ep1(it, jd, k)        + dmdt(k  ) * dqbdldp1(k  )
     .                                       - dmdt(k+1) * dqbdld00(k+1)
      e00(it, jd, k) = e00(it, jd, k)        + dmdt(k  ) * dqbdld00(k  )
     .                                       - dmdt(k+1) * dqbdldm1(k+1)
      em1(it, jd, k) = em1(it, jd, k)        + dmdt(k  ) * dqbdldm1(k  )
     .                                       - dmdt(k+1) * dqbdldm2(k+1)
      em2(it, jd, k) = em2(it, jd, k)        + dmdt(k  ) * dqbdldm2(k  )
c
      ep2(it, jt, k) = ep2(it, jt, k)        - dmdt(k+1) * dqbdltp1(k+1)
      ep1(it, jt, k) = ep1(it, jt, k)        + dmdt(k  ) * dqbdltp1(k  )
     .                                       - dmdt(k+1) * dqbdlt00(k+1)
      e00(it, jt, k) = e00(it, jt, k)        + dmdt(k  ) * dqbdlt00(k  )
     .                                       - dmdt(k+1) * dqbdltm1(k+1)
      em1(it, jt, k) = em1(it, jt, k)        + dmdt(k  ) * dqbdltm1(k  )
     .                                       - dmdt(k+1) * dqbdltm2(k+1)
      em2(it, jt, k) = em2(it, jt, k)        + dmdt(k  ) * dqbdltm2(k  )
c
      ep2(it, je, k) = ep2(it, je, k)        - dmdt(k+1) * dqbdlep1(k+1)
      ep1(it, je, k) = ep1(it, je, k)        + dmdt(k  ) * dqbdlep1(k  )
     .                                       - dmdt(k+1) * dqbdle00(k+1)
      e00(it, je, k) = e00(it, je, k)        + dmdt(k  ) * dqbdle00(k  )
     .                                       - dmdt(k+1) * dqbdlem1(k+1)
      em1(it, je, k) = em1(it, je, k)        + dmdt(k  ) * dqbdlem1(k  )
     .                                       - dmdt(k+1) * dqbdlem2(k+1)
      em2(it, je, k) = em2(it, je, k)        + dmdt(k  ) * dqbdlem2(k  )
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation fe97
c
      rhs(it, k) =                  ( dn(k) * egrn(k) * dvoln(k) 
     .                              - do(k) * egro(k) * dvolo(k) )/dtime
     .
     .                                         - dmdt(k+1) * egrb(k+1) 
     .                                         + dmdt(k  ) * egrb(k  )
     .
     .                                          + rmu(k+1) *  fr (k+1)
     .                                          - rmu(k  ) *  fr (k  ) 
     .
     .                                                     -  df (k+1)
     .                                                     +  df (k  )
     .
     .          + ( pg(k) + er(k) * fedd(k) ) * ( rmu(k+1) *   u (k+1) 
     .                                          - rmu(k  ) *   u (k  ) )
     .
     .          + ( 1.0D0 - 3.0D0 * fedd(k) )
     .            * 0.5D0 *   xmu *   er(k)   *   aur(k  ) * dvol(k  )
     .
     .                              + qf(k)   *  dudr(k  ) * dvol(k  )
   11 continue
c
c-----------------------------------------------------------------------
c     inner boundary condition, all cases; equation bc17
c-----------------------------------------------------------------------
c
      k = ngrs
c
c.......................................................................
c     advection; equation fe25 & fe27 - fe30 & fe32 - fe35 & fe37 - fe40
c.......................................................................
c
      ep1(it, jm, k) = ep1(it, jm, k)   + xmen(k+1) * egrb(k+1) / dtime
c
      ep2(it, jd, k) = ep2(it, jd, k)        - dmdt(k+1) * dqbdldp1(k+1)
      ep1(it, jd, k) = ep1(it, jd, k)        - dmdt(k+1) * dqbdld00(k+1)
      e00(it, jd, k) = e00(it, jd, k)        - dmdt(k+1) * dqbdldm1(k+1)
      em1(it, jd, k) = em1(it, jd, k)        - dmdt(k+1) * dqbdldm2(k+1)
c
      ep2(it, jt, k) = ep2(it, jt, k)        - dmdt(k+1) * dqbdltp1(k+1)
      ep1(it, jt, k) = ep1(it, jt, k)        - dmdt(k+1) * dqbdlt00(k+1)
      e00(it, jt, k) = e00(it, jt, k)        - dmdt(k+1) * dqbdltm1(k+1)
      em1(it, jt, k) = em1(it, jt, k)        - dmdt(k+1) * dqbdltm2(k+1)
c
      ep2(it, je, k) = ep2(it, je, k)        - dmdt(k+1) * dqbdlep1(k+1)
      ep1(it, je, k) = ep1(it, je, k)        - dmdt(k+1) * dqbdle00(k+1)
      e00(it, je, k) = e00(it, je, k)        - dmdt(k+1) * dqbdlem1(k+1)
      em1(it, je, k) = em1(it, je, k)        - dmdt(k+1) * dqbdlem2(k+1)
c
      e00(it, jr, k) = e00(it, jr, k) 
     .                               + rmu  (k) * rn(k) * er(k) / dtime
     .                  - thet * xmu * rmum1(k) * rn(k)     
     .                                       * (urel(k) * er(k) + phil2)
      e00(it, je, k) = e00(it, je, k) 
     .                              - thet * rmu(k) * urel(k) * ern(k)
      e00(it, ju, k) = e00(it, ju, k) 
     .                              - thet * rmu(k) * unom(k) * er (k)
c
c.......................................................................
c     right-hand side; equation bc17
c.......................................................................
c
      rhs(it, k) =                  ( dn(k) * egrn(k) * dvoln(k) 
     .                              - do(k) * egro(k) * dvolo(k) )/dtime
     .
     .                                        -  dmdt(k+1) * egrb(k+1) 
     .          - ( phil2 + er(k) * urel(k) ) *   rmu(k  ) 
     .
     .                                          + rmu(k+1) *  fr (k+1)
     .                                          - rmu(k  ) *  fr (k  ) 
     .
     .                                                     -  df (k+1)
     .                                                     +  df (k  )
     .
     .          + ( pg(k) + er(k) * fedd(k) ) * ( rmu(k+1) *   u (k+1) 
     .                                          - rmu(k  ) *   u (k  ) )
     .
     .          + ( 1.0D0 - 3.0D0 * fedd(k) )
     .                * 0.5D0 * xmu * er(k)   *   aur(k)   * dvol(k  )
     .
     .                              + qf(k)   *  dudr(k)   * dvol(k  )
c
c-----------------------------------------------------------------------
c     outer boundary condition (all cases); equation bc33
c-----------------------------------------------------------------------
c
      k = ngre
c
c.......................................................................
c     advection; equation fe24 & fe26 - fe29 & fe31 - fe34 & fe36 - fe39
c.......................................................................
c
      e00(it, jm, k) = e00(it, jm, k)       - xmen(k) * egrb(k) / dtime
c
      ep1(it, jd, k) = ep1(it, jd, k)            + dmdt(k) * dqbdldp1(k)
      e00(it, jd, k) = e00(it, jd, k)            + dmdt(k) * dqbdld00(k)
      em1(it, jd, k) = em1(it, jd, k)            + dmdt(k) * dqbdldm1(k)
      em2(it, jd, k) = em2(it, jd, k)            + dmdt(k) * dqbdldm2(k)
c
      ep1(it, jt, k) = ep1(it, jt, k)            + dmdt(k) * dqbdltp1(k)
      e00(it, jt, k) = e00(it, jt, k)            + dmdt(k) * dqbdlt00(k)
      em1(it, jt, k) = em1(it, jt, k)            + dmdt(k) * dqbdltm1(k)
      em2(it, jt, k) = em2(it, jt, k)            + dmdt(k) * dqbdltm2(k)
c
      ep1(it, je, k) = ep1(it, je, k)            + dmdt(k) * dqbdlep1(k)
      e00(it, je, k) = e00(it, je, k)            + dmdt(k) * dqbdle00(k)
      em1(it, je, k) = em1(it, je, k)            + dmdt(k) * dqbdlem1(k)
      em2(it, je, k) = em2(it, je, k)            + dmdt(k) * dqbdlem2(k)
c
c     extra terms for lagrangean or nontransmitting eulerian boundary;
c     equation bc33
c
      if (leobc .ne. 3) then
c
            ep1(it, jr, k) = ep1(it, jr, k) 
     .                           - rmu  (k+1) * rn(k+1) * er(k) / dtime
     .              + thet * xmu * rmum1(k+1) * rn(k+1)
     .                                     * (urel(k+1) * er(k) + phir2)
            e00(it, je, k) = e00(it, je, k)
     .                         + thet * rmu(k+1) * urel(k+1) * ern(k)
            ep1(it, ju, k) = ep1(it, ju, k)
     .                         + thet * rmu(k+1) * unom(k+1) * er (k)
      end if
c
c     extra terms for eulerian transmitting outer boundary;
c     equation bc57
c
      if (leobc .eq. 3) then
c
            ep1(it, jm, k) = ep1(it, jm, k) 
     .                                  + xmen(k+1) * egrb(k+1) / dtime
c
            ep2(it, jd, k) = ep2(it, jd, k)  - dmdt(k+1) * dqbdldp1(k+1)
            ep1(it, jd, k) = ep1(it, jd, k)  - dmdt(k+1) * dqbdld00(k+1)
            e00(it, jd, k) = e00(it, jd, k)  - dmdt(k+1) * dqbdldm1(k+1)
            em1(it, jd, k) = em1(it, jd, k)  - dmdt(k+1) * dqbdldm2(k+1)
c
            ep2(it, jt, k) = ep2(it, jt, k)  - dmdt(k+1) * dqbdltp1(k+1)
            ep1(it, jt, k) = ep1(it, jt, k)  - dmdt(k+1) * dqbdlt00(k+1)
            e00(it, jt, k) = e00(it, jt, k)  - dmdt(k+1) * dqbdltm1(k+1)
            em1(it, jt, k) = em1(it, jt, k)  - dmdt(k+1) * dqbdltm2(k+1)
c
            ep2(it, je, k) = ep2(it, je, k)  - dmdt(k+1) * dqbdlep1(k+1)
            ep1(it, je, k) = ep1(it, je, k)  - dmdt(k+1) * dqbdle00(k+1)
            e00(it, je, k) = e00(it, je, k)  - dmdt(k+1) * dqbdlem1(k+1)
            em1(it, je, k) = em1(it, je, k)  - dmdt(k+1) * dqbdlem2(k+1)
c
      end if
c
c.......................................................................
c     right-hand side; equation bc33
c.......................................................................
c
      rhs(it, k) =                  ( dn(k) * egrn(k) * dvoln(k) 
     .                              - do(k) * egro(k) * dvolo(k) )/dtime
     .
     .                                         + dmdt(k  ) * egrb(k  )
     .
     .                                          + rmu(k+1) *  fr (k+1)
     .                                          - rmu(k  ) *  fr (k  ) 
     .
     .                                                     -  df (k+1)
     .                                                     +  df (k  )
     .
     .          + ( pg(k) + er(k) * fedd(k) ) * ( rmu(k+1) *   u (k+1) 
     .                                          - rmu(k  ) *   u (k  ) )
     .
     .          + ( 1.0D0 - 3.0D0 * fedd(k) )
     .            * 0.5D0 *   xmu *   er(k)   *   aur(k  ) * dvol(k  )
     .
     .                              + qf(k)   *  dudr(k  ) * dvol(k  )
c
c     extra terms for lagrangean or nontransmitting eulerian boundary;
c     equation bc33
c
      if (leobc .ne. 3) rhs(it, k) = rhs(it, k)
     .                        + rmu(k+1) * ( phir2 + urel(k+1) * er(k) )
c
c     extra terms for eulerian transmitting outer boundary;
c     equation bc57
c
      if (leobc .eq. 3) rhs(it, k) = rhs(it, k) - dmdt(k+1) * egrb(k+1)
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 12 k = 1, ngrs - 2
      rhs(it,     k) = 0.0D0
      e00(it, jt, k) = 1.0D0
   12 continue
c
      do 13 k = ngre + 2, mgr
      rhs(it,     k) = 0.0D0
      e00(it, jt, k) = 1.0D0
   13 continue
c
c-----------------------------------------------------------------------
c     eulerian inner boundary
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
c     zero flux; equations pz9 & pz10
c
      if (leibc .eq. 1) then
            rhs(it,     k) =   egn(k+1) - egn(k)
            ep1(it, jd, k) =   egn(k+1) * dlegdldn(k+1)
            e00(it, jd, k) = - egn(k  ) * dlegdldn(k  )
            ep1(it, jt, k) =   egn(k+1) * dlegdltn(k+1)
            e00(it, jt, k) = - egn(k  ) * dlegdltn(k  )
      end if
c
c     nonzero flux; equations pz19 & pz20
c
      if (leibc .eq. 2) then
            rhs(it,     k) = 0.0D0
            e00(it, jt, k) = 1.0D0
      end if
c
c-----------------------------------------------------------------------
c     lagrangean inner boundary
c-----------------------------------------------------------------------
c
c     equations pz29 & pz30
c
      if (llibc .eq. 1) then
            rhs(it,     k) =   egn(k+1) - egn(k)
            ep1(it, jd, k) =   egn(k+1) * dlegdldn(k+1)
            e00(it, jd, k) = - egn(k  ) * dlegdldn(k  )
            ep1(it, jt, k) =   egn(k+1) * dlegdltn(k+1)
            e00(it, jt, k) = - egn(k  ) * dlegdltn(k  )
      end if
c
c-----------------------------------------------------------------------
c     eulerian outer boundary 
c-----------------------------------------------------------------------
c
c     zero flux; equations pz47 & pz48
c
      if (leobc .eq. 1) then
	    do 14 k = ngre + 1, ngre + 2
            rhs(it,     k) = - egn(k-1) + egn(k)
            e00(it, jd, k) =   egn(k  ) * dlegdldn(k  ) 
            em1(it, jd, k) = - egn(k-1) * dlegdldn(k-1)
            e00(it, jt, k) =   egn(k  ) * dlegdltn(k  ) 
            em1(it, jt, k) = - egn(k-1) * dlegdltn(k-1)
   14       continue
      end if
c
c     nonzero flux; equations pz57 & pz58
c
      if (leobc .eq. 2) then
	    do 15 k = ngre + 1, ngre + 2
            rhs(it,     k) = 0.0D0
            e00(it, jt, k) = 1.0D0
   15       continue
      end if
c
c     transmitting; equations pz67 & pz68 
c
      if (leobc .eq. 3) then
	    do 16 k = ngre + 1, ngre + 2
            rhs(it,     k) = - egn(k-1) + egn(k)
            e00(it, jd, k) =   egn(k  ) * dlegdldn(k  ) 
            em1(it, jd, k) = - egn(k-1) * dlegdldn(k-1)
            e00(it, jt, k) =   egn(k  ) * dlegdltn(k  ) 
            em1(it, jt, k) = - egn(k-1) * dlegdltn(k-1)
   16       continue
      end if
c
c-----------------------------------------------------------------------
c     lagrangean outer boundary
c-----------------------------------------------------------------------
c
c     equations pz77 & pz78
c
      if (llobc .gt. 0) then
	    do 17 k = ngre + 1, ngre + 2
            rhs(it,     k) = - egn(k-1) + egn(k)
            e00(it, jd, k) =   egn(k  ) * dlegdldn(k  ) 
            em1(it, jd, k) = - egn(k-1) * dlegdldn(k-1)
            e00(it, jt, k) =   egn(k  ) * dlegdltn(k  ) 
            em1(it, jt, k) = - egn(k-1) * dlegdltn(k-1)
   17       continue
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
