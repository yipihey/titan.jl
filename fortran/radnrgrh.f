      subroutine radnrgrh
c
c***********************************************************************
c     radiation energy equation re2
c.......................................................................
c     calling sequence: matgen > radnrgrh > advectc
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq1'
c
c=======================================================================
c     set auxiliary quantities for advection
c=======================================================================
c
      do 2 l = 1, madv
      do 1 k = 1, mgr
      adv(k, l) = 0.0D0
    1 continue
    2 continue
c
c***********************************************************************
c     equilibrium diffusion; equation rd5
c***********************************************************************
c
c++++++++ ltran .ge. 3 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (ltran .ge. 3) then
c
      do 4 k = ngrs, ngre
c
c     equations rd44 & rd45
c
      e00(ie, jt, k) = -4.0D0 * car * tn(k)**4
      e00(ie, je, k) =                           ern(k)
c
c     equation rd46
c
      rhs(ie,     k) =        - car * tn(k)**4 + ern(k)
    4 continue
c
      end if
c
c++++++++ ltran .ge. 3 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c***********************************************************************
c     radiation transport or nonequilibrium diffusion; equation re2
c***********************************************************************
c
c++++++++ ltran .lt. 3 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (ltran .lt. 3) then
c
c=======================================================================
c     set interface energy densities for advection
c=======================================================================
c
      do 5   k  = ngrs - 1, ngre + 2
      qn    (k) = ern (k) / dn(k)
      qo    (k) = ero (k) / do(k)
    5 continue
      do 55  k  = ngrs    , ngre + 1
      qso   (k) = erso(k)
   55 continue
      do 555 k  = ngrs + 1, ngre + 1
      flow  (k) = - dmdt(k)
  555 continue
c
      call advectc
c
      do 6  k  = ngrs + 1, ngre + 1
      erb  (k) = qb (k)
    6 continue
      do 66 k  = ngrs    , ngre + 1
      ersn (k) = qsn(k)
   66 continue
c
c=======================================================================
c     all terms except advection and rhs; note that this loop covers 
c     the entire domain including boundaries
c=======================================================================
c
      do 7 k = ngrs, ngre
c
c-----------------------------------------------------------------------
c     time derivative
c-----------------------------------------------------------------------
c
c     equations re3 - re5
c
      ep1(ie, jr, k) = ep1(ie, jr, k)    + ern(k) * rmup1n(k+1) / dtime
      e00(ie, jr, k) = e00(ie, jr, k)    - ern(k) * rmup1n(k  ) / dtime
c
      e00(ie, je, k) = e00(ie, je, k)    + ern(k) * dvoln (k  ) / dtime
c
c-----------------------------------------------------------------------
c     flux divergence
c-----------------------------------------------------------------------
c
c     equations re6 - re9
c
      ep1(ie, jr, k) = ep1(ie, jr, k) 
     .                  + thet * xmu * rn(k+1) * rmum1(k+1) * fr   (k+1)
      e00(ie, jr, k) = e00(ie, jr, k) 
     .                  - thet * xmu * rn(k  ) * rmum1(k  ) * fr   (k  )
c
      ep1(ie, jf, k) = ep1(ie, jf, k)   + thet * rmu  (k+1) * frnom(k+1)
      e00(ie, jf, k) = e00(ie, jf, k)   - thet * rmu  (k  ) * frnom(k  )
c
c-----------------------------------------------------------------------
c     radiation work
c-----------------------------------------------------------------------
c
c     equations re11 - re15
c
      ep1(ie, jr, k) = ep1(ie, jr, k) + xmu * thet
     .           * rn(k+1) * fedd(k) * er (k) *  rmum1(k+1) * u   (k+1)
      e00(ie, jr, k) = e00(ie, jr, k) - xmu * thet
     .           * rn(k  ) * fedd(k) * er (k) *  rmum1(k  ) * u   (k  )
c
      ep1(ie, ju, k) = ep1(ie, ju, k)
     .              + thet * fedd(k) * er (k) *  rmu  (k+1) * unom(k+1)
      e00(ie, ju, k) = e00(ie, ju, k) 
     .              - thet * fedd(k) * er (k) *  rmu  (k  ) * unom(k  )
c
      e00(ie, je, k) = e00(ie, je, k)
     .              + thet * fedd(k) * ern(k) * (rmu  (k+1) * u   (k+1) 
     .                                         - rmu  (k  ) * u   (k  ))
c
c-----------------------------------------------------------------------
c     radiation anisotropy
c-----------------------------------------------------------------------
c
c     equations re24 - re28
c
      ep1(ie, jr, k) = ep1(ie, jr, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     . *( + thet * rn(k+1) * aur(k) * rmu(k+1) +  dardlrp1(k) * dvol(k))
      e00(ie, jr, k) = e00(ie, jr, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     . *( - thet * rn(k  ) * aur(k) * rmu(k  ) +  dardlr00(k) * dvol(k))
c
      ep1(ie, ju, k) = ep1(ie, ju, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     .                                         *  dardlup1(k) * dvol(k) 
      e00(ie, ju, k) = e00(ie, ju, k) 
     .                  + 0.5D0 * xmu * er (k) * (1.0D0 - 3.0D0*fedd(k))
     .                                         *  dardlu00(k) * dvol(k) 
c
      e00(ie, je, k) = e00(ie, je, k) 
     .           + thet * 0.5D0 * xmu * ern(k) * (1.0D0 - 3.0D0*fedd(k))
     .                                * aur(k) *                dvol(k)
c
c-----------------------------------------------------------------------
c     source - sink terms
c-----------------------------------------------------------------------
c
c     equations re29 - re33
c
      ep1(ie, jr, k) = ep1(ie, jr, k) 
     .       + (          cc * xke (k) * er (k) 
     .         - 4.0D0 * cpi * xkp (k) * plf(k) ) 
     .                                * thet * rn(k+1) * d(k) * rmu(k+1)
      e00(ie, jr, k) = e00(ie, jr, k) 
     .       - (          cc * xke (k) * er (k) 
     .         - 4.0D0 * cpi * xkp (k) * plf(k) ) 
     .                                * thet * rn(k  ) * d(k) * rmu(k  )
c
      e00(ie, jd, k) = e00(ie, jd, k)  
     .       + (          cc * xke (k) * er (k) 
     .         - 4.0D0 * cpi * xkp (k) * plf(k) ) 
     .                                          * thet * dn(k) * dvol(k)
     .       + (          cc * xken(k) * er (k) *  dlkedldn(k)
     .         - 4.0D0 * cpi * xkpn(k) * plf(k) *  dlkpdldn(k) )
     .                                          * thet * d (k) * dvol(k)
c
      e00(ie, jt, k) = e00(ie, jt, k)  
     .       + (          cc * xken(k) * er (k) *  dlkedltn(k)
     .         - 4.0D0 * cpi * xkpn(k) * plf(k) *  dlkpdltn(k) 
     .         - 4.0D0 * cpi * xkp (k) *           dplfdltn(k) )
     .                                          * thet * d (k) * dvol(k)
c
      e00(ie, je, k) = e00(ie, je, k)
     .                  + cc * xke (k) * ern(k) * thet * d (k) * dvol(k) 
    7 continue
c
c=======================================================================
c     advection and right-hand side
c=======================================================================
c
c-----------------------------------------------------------------------
c     interior zones
c-----------------------------------------------------------------------
c
      do 8 k = ngrs + 1, ngre - 1
c
c.......................................................................
c     advection
c.......................................................................
c
c     equations re38 & re39
c
      e00(ie, jm, k) = e00(ie, jm, k)   - erb(k  ) * xmen(k  ) / dtime
      ep1(ie, jm, k) = ep1(ie, jm, k)   + erb(k+1) * xmen(k+1) / dtime
c
c     equations re40 - re49
c
      ep2(ie, jd, k) = ep2(ie, jd, k)        + dmdt(k+1) * dqbdlqp1(k+1)
      ep1(ie, jd, k) = ep1(ie, jd, k)        - dmdt(k  ) * dqbdlqp1(k  )
     .                                       + dmdt(k+1) * dqbdlq00(k+1)
      e00(ie, jd, k) = e00(ie, jd, k)        - dmdt(k  ) * dqbdlq00(k  )
     .                                       + dmdt(k+1) * dqbdlqm1(k+1)
      em1(ie, jd, k) = em1(ie, jd, k)        - dmdt(k  ) * dqbdlqm1(k  )
     .                                       + dmdt(k+1) * dqbdlqm2(k+1)
      em2(ie, jd, k) = em2(ie, jd, k)        - dmdt(k  ) * dqbdlqm2(k  )
c
      ep2(ie, je, k) = ep2(ie, je, k)        - dmdt(k+1) * dqbdlqp1(k+1)
      ep1(ie, je, k) = ep1(ie, je, k)        + dmdt(k  ) * dqbdlqp1(k  )
     .                                       - dmdt(k+1) * dqbdlq00(k+1)
      e00(ie, je, k) = e00(ie, je, k)        + dmdt(k  ) * dqbdlq00(k  )
     .                                       - dmdt(k+1) * dqbdlqm1(k+1)
      em1(ie, je, k) = em1(ie, je, k)        + dmdt(k  ) * dqbdlqm1(k  )
     .                                       - dmdt(k+1) * dqbdlqm2(k+1)
      em2(ie, je, k) = em2(ie, je, k)        + dmdt(k  ) * dqbdlqm2(k  )
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation re50
c
      rhs(ie, k) =                        ( ern(k) * dvoln(k) 
     .                                    - ero(k) * dvolo(k) ) / dtime 
     .
     .                                          - dmdt(k+1) * erb (k+1)
     .                                          + dmdt(k  ) * erb (k  )
     .
     .                                           + rmu(k+1) * fr  (k+1)
     .                                           - rmu(k  ) * fr  (k  )
     .
     .                      + fedd(k) * er (k) * ( rmu(k+1) *  u  (k+1) 
     .                                           - rmu(k  ) *  u  (k  ))
     .
     .            + ( 1.0D0 - fedd(k) * 3.0D0 )
     .                  * 0.5D0 * xmu * er (k)   * aur(k)   * dvol(k  )
     .
     .       + (          cc * xke(k) * er (k) 
     .         - 4.0D0 * cpi * xkp(k) * plf(k) ) * d  (k)   * dvol(k  )
    8 continue
c
c-----------------------------------------------------------------------
c     inner boundary condition (all cases); equation bc18
c-----------------------------------------------------------------------
c
      k = ngrs
c
c.......................................................................
c     advection; equations re38 & re41 - re44 & re46 - re49
c.......................................................................
c
      e00(ie, jr, k) = e00(ie, jr, k) 
     .                - thet * xmu * rmum1(k) *  rn(k) * urel(k) * er(k)
     .                             + rmu  (k) * (rn(k) / dtime)  * er(k)
c
      ep1(ie, jm, k) = ep1(ie, jm, k)    + erb(k+1) * xmen(k+1) / dtime
c
      ep2(ie, jd, k) = ep2(ie, jd, k)        + dmdt(k+1) * dqbdlqp1(k+1)
      ep1(ie, jd, k) = ep1(ie, jd, k)        + dmdt(k+1) * dqbdlq00(k+1)
      e00(ie, jd, k) = e00(ie, jd, k)        + dmdt(k+1) * dqbdlqm1(k+1)
      em1(ie, jd, k) = em1(ie, jd, k)        + dmdt(k+1) * dqbdlqm2(k+1)
c
      e00(ie, ju, k) = e00(ie, ju, k) - thet * rmu(k) * unom(k) * er (k)
c
      e00(ie, je, k) = e00(ie, je, k) - thet * rmu(k) * urel(k) * ern(k)
c
      ep2(ie, je, k) = ep2(ie, je, k)        - dmdt(k+1) * dqbdlqp1(k+1)
      ep1(ie, je, k) = ep1(ie, je, k)        - dmdt(k+1) * dqbdlq00(k+1)
      e00(ie, je, k) = e00(ie, je, k)        - dmdt(k+1) * dqbdlqm1(k+1)
      em1(ie, je, k) = em1(ie, je, k)        - dmdt(k+1) * dqbdlqm2(k+1)
c
c.......................................................................
c     right-hand side; equation bc18
c.......................................................................
c
      rhs(ie, k) =                        ( ern(k) * dvoln(k)  
     .                                    - ero(k) * dvolo(k) ) / dtime 
     .
     .                                          - dmdt(k+1) * erb (k+1)
     .                                 - urel(k) * rmu(k  ) * er  (k  )
     .
     .                                           + rmu(k+1) * fr  (k+1)
     .                                           - rmu(k  ) * fr  (k  )
     .
     .                      + fedd(k) * er (k) * ( rmu(k+1) * u   (k+1) 
     .                                           - rmu(k  ) * u   (k  ))
     .
     .            + ( 1.0D0 - fedd(k) * 3.0D0 )
     .                  * 0.5D0 * xmu * er (k)   * aur(k)   * dvol(k  )
     .
     .       + (          cc * xke(k) * er (k) 
     .         - 4.0D0 * cpi * xkp(k) * plf(k) ) * d  (k)   * dvol(k  )
c
c-----------------------------------------------------------------------
c     outer boundary condition (all cases); equation bc34
c-----------------------------------------------------------------------
c
      k = ngre
c
c.......................................................................
c     advection; equations re38 & re40 - re43 & re45 - re48
c.......................................................................
c
      e00(ie, jm, k) = e00(ie, jm, k)        - erb(k) * xmen(k) / dtime
c
      ep1(ie, jd, k) = ep1(ie, jd, k)            - dmdt(k) * dqbdlqp1(k)
      e00(ie, jd, k) = e00(ie, jd, k)            - dmdt(k) * dqbdlq00(k)
      em1(ie, jd, k) = em1(ie, jd, k)            - dmdt(k) * dqbdlqm1(k)
      em2(ie, jd, k) = em2(ie, jd, k)            - dmdt(k) * dqbdlqm2(k)
c
      ep1(ie, je, k) = ep1(ie, je, k)            + dmdt(k) * dqbdlqp1(k)
      e00(ie, je, k) = e00(ie, je, k)            + dmdt(k) * dqbdlq00(k)
      em1(ie, je, k) = em1(ie, je, k)            + dmdt(k) * dqbdlqm1(k)
      em2(ie, je, k) = em2(ie, je, k)            + dmdt(k) * dqbdlqm2(k)
c
c     extra terms for lagrangean or nontransmitting eulerian boundary;
c     equation bc34
c
      if (leobc .ne. 3) then
c
            ep1(ie, jr, k) = ep1(ie, jr, k)
     .         + thet * xmu * rmum1(k+1) *  rn(k+1) * urel(k+1) * er (k)
     .                      - rmu  (k+1) * (rn(k+1) / dtime)    * er (k)
c
            ep1(ie, ju, k) = ep1(ie, ju, k) 
     .                            + thet * rmu(k+1) * unom(k+1) * er (k)
c
            e00(ie, je, k) = e00(ie, je, k) 
     .                            + thet * rmu(k+1) * urel(k+1) * ern(k) 
      end if
c
c     extra terms for transmitting eulerian boundary; equation bc58
c
      if (leobc .eq. 3) then
c
            ep1(ie, jm, k) = ep1(ie, jm, k) 
     .                                  + erb(k+1) * xmen(k+1) / dtime
c
            ep2(ie, jd, k) = ep2(ie, jd, k) + dmdt(k+1) * dqbdlqp1(k+1)
            ep1(ie, jd, k) = ep1(ie, jd, k) + dmdt(k+1) * dqbdlq00(k+1)
            e00(ie, jd, k) = e00(ie, jd, k) + dmdt(k+1) * dqbdlqm1(k+1)
            em1(ie, jd, k) = em1(ie, jd, k) + dmdt(k+1) * dqbdlqm2(k+1)
c
            ep2(ie, je, k) = ep2(ie, je, k) - dmdt(k+1) * dqbdlqp1(k+1)
            ep1(ie, je, k) = ep1(ie, je, k) - dmdt(k+1) * dqbdlq00(k+1)
            e00(ie, je, k) = e00(ie, je, k) - dmdt(k+1) * dqbdlqm1(k+1)
            em1(ie, je, k) = em1(ie, je, k) - dmdt(k+1) * dqbdlqm2(k+1)
      end if
c
c.......................................................................
c     right-hand side; equation bc34
c.......................................................................
c
      rhs(ie, k) =                         ( ern(k) * dvoln(k) 
     .                                     - ero(k) * dvolo(k) ) / dtime
     .
     .                                          + dmdt(k  ) * erb (k  )
     .
     .                                           + rmu(k+1) *  fr (k+1)
     .                                           - rmu(k  ) *  fr (k  )
     .
     .                      + fedd(k) * er (k) * ( rmu(k+1) *   u (k+1) 
     .                                           - rmu(k  ) *   u (k  ))
     .
     .            + ( 1.0D0 - fedd(k) * 3.0D0 )
     .                  * 0.5D0 * xmu * er (k)   * aur(k)   * dvol(k  )
     .
     .       + (          cc * xke(k) * er (k) 
     .         - 4.0D0 * cpi * xkp(k) * plf(k) ) *   d(k)   * dvol(k  )
c
c     extra terms for lagrangean or nontransmitting eulerian boundary;
c     equation bc34
c
      if (leobc .ne. 3) rhs(ie, k) = rhs(ie, k) 
     .                                    + urel(k+1) * rmu(k+1) * er(k)
c
c     extra terms for transmitting eulerian boundary; equation bc58
c
      if (leobc .eq. 3) rhs(ie, k) = rhs(ie, k) - dmdt(k+1) * erb(k+1)
c
      end if
c
c++++++++ ltran .lt. 3 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 9 k = 1, ngrs - 2
      rhs(ie,     k) = 0.0D0
      e00(ie, je, k) = 1.0D0
    9 continue
c
      do 10 k = ngre + 3, mgr
      rhs(ie,     k) = 0.0D0
      e00(ie, je, k) = 1.0D0
   10 continue
c
c-----------------------------------------------------------------------
c     inner boundary; equations pz31 & pz32
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
      rhs(ie,     k) = + ern(k+1) - ern(k)
      ep1(ie, je, k) = + ern(k+1)
      e00(ie, je, k) =            - ern(k)
c
c-----------------------------------------------------------------------
c     outer boundary
c-----------------------------------------------------------------------
c
      k = ngre + 1
c
c     optically transmitting; equations pz79 & pz80
c
      if (lrobc .eq. 1) then
            rhs(ie,     k  ) = - ern(k-1) + ern(k  )
            e00(ie, je, k  ) =              ern(k  )
            em1(ie, je, k  ) = - ern(k-1)
            rhs(ie,     k+1) = - ern(k  ) + ern(k+1)
            e00(ie, je, k+1) =              ern(k+1)
            em1(ie, je, k+1) = - ern(k  )
      end if
c
c     optically reflecting; equations pz83 & pz84
c
      if (lrobc .eq. 2) then
            rhs(ie,     k  ) = - ern(k-1) + ern(k  )
            e00(ie, je, k  ) =              ern(k  )
            em1(ie, je, k  ) = - ern(k-1)
            rhs(ie,     k+1) = - ern(k-2) + ern(k+1)
            e00(ie, je, k+1) =              ern(k+1)
            em1(ie, je, k+1) = - ern(k-2)
      end if
c
c     imposed flux; equations pz87 & pz88
c
      if (lrobc .eq. 3) then
             rhs(ie,     k  ) = - ern(k-1) + ern(k  )
             e00(ie, je, k  ) =              ern(k  )
             em1(ie, je, k  ) = - ern(k-1)
             rhs(ie,     k+1) = - ern(k  ) + ern(k+1)
             e00(ie, je, k+1) =              ern(k+1)
             em1(ie, je, k+1) = - ern(k  )
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
