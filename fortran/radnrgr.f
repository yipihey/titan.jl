      subroutine radnrgr
c
c***********************************************************************
c     radiation energy equation re2
c     static material, time-dependent radiation.
c.......................................................................
c     calling sequence: matgen > radnrgr > advectc
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq1'
c
c***********************************************************************
c     set auxiliary quantities
c***********************************************************************
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
c++++++++ ltran .eq. 3 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (ltran .eq. 3) then
c
c     equations rd44 & rd45
c
      do 3 k = ngrs, ngre
      e00(ie, jt, k) = -4.0D0 * car * tn(k)**4
      e00(ie, je, k) =                           ern(k)
c
c     equation rd46
c
      rhs(ie,     k) =        - car * tn(k)**4 + ern(k)
    3 continue
c
      end if
c
c++++++++ ltran .eq. 3 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c***********************************************************************
c     radiation transport or nonequilibrium diffusion; equation re2
c     (no velocity terms)
c***********************************************************************
c
c++++++++ ltran .lt. 3 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (ltran .lt. 3) then
c
c=======================================================================
c     set interface energy densities for advection
c=======================================================================
c
      do 4   k  = ngrs - 1, ngre + 2
      qn    (k) = ern (k) / dn(k)
      qo    (k) = ero (k) / do(k)
    4 continue
      do 44  k  = ngrs    , ngre + 1
      qso   (k) = erso(k)
   44 continue
      do 444 k  = ngrs + 1, ngre + 1
      flow  (k) = - drdt(k)*rmu(k)
  444 continue
c
      call advectc
c
      do 5  k  = ngrs + 1, ngre + 1
      erb  (k) = qb (k)
    5 continue
      do 55 k  = ngrs    , ngre + 1
      ersn (k) = qsn(k)
   55 continue
c
c=======================================================================
c     all terms except advection and rhs; note that this loop covers
c     the entire domain including boundaries
c=======================================================================
c
      do 6 k = ngrs, ngre
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
     .                    + thet * xmu * rn(k+1) * rmum1(k+1) * fr (k+1)
      e00(ie, jr, k) = e00(ie, jr, k) 
     .                    - thet * xmu * rn(k  ) * rmum1(k  ) * fr (k  )
c
      ep1(ie, jf, k) = ep1(ie, jf, k)   + thet * rmu  (k+1) * frnom(k+1)
      e00(ie, jf, k) = e00(ie, jf, k)   - thet * rmu  (k  ) * frnom(k  )
c
c-----------------------------------------------------------------------
c     source - sink terms
c-----------------------------------------------------------------------
c
c     equations re29 - re30 & re32 - re33
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
      e00(ie, jt, k) = e00(ie, jt, k)  
     .       + (          cc * xken(k) * er (k) *  dlkedltn(k)
     .         - 4.0D0 * cpi * xkpn(k) * plf(k) *  dlkpdltn(k) 
     .         - 4.0D0 * cpi * xkp (k) *           dplfdltn(k) )
     .                                          * thet * d (k) * dvol(k)
c
      e00(ie, je, k) = e00(ie, je, k)
     .                  + cc * xke (k) * ern(k) * thet * d (k) * dvol(k) 
    6 continue
c
c=======================================================================
c     advection and right-hand side
c=======================================================================
c
c-----------------------------------------------------------------------
c     interior zones
c-----------------------------------------------------------------------
c
      do 7 k = ngrs + 1, ngre - 1
c
c.......................................................................
c     advection
c.......................................................................
c
      ep1(ie, jr, k) = ep1(ie, jr, k) 
     .- thet*xmu * dn(k+1) * rn(k+1) * rmum1(k+1) * drdt(k+1) * erb(k+1)
     .           - dn(k+1) * rn(k+1) / dtime      * rmu (k+1) * erb(k+1)
      e00(ie, jr, k) = e00(ie, jr, k)
     .+ thet*xmu * dn(k  ) * rn(k  ) * rmum1(k  ) * drdt(k  ) * erb(k  )
     .           + dn(k  ) * rn(k  ) / dtime      * rmu (k  ) * erb(k  )
c
      ep2(ie, je, k) = ep2(ie, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlqp1(k+1)
      ep1(ie, je, k) = ep1(ie, je, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlqp1(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlq00(k+1)
      e00(ie, je, k) = e00(ie, je, k) 
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlq00(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlqm1(k+1)
      em1(ie, je, k) = em1(ie, je, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlqm1(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlqm2(k+1)
      em2(ie, je, k) = em2(ie, je, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlqm2(k  )
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation re50 without velocity terms
c
      rhs(ie, k) =                        ( ern(k) * dvoln(k) 
     .                                    - ero(k) * dvolo(k) ) / dtime 
     .
     .                      - dn(k+1) * drdt(k+1) * rmu(k+1) * erb(k+1)
     .                      + dn(k  ) * drdt(k  ) * rmu(k  ) * erb(k  ) 
     .
     .                                            + rmu(k+1) * fr (k+1)
     .                                            - rmu(k  ) * fr (k  )
     .
     .           + (          cc * xke(k) * er (k) 
     .             - 4.0D0 * cpi * xkp(k) * plf(k) ) * d(k) * dvol(k  )
    7 continue
c
c-----------------------------------------------------------------------
c     inner boundary condition (all cases); equation bc18
c-----------------------------------------------------------------------
c
      k = ngrs
c
c.......................................................................
c     advection
c.......................................................................
c
      ep1(ie, jr, k) = ep1(ie, jr, k) 
     .- thet*xmu * dn(k+1) * rn(k+1) * rmum1(k+1) * drdt(k+1) * erb(k+1)
     .           - dn(k+1) * rn(k+1) / dtime      * rmu (k+1) * erb(k+1)
c
      ep2(ie, je, k) = ep2(ie, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlqp1(k+1)
      ep1(ie, je, k) = ep1(ie, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlq00(k+1)
      e00(ie, je, k) = e00(ie, je, k) 
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlqm1(k+1)
      em1(ie, je, k) = em1(ie, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlqm2(k+1)
c
c.......................................................................
c     right-hand side; equation bc18
c.......................................................................
c
      rhs(ie, k) =                        ( ern(k) * dvoln(k)  
     .                                    - ero(k) * dvolo(k) ) / dtime 
     .
     .                      - dn(k+1) * drdt(k+1) * rmu(k+1) * erb(k+1)
     .
     .                                            + rmu(k+1) * fr (k+1)
     .                                            - rmu(k  ) * fr (k  )
     .
     .              + (          cc * xke(k) * er (k) 
     .                - 4.0D0 * cpi * xkp(k) * plf(k) ) * d(k) * dvol(k)
c
c-----------------------------------------------------------------------
c     outer boundary condition (all cases); equation bc34
c-----------------------------------------------------------------------
c
      k = ngre
c
c.......................................................................
c     advection
c.......................................................................
c
      e00(ie, jr, k) = e00(ie, jr, k) 
     .        + thet * xmu * dn(k) * rn(k) * rmum1(k) * drdt(k) * erb(k)
     .                     + dn(k) * rn(k) / dtime    * rmu (k) * erb(k)
c
      ep1(ie, je, k) = ep1(ie, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlqp1(k)
      e00(ie, je, k) = e00(ie, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlq00(k)
      em1(ie, je, k) = em1(ie, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlqm1(k)
      em2(ie, je, k) = em2(ie, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlqm2(k)
c
c.......................................................................
c     right-hand side, equation bc34
c.......................................................................
c
      rhs(ie, k) =                        ( ern(k) * dvoln(k)  
     .                                    - ero(k) * dvolo(k) ) / dtime 
     .
     .                        + dn(k) * drdt(k  ) * rmu(k  ) * erb(k  )
     .
     .                                            + rmu(k+1) * fr (k+1)
     .                                            - rmu(k  ) * fr (k  )
     .
     .              + (          cc * xke(k) * er (k) 
     .                - 4.0D0 * cpi * xkp(k) * plf(k) ) * d(k) * dvol(k)
c
      end if
c
c++++++++ ltran .lt. 3 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 8 k = 1, ngrs - 2
      rhs(ie,     k) = 0.0D0
      e00(ie, je, k) = 1.0D0
    8 continue
c
      do 9 k = ngre + 3, mgr
      rhs(ie,     k) = 0.0D0
      e00(ie, je, k) = 1.0D0
    9 continue
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
