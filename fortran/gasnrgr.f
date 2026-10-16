      subroutine gasnrgr
c
c***********************************************************************
c     radiating fluid energy equation, equation fe3 (no velocity terms) 
c     static material, time-dependent radiation.
c.......................................................................
c     calling sequence: matgen > gasnrgr > advectc
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq1'
      include'titan.eq3'
c
c=======================================================================
c     set interface energy densities for advection
c=======================================================================
c
      do 2 l = 1, madv
      do 1 k = 1, mgr
      adv(k, l) = 0.0D0
    1 continue
    2 continue
c
      do 3   k  = ngrs - 1, ngre + 2
      qn    (k) = egrn (k)
      qo    (k) = egro (k)
    3 continue
      do 4   k  = ngrs    , ngre + 1
      qso   (k) = egrso(k)
    4 continue
      do 44  k  = ngrs + 1, ngre + 1
      flow  (k) = - drdt(k)*rmu(k)
   44 continue
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
c     equations fe9 & fe10
c
      dlqdlt(k) = ( egn(k) * dlegdltn(k)                  ) / egrn(k)
      dlqdle(k) = (                        ern(k) / dn(k) ) / egrn(k)
    6 continue
c
      do 7 k = ngrs + 1, ngre + 1
c
c     equations fe12 - fe19
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
    7 continue
c
c=======================================================================
c     all terms except advection and rhs, note that this loop covers
c     the entire domain, including boundaries
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
      e00(it, jt, k) = e00(it, jt, k)
     .                + dlegdltn(k)  * dn(k) * egn(k) * dvoln(k) / dtime
c
      e00(it, je, k) = e00(it, je, k) +        ern(k) * dvoln(k) / dtime
c
c-----------------------------------------------------------------------
c     flux divergence
c-----------------------------------------------------------------------
c
c     equations re6 - re7
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .                     + thet * xmu * rn(k+1) * rmum1(k+1) * fr(k+1)
      e00(it, jr, k) = e00(it, jr, k) 
     .                     - thet * xmu * rn(k  ) * rmum1(k  ) * fr(k  )
c
c     equations re8 - re9
c
      ep1(it, jf, k) = ep1(it, jf, k)     + thet * frnom(k+1) * rmu(k+1)
      e00(it, jf, k) = e00(it, jf, k)     - thet * frnom(k  ) * rmu(k  )
c
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
      ep1(it, jr, k) = ep1(it, jr, k) 
     .- thet*xmu* dn(k+1) * rn(k+1) * rmum1(k+1) * drdt(k+1) * egrb(k+1)
     .          - dn(k+1) * rn(k+1) / dtime      * rmu (k+1) * egrb(k+1)
      e00(it, jr, k) = e00(it, jr, k)
     .+ thet*xmu* dn(k  ) * rn(k  ) * rmum1(k  ) * drdt(k  ) * egrb(k  )
     .          + dn(k  ) * rn(k  ) / dtime      * rmu (k  ) * egrb(k  )
c
      ep2(it, jt, k) = ep2(it, jt, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdltp1(k+1)
      ep1(it, jt, k) = ep1(it, jt, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdltp1(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlt00(k+1)
      e00(it, jt, k) = e00(it, jt, k) 
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlt00(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdltm1(k+1)
      em1(it, jt, k) = em1(it, jt, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdltm1(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdltm2(k+1)
      em2(it, jt, k) = em2(it, jt, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdltm2(k  )
c
      ep2(it, je, k) = ep2(it, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlep1(k+1)
      ep1(it, je, k) = ep1(it, je, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlep1(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdle00(k+1)
      e00(it, je, k) = e00(it, je, k) 
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdle00(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlem1(k+1)
      em1(it, je, k) = em1(it, je, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlem1(k  )
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlem2(k+1)
      em2(it, je, k) = em2(it, je, k)
     .                  + dn(k  ) * rmu(k  ) * drdt(k  ) * dqbdlem2(k  )
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation fe97 (without velocity terms)
c
      rhs(it, k) =                  ( dn(k) * egrn(k) * dvoln(k) 
     .                              - do(k) * egro(k) * dvolo(k) )/dtime
     .
     .                     - dn(k+1) * drdt(k+1) * rmu(k+1) * egrb(k+1)
     .                     + dn(k  ) * drdt(k  ) * rmu(k  ) * egrb(k  ) 
     .
     .                                           + rmu(k+1) *  fr (k+1)
     .                                           - rmu(k  ) *  fr (k  ) 
   11 continue
c
c-----------------------------------------------------------------------
c     inner boundary condition, all cases; equation bc17
c-----------------------------------------------------------------------
c
      k = ngrs
c
c.......................................................................
c     advection
c.......................................................................
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .- thet * xmu
     .         * dn(k+1) *  rn(k+1) * rmum1(k+1) * drdt(k+1) * egrb(k+1)
     .         - dn(k+1) * (rn(k+1) / dtime)     * rmu (k+1) * egrb(k+1)
c
      ep2(it, jt, k) = ep2(it, jt, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdltp1(k+1)
      ep1(it, jt, k) = ep1(it, jt, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlt00(k+1)
      e00(it, jt, k) = e00(it, jt, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdltm1(k+1)
      em1(it, jt, k) = em1(it, jt, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdltm2(k+1)
c
      ep2(it, je, k) = ep2(it, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlep1(k+1)
      ep1(it, je, k) = ep1(it, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdle00(k+1)
      e00(it, je, k) = e00(it, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlem1(k+1)
      em1(it, je, k) = em1(it, je, k)
     .                  - dn(k+1) * rmu(k+1) * drdt(k+1) * dqbdlem2(k+1)
c
      e00(it, jr, k) = e00(it, jr, k) 
     .                          - thet * rn (k) * xmu * rmum1(k) * phil2
c
c.......................................................................
c     right-hand side; equation bc17
c.......................................................................
c
      rhs(it, k) =                  ( dn(k) * egrn(k) * dvoln(k) 
     .                              - do(k) * egro(k) * dvolo(k) )/dtime
     .
     .                     - dn(k+1) * drdt(k+1) * rmu (k+1) * egrb(k+1)
     .                               - phil2     * rmu (k  )
     .
     .                                           + rmu (k+1) *  fr (k+1)
     .                                           - rmu (k  ) *  fr (k  )
c
c-----------------------------------------------------------------------
c     outer boundary condition, all cases; equation bc33
c-----------------------------------------------------------------------
c
      k = ngre
c
c.......................................................................
c     advection
c.......................................................................
c
      e00(it, jr, k) = e00(it, jr, k)
     .      + thet * xmu * dn(k) *  rn(k) * rmum1(k) * drdt(k) * egrb(k)
     .                   + dn(k) * (rn(k) / dtime)   * rmu (k) * egrb(k)
c.......................................................................
c
      ep1(it, jt, k) = ep1(it, jt, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdltp1(k)
      e00(it, jt, k) = e00(it, jt, k)
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlt00(k)
      em1(it, jt, k) = em1(it, jt, k)
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdltm1(k)
      em2(it, jt, k) = em2(it, jt, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdltm2(k)
c
      ep1(it, je, k) = ep1(it, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlep1(k)
      e00(it, je, k) = e00(it, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdle00(k)
      em1(it, je, k) = em1(it, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlem1(k)
      em2(it, je, k) = em2(it, je, k) 
     .                          + dn(k) * rmu(k) * drdt(k) * dqbdlem2(k)
c
      ep1(it, jr, k) = ep1(it, jr, k) 
     .                      + thet * rn (k+1) * xmu * rmum1(k+1) * phir2
c
c.......................................................................
c     right-hand side; equation bc33
c.......................................................................
c
      rhs(it, k) =                  ( dn(k) * egrn(k) * dvoln(k) 
     .                              - do(k) * egro(k) * dvolo(k) )/dtime
     .
     .                              + phir2     * rmu (k+1)
     .                      + dn(k) * drdt(k  ) * rmu (k  ) * egrb(k  ) 
     .
     .                                          + rmu (k+1) *  fr (k+1)
     .                                          - rmu (k  ) *  fr (k  ) 
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
c     zero flux, equations pz9 & pz10
c
      if (leibc .eq. 1) then
            rhs(it,     k) =   egn(k+1) - egn(k)
            ep1(it, jt, k) =   egn(k+1) * dlegdltn(k+1)
            e00(it, jt, k) = - egn(k  ) * dlegdltn(k  )
      end if
c
c     nonzero flux, equations pz19 & pz20
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
            ep1(it, jt, k) =   egn(k+1) * dlegdltn(k+1)
            e00(it, jt, k) = - egn(k  ) * dlegdltn(k  )
      end if
c
c-----------------------------------------------------------------------
c     eulerian outer boundary 
c-----------------------------------------------------------------------
c
c     zero flux, equations pz47 & pz48
c
      if (leobc .eq. 1) then
	    do 14 k = ngre + 1, ngre + 2
            rhs(it,     k) = - egn(k-1) + egn(k)
            e00(it, jt, k) =   egn(k  ) * dlegdltn(k  ) 
            em1(it, jt, k) = - egn(k-1) * dlegdltn(k-1)
   14       continue
      end if
c
c     nonzero flux, equations pz57 & pz58
c
      if (leobc .eq. 2) then
	    do 15 k = ngre + 1, ngre + 2
            rhs(it,     k) = 0.0D0
            e00(it, jt, k) = 1.0D0
   15       continue
      end if
c
c     transmitting, equations pz67 & pz68
c
      if (leobc .eq. 3) then
	    do 16 k = ngre + 1, ngre + 2
            rhs(it,     k) = - egn(k-1) + egn(k)
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
            e00(it, jt, k) =   egn(k  ) * dlegdltn(k  ) 
            em1(it, jt, k) = - egn(k-1) * dlegdltn(k-1)
   17       continue
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
