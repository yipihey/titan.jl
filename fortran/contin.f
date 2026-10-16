      subroutine contin
c
c***********************************************************************
c     equation of continuity c2
c.......................................................................
c     calling sequence: matgen > contin > advectc
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq1'
      include'titan.eq3'
c
c=======================================================================
c     set interface densities for advection
c=======================================================================
c
      do 2 l = 1, madv
      do 1 k = 1, mgr
      adv(k, l) = 0.0D0
    1 continue
    2 continue
c
      do 3   k  = ngrs - 1, ngre + 2
      qn    (k) = dn(k)
      qo    (k) = do(k)
    3 continue
      do 33  k  = ngrs    , ngre + 1
      qso   (k) = dso(k)
   33 continue
      do 333 k  = ngrs + 1, ngre + 1
      flow  (k) = urel(k)
  333 continue
c
      call advectc
c
      do 4   k  = ngrs + 1, ngre + 1
      db    (k) = qb(k) 
    4 continue
      do 44  k  = ngrs    , ngre + 1
      dsn   (k) = qsn(k)
   44 continue
c
c=======================================================================
c     set mass diffusion flux
c=======================================================================
c
      do 6 l = 1, mdfz
      do 5 k = 1, mgr
      dfz(k, l) = 0.0D0
    5 continue
    6 continue
c
      sig = sigd
      zet = 0.0D0
c
      call diffuse
c
c=======================================================================
c     time derivative and diffusion at all points, including boundaries
c=======================================================================
c
      do 8 k = ngrs, ngre
c
c-----------------------------------------------------------------------
c     time derivative
c-----------------------------------------------------------------------
c
c     equations c13 - c15
c
      ep1(id, jr, k) = ep1(id, jr, k)      + dn(k) * rmup1n(k+1) / dtime
      e00(id, jr, k) = e00(id, jr, k)      - dn(k) * rmup1n(k  ) / dtime
c
      e00(id, jd, k) = e00(id, jd, k)      + dn(k) * dvoln (k  ) / dtime
c
c-----------------------------------------------------------------------
c     diffusion
c-----------------------------------------------------------------------
c
c     equations c50 - c56
c
      ep2(id, jr, k) = ep2(id, jr, k)                    - ddfdlrp1(k+1)
      ep1(id, jr, k) = ep1(id, jr, k)      + ddfdlrp1(k) - ddfdlr00(k+1)
      e00(id, jr, k) = e00(id, jr, k)      + ddfdlr00(k) - ddfdlrm1(k+1)
      em1(id, jr, k) = em1(id, jr, k)      + ddfdlrm1(k)
c
      ep1(id, jd, k) = ep1(id, jd, k)                    - ddfdld00(k+1)
      e00(id, jd, k) = e00(id, jd, k)      + ddfdld00(k) - ddfdldm1(k+1)
      em1(id, jd, k) = em1(id, jd, k)      + ddfdldm1(k)
    8 continue
c
c=======================================================================
c     advection and right-hand side
c=======================================================================
c
c-----------------------------------------------------------------------
c     interior zones
c-----------------------------------------------------------------------
c
      do 9 k = ngrs + 1, ngre - 1
c
c.......................................................................
c     advection
c.......................................................................
c
c     equations c16 & c17
c
      ep1(id, jr, k) = ep1(id, jr, k) 
     .         + thet * xmu * rn(k+1) * rmum1(k+1) * urel(k+1) * db(k+1)
     .                      - rn(k+1) / dtime      * rmu (k+1) * db(k+1)
      e00(id, jr, k) = e00(id, jr, k)
     .         - thet * xmu * rn(k  ) * rmum1(k  ) * urel(k  ) * db(k  )
     .                      + rn(k  ) / dtime      * rmu (k  ) * db(k  )
c
c     equations c18 & c19
c
      ep1(id, ju, k) = ep1(id, ju, k) 
     .                         + thet * rmu  (k+1) * unom(k+1) * db(k+1) 
      e00(id, ju, k) = e00(id, ju, k) 
     .                         - thet * rmu  (k  ) * unom(k  ) * db(k  ) 
c
c     equations c37 - c41
c
      ep2(id, jd, k) = ep2(id, jd, k)
     .                            + rmu(k+1) * urel(k+1) * dqbdlqp1(k+1)
      ep1(id, jd, k) = ep1(id, jd, k)
     .                            + rmu(k+1) * urel(k+1) * dqbdlq00(k+1)
     .                            - rmu(k  ) * urel(k  ) * dqbdlqp1(k  )
      e00(id, jd, k) = e00(id, jd, k) 
     .                            + rmu(k+1) * urel(k+1) * dqbdlqm1(k+1)
     .                            - rmu(k  ) * urel(k  ) * dqbdlq00(k  )
      em1(id, jd, k) = em1(id, jd, k)
     .                            + rmu(k+1) * urel(k+1) * dqbdlqm2(k+1)
     .                            - rmu(k  ) * urel(k  ) * dqbdlqm1(k  )
      em2(id, jd, k) = em2(id, jd, k)
     .                            - rmu(k  ) * urel(k  ) * dqbdlqm2(k  )
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation c57
c
      rhs(id, k) =                            ( dn(k)*dvoln(k) 
     .                                        - do(k)*dvolo(k) ) / dtime
     .          
     .                                  + rmu(k+1) * urel(k+1) * db(k+1)
     .                                  - rmu(k  ) * urel(k  ) * db(k  ) 
     .          
     .                                                         - df(k+1) 
     .                                                         + df(k  )
    9 continue
c
c-----------------------------------------------------------------------
c     inner boundary condition
c-----------------------------------------------------------------------
c
      k = ngrs
c
c.......................................................................
c     advection
c.......................................................................
c
c     equation c16 & c17
c
      ep1(id, jr, k) = ep1(id, jr, k)
     .         + thet * xmu * rn(k+1) * rmum1(k+1) * urel(k+1) * db(k+1)
     .                      - rn(k+1) / dtime      * rmu (k+1) * db(k+1)
      e00(id, jr, k) = e00(id, jr, k) 
     .         - thet * xmu * rn(k  ) * rmum1(k  ) * phil0
c
c     equations c19
c
      ep1(id, ju, k) = ep1(id, ju, k) 
     .                         + thet * rmu  (k+1) * unom(k+1) * db(k+1) 
c
c     equations c38 - c41
c
      ep2(id, jd, k) = ep2(id, jd, k)
     .                          + rmu  (k+1) * urel(k+1) * dqbdlqp1(k+1)
      ep1(id, jd, k) = ep1(id, jd, k)
     .                          + rmu  (k+1) * urel(k+1) * dqbdlq00(k+1)
      e00(id, jd, k) = e00(id, jd, k) 
     .                          + rmu  (k+1) * urel(k+1) * dqbdlqm1(k+1)
      em1(id, jd, k) = em1(id, jd, k)
     .                          + rmu  (k+1) * urel(k+1) * dqbdlqm2(k+1)
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation bc15
c
      rhs(id, k) =                            ( dn(k)*dvoln(k) 
     .                                        - do(k)*dvolo(k) ) / dtime
     .
     .                                  + rmu(k+1) * urel(k+1) * db(k+1)
     .                                  - rmu(k  ) * phil0
     .
     .                                                         - df(k+1)
     .                                                         + df(k  )
c
c-----------------------------------------------------------------------
c     outer boundary condition
c-----------------------------------------------------------------------
c
      k = ngre
c
c.......................................................................
c     advection
c.......................................................................
c
c     equation c16 & c17
c
      ep1(id, jr, k) = ep1(id, jr, k) 
     .         + thet * xmu * rn(k+1) * rmum1(k+1) * phir0
      e00(id, jr, k) = e00(id, jr, k)
     .         - thet * xmu * rn(k  ) * rmum1(k  ) * urel(k  ) * db(k  )
     .                      + rn(k  ) / dtime      * rmu (k  ) * db(k  )
c
c     equation c18
c
      e00(id, ju, k) = e00(id, ju, k)   
     .                         - thet * rmu  (k  ) * unom(k  ) * db(k  )
c
c     equation c37 - c40
c
      ep1(id, jd, k) = ep1(id, jd, k)   - rmu(k) * urel(k) * dqbdlqp1(k)
      e00(id, jd, k) = e00(id, jd, k)   - rmu(k) * urel(k) * dqbdlq00(k)
      em1(id, jd, k) = em1(id, jd, k)   - rmu(k) * urel(k) * dqbdlqm1(k)
      em2(id, jd, k) = em2(id, jd, k)   - rmu(k) * urel(k) * dqbdlqm2(k)
c
c     extra terms for eulerian transmitting boundary;
c     equation c16, c19, c38 - c41
c
      if( leobc .eq. 3 ) then
c
            ep1(id, jr, k) = ep1(id, jr, k) 
     .         + thet * xmu * rn(k+1) * rmum1(k+1) * urel(k+1) * db(k+1)
     .                      - rn(k+1) / dtime      * rmu (k+1) * db(k+1)
c
            ep1(id, ju, k) = ep1(id, ju, k) 
     .                         + thet * rmu  (k+1) * unom(k+1) * db(k+1)
c
            ep2(id, jd, k) = ep2(id, jd, k)
     .                          + rmu  (k+1) * urel(k+1) * dqbdlqp1(k+1)
            ep1(id, jd, k) = ep1(id, jd, k)
     .                          + rmu  (k+1) * urel(k+1) * dqbdlq00(k+1)
            e00(id, jd, k) = e00(id, jd, k) 
     .                          + rmu  (k+1) * urel(k+1) * dqbdlqm1(k+1)
            em1(id, jd, k) = em1(id, jd, k)
     .                          + rmu  (k+1) * urel(k+1) * dqbdlqm2(k+1)
c
      end if
c
c.......................................................................
c     right-hand side
c.......................................................................
c
c     equation bc31
c
      rhs(id, k) =                            ( dn(k)*dvoln(k) 
     .                                        - do(k)*dvolo(k) ) / dtime
     .
     .                                  + rmu(k+1) * phir0
     .                                  - rmu(k  ) * urel(k  ) * db(k  )
     .
     .                                                         - df(k+1)
     .                                                         + df(k  )
c
c     extra term for eulerian transmitting outer boundary; equation bc38
c
      if( leobc .eq. 3 ) rhs(id, k) = rhs(id, k)     
     .                                  + rmu(k+1) * urel(k+1) * db(k+1)
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 11 k = 1, ngrs - 2
      rhs(id,     k) = 0.0D0
      e00(id, jd, k) = 1.0D0
   11 continue
c
      do 12 k = ngre + 2, mgr
      rhs(id,     k) = 0.0D0
      e00(id, jd, k) = 1.0D0
   12 continue
c
c-----------------------------------------------------------------------
c     eulerian inner boundary
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
c     zero flux; equations pz5 & pz6
c
      if(leibc .eq. 1) then
            rhs(id,     k) =   dn(k+1) - dn(k)
            ep1(id, jd, k) =   dn(k+1)
            e00(id, jd, k) =           - dn(k)
      end if
c
c     nonzero flux; equations pz15 & pz16
c
      if(leibc .eq. 2) then
            rhs(id,     k) = 0.0D0
            e00(id, jd, k) = 1.0D0
      end if
c
c-----------------------------------------------------------------------
c     lagrangean inner boundary
c-----------------------------------------------------------------------
c
c     equations pz25 & pz26
c
      if( llibc .gt. 0 ) then
            rhs(id,     k) =   dn(k+1) - dn(k)
            ep1(id, jd, k) =   dn(k+1)
            e00(id, jd, k) =           - dn(k)
      end if
c
c-----------------------------------------------------------------------
c     eulerian outer boundary 
c-----------------------------------------------------------------------
c
c     zero flux; equations pz43 & pz44
c
      if(leobc .eq. 1) then
            do 13 k = ngre + 1, ngre + 2
            rhs(id,     k) = - dn(k-1) + dn(k)
            e00(id, jd, k) =             dn(k)
            em1(id, jd, k) = - dn(k-1)
   13       continue
      end if
c
c     nonzero flux; equations pz53 & pz54
c
      if(leobc .eq. 2) then
            do 14 k = ngre + 1, ngre + 2
            rhs(id,     k) = 0.0D0
            e00(id, jd, k) = 1.0D0
   14       continue
      end if
c
c     transmitting; equations pz63 & pz64
c
      if( leobc .eq. 3 ) then
            do 15 k = ngre + 1, ngre + 2
            rhs(id,     k) = - dn(k-1) + dn(k)
            e00(id, jd, k) =             dn(k)
            em1(id, jd, k) = - dn(k-1)
   15       continue
      end if
c
c-----------------------------------------------------------------------
c     lagrangean outer boundary 
c-----------------------------------------------------------------------
c
c     equations pz73 & pz74
c
      if( llobc .gt. 0 ) then
            do 16 k = ngre + 1, ngre + 2
            rhs(id,     k) = - dn(k-1) + dn(k)
            e00(id, jd, k) =             dn(k)
            em1(id, jd, k) = - dn(k-1)
   16       continue
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
