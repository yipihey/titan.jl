      subroutine gasmomh
c
c***********************************************************************
c     gas momentum equation gm4
c.......................................................................
c     calling sequence: matgen > gasmomh > advecti
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq2'
c
c=======================================================================
c     define various quantities and their derivatives
c=======================================================================
c
c-----------------------------------------------------------------------
c     f0 = <r>**(3 mu/2); equations gm39 - gm40
c-----------------------------------------------------------------------
c
      f0 (xk,xp) =             ( (xk+xp)/2. )**( 1.5D0*xmu )
      f0r(xk,xp) = .75D0*xmu * ( (xk+xp)/2. )**( 1.5D0*xmu - 1.0D0 )
c
c-----------------------------------------------------------------------
c     f1 = r**(-mu/2); equations gm41 - gm42
c-----------------------------------------------------------------------
c
      f1  (xk) =               xk**( - 0.5D0*xmu )
      f1r (xk) = - 0.5D0*xmu * xk**( - 0.5D0*xmu - 1.0D0 )
c
c=======================================================================
c     set auxiliary quantities for viscosity
c=======================================================================
c
      do 1 k = ngrs, ngre
c
c     equations gm43 - gm45
c
       r3     (k) = f0  (r(k),r(k+1))
      dr3dlrp1(k) = f0r (r(k),r(k+1)) * thet * rn(k+1)
      dr3dlr00(k) = f0r (r(k),r(k+1)) * thet * rn(k  )
c
c     equations gm46 & gm47
c
       r1     (k) = f1  (r(k))
      dr1dlr00(k) = f1r (r(k))        * thet * rn(k  )
    1 continue
c
c=======================================================================
c     set cell-interface quantities for advection
c=======================================================================
c
      do 3 l = 1, madv
      do 2 k = 1, mgr
      adv(k, l) = 0.0D0
    2 continue
    3 continue
c
      do 4  k  = ngrs - 1, ngre + 2
      qn   (k) = un (k)
      qo   (k) = uo (k)
    4 continue
      do 44 k  = ngrs, ngre
      qso  (k) = uso(k)
      flow (k) = - ( dmdt(k) + dmdt(k+1) )
   44 continue
c
      call advecti
c
      do 5 k  = ngrs, ngre
      ub  (k) = qb (k) 
      usn (k) = qsn(k)
    5 continue
c
c=======================================================================
c     interior zones
c=======================================================================
c
      do 6 k = ngrs + 1, ngre
c
c-----------------------------------------------------------------------
c     time derivative
c-----------------------------------------------------------------------
c
c     equations gm12 - gm17
c
      ep1(iu, jr, k) = ep1(iu, jr, k) 
     .               + un(k) *  dn(k  ) * rmup1n(k+1) / (2.0D0 * dtime)
      e00(iu, jr, k) = e00(iu, jr, k)
     .               - un(k) *  dn(k  ) * rmup1n(k  ) / (2.0D0 * dtime)
     .               + un(k) *  dn(k-1) * rmup1n(k  ) / (2.0D0 * dtime)
      em1(iu, jr, k) = em1(iu, jr, k) 
     .               - un(k) *  dn(k-1) * rmup1n(k-1) / (2.0D0 * dtime)
c
      e00(iu, jd, k) = e00(iu, jd, k)
     .               + un(k) *  dn(k  ) * dvoln (k  ) / (2.0D0 * dtime)
      em1(iu, jd, k) = em1(iu, jd, k)
     .               + un(k) *  dn(k-1) * dvoln (k-1) / (2.0D0 * dtime)
c
      e00(iu, ju, k) = e00(iu, ju, k)
     .             + unom(k) * (dn(k-1) * dvoln (k-1) 
     .                         +dn(k  ) * dvoln (k  ))/ (2.0D0 * dtime)
c
c-----------------------------------------------------------------------
c     advection
c-----------------------------------------------------------------------
c
c     equations gm20 - gm27
c
      ep1(iu, jm, k) = ep1(iu, jm, k) 
     .                             + 0.5D0 * xmen(k+1) * ub(k  ) / dtime
      e00(iu, jm, k) = e00(iu, jm, k) 
     .                             + 0.5D0 * xmen(k  ) * ub(k  ) / dtime
     .                             - 0.5D0 * xmen(k  ) * ub(k-1) / dtime
      em1(iu, jm, k) = em1(iu, jm, k) 
     .                             - 0.5D0 * xmen(k-1) * ub(k-1) / dtime
c
      ep2(iu, ju, k) = ep2(iu, ju, k) 
     .      - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdqp2(k  ) * unom(k+2)
      ep1(iu, ju, k) = ep1(iu, ju, k) 
     .      + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdqp2(k-1) * unom(k+1)
     .      - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdqp1(k  ) * unom(k+1)
      e00(iu, ju, k) = e00(iu, ju, k) 
     .      + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdqp1(k-1) * unom(k  )
     .      - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdq00(k  ) * unom(k  )
      em1(iu, ju, k) = em1(iu, ju, k) 
     .      + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdq00(k-1) * unom(k-1)
     .      - 0.5D0 * (dmdt(k  ) + dmdt(k+1)) * dqbdqm1(k  ) * unom(k-1)
      em2(iu, ju, k) = em2(iu, ju, k) 
     .      + 0.5D0 * (dmdt(k-1) + dmdt(k  )) * dqbdqm1(k-1) * unom(k-2)
c
c-----------------------------------------------------------------------
c     pressure gradient
c-----------------------------------------------------------------------
c
c     equation gm28 - gm32
c
      e00(iu, jr, k) = e00(iu, jr, k)
     .               + xmu * thet * rn(k) * rmum1(k) * (pg(k) - pg(k-1))
c
      e00(iu, jd, k) = e00(iu, jd, k)
     .                        + thet * rmu(k) * pgn(k  ) * dlpgdldn(k  )
      em1(iu, jd, k) = em1(iu, jd, k) 
     .                        - thet * rmu(k) * pgn(k-1) * dlpgdldn(k-1)
c
      e00(iu, jt, k) = e00(iu, jt, k) 
     .                        + thet * rmu(k) * pgn(k  ) * dlpgdltn(k  )
      em1(iu, jt, k) = em1(iu, jt, k) 
     .                        - thet * rmu(k) * pgn(k-1) * dlpgdltn(k-1)
c
c-----------------------------------------------------------------------
c     gravity
c-----------------------------------------------------------------------
c
c     equations gm33 - gm35
c
      ep1(iu, jr, k) = ep1(iu, jr, k)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu  (k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*  0.5D0 *  thet * (+ d(k  ) * rmu(k+1) ) * rn(k+1) 
c
      e00(iu, jr, k) = e00(iu, jr, k) 
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu  (k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*  0.5D0 *  thet * (- d(k  ) * rmu(k  ) 
     .                              + d(k-1) * rmu(k  ) ) * rn(k  )
     .      -  (        2.0D0 * xmu**2 * cpi * cgrav * xm (k)
     .                                        * thet * rn (k) / rmup1(k)
     .         )*  0.5D0 * (d (k-1) * dvol(k-1) + d (k) * dvol(k))
c
      em1(iu, jr, k) = em1(iu, jr, k)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu  (k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*  0.5D0 *  thet * (- d(k-1) * rmu(k-1) ) * rn(k-1) 
c
c     equation gm36
c
      e00(iu, jm, k) = e00(iu, jm, k) 
     .      -  (    thet * 2.0D0 * xmu * cpi * cgrav * xmen(k) / rmu(k)
     .         )*  0.5D0 * (d (k-1) * dvol(k-1) + d (k) * dvol(k))
c
c     equations gm37 & gm38
c
      e00(iu, jd, k) = e00(iu, jd, k)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu(k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*  0.5D0 *  thet * dn(k  ) * dvol(k  ) 
c
      em1(iu, jd, k) = em1(iu, jd, k)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu(k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*  0.5D0 *  thet * dn(k-1) * dvol(k-1) 
c
c-----------------------------------------------------------------------
c     (artificial) viscous force
c-----------------------------------------------------------------------
c
c     equations gm55 - gm64
c
      ep1(iu, jr, k) = ep1(iu, jr, k) 
     .                  +  r1     (k) * ( dr3dlrp1(k  ) *  qf     (k  ))
      e00(iu, jr, k) = e00(iu, jr, k) 
     .                  + dr1dlr00(k) * (  r3     (k  ) *  qf     (k  ) 
     .                                  -  r3     (k-1) *  qf     (k-1))
     .                  +  r1     (k) * ( dr3dlr00(k  ) *  qf     (k  ) 
     .                                  - dr3dlrp1(k-1) *  qf     (k-1))
      em1(iu, jr, k) = em1(iu, jr, k) 
     .                  -  r1     (k) * ( dr3dlr00(k-1) *  qf     (k-1))
c
      ep1(iu, jr, k) = ep1(iu, jr, k) 
     .                  +  r1     (k) * (  r3     (k  ) * dqfdlrp1(k  ))
      e00(iu, jr, k) = e00(iu, jr, k) 
     .                  +  r1     (k) * (  r3     (k  ) * dqfdlr00(k  ) 
     .                                  -  r3     (k-1) * dqfdlrp1(k-1))
      em1(iu, jr, k) = em1(iu, jr, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdlr00(k-1))
c
      ep1(iu, ju, k) = ep1(iu, ju, k) 
     .                  +  r1     (k) * (  r3     (k  ) * dqfdlup1(k  ))
      e00(iu, ju, k) = e00(iu, ju, k) 
     .                  +  r1     (k) * (  r3     (k  ) * dqfdlu00(k  ) 
     .                                  -  r3     (k-1) * dqfdlup1(k-1))
      em1(iu, ju, k) = em1(iu, ju, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdlu00(k-1))
c
      e00(iu, jt, k) = e00(iu, jt, k) 
     .                  +  r1     (k) * (  r3     (k  ) * dqfdlt00(k  ))
      em1(iu, jt, k) = em1(iu, jt, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdlt00(k-1))
c
      e00(iu, jd, k) = e00(iu, jd, k) 
     .                  +  r1     (k) * (  r3     (k  ) * dqfdld00(k  ))
      em1(iu, jd, k) = em1(iu, jd, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdld00(k-1))
c
c-----------------------------------------------------------------------
c     right hand side
c-----------------------------------------------------------------------
c
c     equation gm65
c
      rhs(iu, k) = 
     .             ( un(k) * (dn(k-1) * dvoln(k-1) + dn(k) * dvoln(k))
     .             - uo(k) * (do(k-1) * dvolo(k-1) + do(k) * dvolo(k)) )
     .                                                 / (2.0D0 * dtime)
     .
     .                 - 0.5D0 * ( (dmdt(k  ) + dmdt(k+1)) *   ub(k  )
     .                           - (dmdt(k-1) + dmdt(k  )) *   ub(k-1) )
     .
     .                                          + rmu(k) * (   pg(k  ) 
     .                                                     -   pg(k-1) )
     .
     .         +  (          2.0D0 * xmu * cpi * cgrav * xm(k) / rmu(k)
     .            + (1.0D0 - 0.5D0 * xmu)* g 
     .            )* 0.5D0 *   ( d (k) * dvol(k) + d (k-1) * dvol(k-1) )
     .
     .                                 + r1(k) * ( r3(k  ) *   qf(k  ) 
     .                                           - r3(k-1) *   qf(k-1) )
    6 continue
c
c=======================================================================
c     inner boundary condition
c=======================================================================
c
      k = ngrs
c
c-----------------------------------------------------------------------
c     eulerian or lagrangean boundary; equations bc16 or bc66
c-----------------------------------------------------------------------
c
      if (leibc .eq. 1 .or. leibc .eq. 2 .or. llibc .eq. 1) then
            rhs(iu,     k) = 0.0D0
            e00(iu, ju, k) = 1.0D0
      end if
c
c=======================================================================
c     outer boundary condition
c=======================================================================
c
      k = ngre + 1
c
c-----------------------------------------------------------------------
c     eulerian or lagrangean nontransmitting boundary;
c     equations bc32 or bc69
c-----------------------------------------------------------------------
c
      if (leobc .eq. 1 .or. leobc .eq. 2 .or. llobc .eq. 1) then
            rhs(iu,     k) = 0.0D0
            e00(iu, ju, k) = 1.0D0
      end if
c
c-----------------------------------------------------------------------
c     lagrangean boundary with specified external pressure;
c     equations bc73 and bc76
c-----------------------------------------------------------------------
c
c++++++++ llobc .eq. 2 +++++++++++++++++++++++++++++++++++++++++++ begin
      if (llobc .eq. 2 ) then
	    rx = (omega - 1.0D0) / (omega + 1.0D0)
            if   (omega .ge. 1.D+30)   rx = 1.0D0
c
c.......................................................................
c           time derivative; equations bc77 - bc80
c.......................................................................
c
            e00(iu, jr, k) = e00(iu, jr, k)
     .                           + un(k) * dn(k-1) * rmup1n(k  ) / dtime
            em1(iu, jr, k) = em1(iu, jr, k)
     .                           - un(k) * dn(k-1) * rmup1n(k-1) / dtime
c
            em1(iu, jd, k) = em1(iu, jd, k)
     .                           + un(k) * dn(k-1) * dvoln (k-1) / dtime
c
            e00(iu, ju, k) = e00(iu, ju, k)
     .                         + unom(k) * dn(k-1) * dvoln (k-1) / dtime
c
c.......................................................................
c           advection; equations bc81 - bc85
c.......................................................................
c
            em1(iu, jm, k) = em1(iu, jm, k) 
     .                        - rx * xmen (k-1) * ub(k-1) / dtime
c
            ep1(iu, ju, k) = ep1(iu, ju, k) 
     .                      + rx * dmdt(k-1) * unom(k+1) * dqbdqp2(k-1)
            e00(iu, ju, k) = e00(iu, ju, k) 
     .                      + rx * dmdt(k-1) * unom(k  ) * dqbdqp1(k-1)
            em1(iu, ju, k) = em1(iu, ju, k) 
     .                      + rx * dmdt(k-1) * unom(k-1) * dqbdq00(k-1)
            em2(iu, ju, k) = em2(iu, ju, k) 
     .                      + rx * dmdt(k-1) * unom(k-2) * dqbdqm1(k-1)
c
c.......................................................................
c           pressure gradient; equations bc86 - bc88
c.......................................................................
c
	    e00(iu, jr, k) = e00(iu, jr, k) 
     . + thet * xmu * rn(k) * 2.0D0 * rx * rmum1(k) * (pextr - pg(k-1) )
c
	    em1(iu, jd, k) = em1(iu, jd, k) 
     . - thet         * 2.0D0 * rx * rmu(k) * pgn(k-1) * dlpgdldn(k-1)
c
	    em1(iu, jt, k) = em1(iu, jt, k) 
     . - thet         * 2.0D0 * rx * rmu(k) * pgn(k-1) * dlpgdltn(k-1)
c
c.......................................................................
c           gravity; equations bc89 - bc92
c.......................................................................
c
            e00(iu, jr, k) = e00(iu, jr, k) 
     .      +  (         - 2.0D0 * xmu * cpi * cgrav * xm (k)
     .                                 * xmu * thet  * rn (k) / rmup1(k)
     .         )*               d (k-1) * dvol(k-1)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu  (k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*     thet * (+ d (k-1) *  rmu(k  )) * rn(k  )
c
            em1(iu, jr, k) = em1(iu, jr, k)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu  (k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*     thet * (- d (k-1) *  rmu(k-1)) * rn(k-1)
c
            e00(iu, jm, k) = e00(iu, jm, k) 
     .      -  (    thet * 2.0D0 * xmu * cpi * cgrav *xmen(k) / rmu  (k)
     .         )*               d (k-1) * dvol(k-1) 
c
            em1(iu, jd, k) = em1(iu, jd, k)
     .      +  (           2.0D0 * xmu * cpi * cgrav * xm (k) / rmu  (k)
     .          + (1.0D0 - 0.5D0 * xmu)* g 
     .         )*     thet *    dn(k-1) * dvol(k-1) 
c
c.......................................................................
c           viscous force; equations bc99 - bc103
c.......................................................................
c
            e00(iu, jr, k) = e00(iu, jr, k) 
     .                  - dr1dlr00(k) * (  r3     (k-1) *  qf     (k-1))
     .                  -  r1     (k) * ( dr3dlrp1(k-1) *  qf     (k-1) 
     .                                  +  r3     (k-1) * dqfdlrp1(k-1))
            em1(iu, jr, k) = em1(iu, jr, k) 
     .                  -  r1     (k) * ( dr3dlr00(k-1) *  qf     (k-1)
     .                                  +  r3     (k-1) * dqfdlr00(k-1))
            e00(iu, ju, k) = e00(iu, ju, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdlup1(k-1))
            em1(iu, ju, k) = em1(iu, ju, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdlu00(k-1))
c
            em1(iu, jt, k) = em1(iu, jt, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdlt00(k-1))
c
            em1(iu, jd, k) = em1(iu, jd, k) 
     .                  -  r1     (k) * (  r3     (k-1) * dqfdld00(k-1))
c
c.......................................................................
c           right-hand side; equation bc104
c.......................................................................
c
            rhs(iu, k) =         ( un(k) * dn(k-1) * dvoln(k-1) 
     .                           - uo(k) * do(k-1) * dvolo(k-1)) / dtime
     .
     .                                     + rx * dmdt(k-1) *   ub(k-1)
     .
     .                       + 2.0D0 * rx * rmu(k) * (pextr -   pg(k-1))
     .
     .         +  (          2.0D0 * xmu * cpi * cgrav * xm(k) / rmu(k)
     .            + (1.0D0 - 0.5D0 * xmu)* g 
     .            )                             *    d(k-1) * dvol(k-1)
     .
     .                                  - r1(k) *   r3(k-1) *   qf(k-1)
c
c++++++++ llobc .eq. 2 +++++++++++++++++++++++++++++++++++++++++++++ end
c
      end if
c
c-----------------------------------------------------------------------
c     eulerian or lagrangean transmitting boundary; equation bc41
c-----------------------------------------------------------------------
c
c++++++++ llobc .eq. 3 .or. llobc .eq. 3 +++++++++++++++++++++++++ begin
      if (leobc .eq. 3 .or. llobc .eq. 3 ) then
c
c.......................................................................
c           time derivative; equations bc42 - bc46
c.......................................................................
c
            e00(iu, jr, k) = e00(iu, jr, k)
     .     + (un(k-1) + un(k)) * dn(k-1) * rmup1n(k  ) / (2.0D0 * dtime)
            em1(iu, jr, k) = em1(iu, jr, k)
     .     - (un(k-1) + un(k)) * dn(k-1) * rmup1n(k-1) / (2.0D0 * dtime)
c
            em1(iu, jd, k) = em1(iu, jd, k)
     .     + (un(k-1) + un(k)) * dn(k-1) * dvoln (k-1) / (2.0D0 * dtime)
c
            e00(iu, ju, k) = e00(iu, ju, k)
     .     + (        unom(k)) * dn(k-1) * dvoln (k-1) / (2.0D0 * dtime)
            em1(iu, ju, k) = em1(iu, ju, k)
     .     + (unom(k-1)      ) * dn(k-1) * dvoln (k-1) / (2.0D0 * dtime)
c
c.......................................................................
c           advection; equations bc47 - bc50
c.......................................................................
c
            e00(iu, jm, k) = e00(iu, jm, k) + xmen(k  ) * u (k  )/ dtime
            em1(iu, jm, k) = em1(iu, jm, k) - xmen(k-1) * u (k-1)/ dtime
c
            e00(iu, ju, k) = e00(iu, ju, k) - dmdt(k  ) * unom(k  )*thet
            em1(iu, ju, k) = em1(iu, ju, k) + dmdt(k-1) * unom(k-1)*thet
c
c.......................................................................
c           sound speed * velocity gradient; equations  bc51 - bc56
c.......................................................................
c
            e00(iu, jr, k) = e00(iu, jr, k) 
     .             + 0.5D0 * thet * xmu * rn(k  ) * rmum1(k  ) * d (k-1)
     .                                  * (u  (k) -   u  (k-1))* as(k-1)
            em1(iu, jr, k) = em1(iu, jr, k) 
     .             + 0.5D0 * thet * xmu * rn(k-1) * rmum1(k-1) * d (k-1)
     .                                  * (u  (k) -   u  (k-1))* as(k-1)
c
            e00(iu, ju, k) = e00(iu, ju, k) 
     .                           + 0.5D0 * (rmu(k) + rmu(k-1)) * d (k-1)
     .                           *  thet * (unom(k)          ) * as(k-1)
            em1(iu, ju, k) = em1(iu, ju, k) 
     .                           + 0.5D0 * (rmu(k) + rmu(k-1)) * d (k-1)
     .                           *  thet * (      - unom(k-1)) * as(k-1)
c
            em1(iu, jt, k) = em1(iu, jt, k) 
     .                           + 0.5D0 * (rmu(k) + rmu(k-1)) * d (k-1)
     .                           * thet  * (u  (k) - u  (k-1)) * as(k-1)
     .           * 0.5D0 * thet *   (pgn(k-1) / pg(k-1)) * dlpgdltn(k-1)
c
            em1(iu, jd, k) = em1(iu, jd, k)
     .                           + 0.5D0 * (rmu(k) + rmu(k-1)) * d (k-1)
     .                           * thet  * (u  (k) - u  (k-1)) * as(k-1)
     .           * 0.5D0 * thet * ( (pgn(k-1) / pg(k-1)) * dlpgdldn(k-1)
     .                            +  dn (k-1) / d (k-1) )
c
c.......................................................................
c           right-hand side; equation bc41
c.......................................................................
c
	    rhs(iu, k) = 
     .                      ( (un(k-1) + un(k)) * dn(k-1) * dvoln(k-1)
     .                      - (uo(k-1) + uo(k)) * do(k-1) * dvolo(k-1) )
     .                                                 / (2.0D0 * dtime)
     .
     .                                           - dmdt(k  )  *  u(k  ) 
     .                                           + dmdt(k-1)  *  u(k-1)
     .
     .                         + 0.5D0 * (rmu(k) +  rmu(k-1)) * d (k-1)
     .                                 * (u  (k) -  u  (k-1)) * as(k-1)
c
c++++++++ llobc .eq. 3 .or. llobc .eq. 3 +++++++++++++++++++++++++++ end
c
      end if
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 7 k = 1, ngrs - 2
      rhs(iu,     k) = 0.0D0
      e00(iu, ju, k) = 1.0D0
    7 continue
c
      do 8 k = ngre + 3, mgr
      rhs(iu,     k) = 0.0D0
      e00(iu, ju, k) = 1.0D0
    8 continue
c
c-----------------------------------------------------------------------
c     eulerian inner boundary
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
c     zero flux; equations pz7 & pz8
c
      if (leibc .eq. 1) then
            rhs(iu,     k) = un  (k+2) + un  (k)
            ep2(iu, ju, k) = unom(k+2)
            e00(iu, ju, k) =             unom(k)
      end if
c
c     nonzero flux; equations pz17 & pz18
c
      if (leibc .eq. 2) then
            rhs(iu,     k) = 0.0D0
            e00(iu, ju, k) = 1.0D0
      end if
c
c-----------------------------------------------------------------------
c     lagrangean inner boundary
c-----------------------------------------------------------------------
c
c     equations pz27 & pz28
c
      if (llibc .eq. 1) then
            rhs(iu,     k) = un  (k) - un  (k+1)
            ep1(iu, ju, k) =         - unom(k+1)
            e00(iu, ju, k) = unom(k)
      end if
c
c-----------------------------------------------------------------------
c     eulerian outer boundary 
c-----------------------------------------------------------------------
c
      k = ngre + 2
c
c     zero flux; equations pz45 & pz46
c
      if (leobc .eq. 1) then
            rhs(iu,     k) = un  (k-2) + un  (k)
            e00(iu, ju, k) =           + unom(k)
            em2(iu, ju, k) = unom(k-2)
      end if
c
c     nonzero flux; equations pz55 & pz56
c
      if (leobc .eq. 2) then
            rhs(iu,     k) = 0.0D0
            e00(iu, ju, k) = 1.0D0
      end if
c
c     transmitting boundary; equations pz65 & pz66
c
      if (leobc .eq. 3) then
            rhs(iu,     k) = - un  (k-1) + un  (k)
            e00(iu, ju, k) =             + unom(k)
            em1(iu, ju, k) = - unom(k-1)
      end if
c
c-----------------------------------------------------------------------
c     lagrangean outer boundary
c-----------------------------------------------------------------------
c
c     equations pz75 & pz76
c
      if (llobc .gt. 0) then
c           rhs(iu,     k) = u   (k) - ( rn(k) - ro(k) ) / dtime
c           e00(iu, jr, k) =         -   rn(k)           / dtime
c           e00(iu, ju, k) = unom(k) * thet
            rhs(iu,     k) = - un  (k-1) + un  (k)
            e00(iu, ju, k) =             + unom(k)
            em1(iu, ju, k) = - unom(k-1)
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
