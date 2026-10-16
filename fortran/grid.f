      subroutine grid
c
c***********************************************************************
c     adaptive grid definition
c.......................................................................
c     called by: matgen, gridinit
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c-----------------------------------------------------------------------
c
      alpha = alph * ( 1.0D0 + alph )
      gamma = ( tau / dtime )**ibet
c
c***********************************************************************
c     lagrangean grid; equations lg1 & ag15 - ag18
c***********************************************************************
c
      if ( lgrid .eq. 3 ) then
c
             do 301 k = ngrs, ngre + 1
             rhs( ir,     k ) = (rn(k) - ro(k)) / dtime - u   (k)
             e00( ir, jr, k ) =  rn(k)          / dtime
             e00( ir, ju, k ) =                         - unom(k) * thet
  301        continue
             end if
c
c***********************************************************************
c     eulerian grid; equations eg1 & ag12 - ag14
c***********************************************************************
c
      if ( lgrid .eq. 2 ) then
c
             do 201 k = ngrs, ngre + 1
             rhs( ir,     k ) = rn(k) - ro(k)
             e00( ir, jr, k ) = rn(k)
  201        continue
             end if
c
c***********************************************************************
c     adaptive grid; equations ag1 - ag6                          
c***********************************************************************
c
c+++++++++ lgrid .eq. 1 +++++++++++++++++++++++++++++++++++++++++++begin
c
      if ( lgrid .eq. 1 ) then
c
c=======================================================================
c     abscissa
c=======================================================================
c-----------------------------------------------------------------------
c     linear grid concentration	
c-----------------------------------------------------------------------
c
c     equations ag7 & ag19 - ag21
c
      if (ladx .eq. 1) then
           do 1 k = ngrs, ngre
           xnu     (k) =   xscale            / ( rn(k+1) - rn(k) )
           dnudlr00(k) = + xscale * rn(k  )  / ( rn(k+1) - rn(k) )**2
           dnudlrp1(k) = - xscale * rn(k+1)  / ( rn(k+1) - rn(k) )**2
           xnuo    (k) =   xscale            / ( ro(k+1) - ro(k) )
    1      continue
      end if
c
c-----------------------------------------------------------------------
c     logarithmic grid concentration	
c-----------------------------------------------------------------------
c
c     equations ag8 & ag22 - ag24
c
      if (ladx .eq. 2) then
           do 2 k = ngrs, ngre
           xnu     (k) = ( rn(k+1) + rn(k) ) / ( rn(k+1) - rn(k) ) /2.D0
           dnudlr00(k) = + rn(k+1) * rn(k)   / ( rn(k+1) - rn(k) )**2 
           dnudlrp1(k) = - rn(k+1) * rn(k)   / ( rn(k+1) - rn(k) )**2 
           xnuo    (k) = ( ro(k+1) + ro(k) ) / ( ro(k+1) - ro(k) ) /2.D0
    2      continue
      end if
c
c-----------------------------------------------------------------------
c     spatial operator; equations ag5 & ag32, ag34 - ag37
c-----------------------------------------------------------------------
c
      do 3 k = ngrs + 1, ngre - 1
c
      xnt     (k) =                        xnu     (k) - alpha * 
     .            ( xnu     (k-1) - 2.D0 * xnu     (k) + xnu     (k+1) )
c
      dntdlrm1(k) = dnudlr00(k-1)
      dntdlr00(k) = dnudlrp1(k-1) - 2.D0 * dnudlr00(k) 
      dntdlrp1(k) =               - 2.D0 * dnudlrp1(k) + dnudlr00(k+1)
      dntdlrp2(k) =                                    + dnudlrp1(k+1)
c
      dntdlrm1(k) = - alpha * dntdlrm1(k)
      dntdlr00(k) = - alpha * dntdlr00(k) + dnudlr00(k) 
      dntdlrp1(k) = - alpha * dntdlrp1(k) + dnudlrp1(k)
      dntdlrp2(k) = - alpha * dntdlrp2(k)
c
      xnto    (k) =                        xnuo    (k) - alpha * 
     .            ( xnuo    (k-1) - 2.D0 * xnuo    (k) + xnuo    (k+1) )
c
c-----------------------------------------------------------------------
c     temporal operator; equations ag4 & ag33, ag38 - ag41
c-----------------------------------------------------------------------
c
      xnc     (k) = ( 1.0D0 + gamma ) * xnt     (k)
     .                      - gamma   * xnto    (k)
c
      dncdlrm1(k) = ( 1.0D0 + gamma ) * dntdlrm1(k)
      dncdlr00(k) = ( 1.0D0 + gamma ) * dntdlr00(k)
      dncdlrp1(k) = ( 1.0D0 + gamma ) * dntdlrp1(k)
      dncdlrp2(k) = ( 1.0D0 + gamma ) * dntdlrp2(k)
    3 continue
c
c=======================================================================
c     ordinates:
c     l = 1 => xm     l = 2 => d      l = 3 => t          l = 4 => er
c     l = 5 => pg     l = 6 => eg     l = 7 => chi        l = 8 => qx
c     l = 9 => fr     l =10 => bri
c=======================================================================
c-----------------------------------------------------------------------
c     initialization
c-----------------------------------------------------------------------
c
      do 4 k = 1, mgr
      ss(k) = 0.0D0
    4 continue
c
      do 8 m = 1, meqn
      do 5 k = 1, mgr
      dssdlx00(k, m) = 0.0D0
      dssdlxp1(k, m) = 0.0D0
      dssdlxp2(k, m) = 0.0D0
    5 continue
c
      do 7 l = 1, mad
      do 6 k = 1, mgr
       cs     (k, l   ) = 0.0D0
      dcsdlx00(k, l, m) = 0.0D0
      dcsdlxp1(k, l, m) = 0.0D0
      dcsdlxp2(k, l, m) = 0.0D0
    6 continue
    7 continue
    8 continue
c
c-----------------------------------------------------------------------
c     linear resolution
c-----------------------------------------------------------------------
c.......................................................................
c     mass; equations ag59 - ag61
c.......................................................................
c
      if (wt(1) .gt. 0.0D0 .and. lady(1) .eq. 1 .or.
     .    wt(11).gt. 0.0D0 .and. lady(11).eq. 1 ) then
            l = 1
            do 11 k = ngrs, ngre - 1
             cs     (k,l   ) = ( xmen(k+1) - xmen(k) ) / yscl(l)
            dcsdlx00(k,l,im) = (           - xmen(k) ) / yscl(l)
            dcsdlxp1(k,l,im) = ( xmen(k+1)           ) / yscl(l)
   11       continue
      end if
c
c.......................................................................
c     density; equations ag62 - ag64
c.......................................................................
c
      if (wt(2) .gt. 0.0D0 .and. lady(2) .eq. 1) then
            l = 2
            do 12 k = ngrs, ngre - 1
             cs     (k,l   ) = ( dn(k+1) - dn(k) ) / yscl(l)
            dcsdlx00(k,l,id) = (         - dn(k) ) / yscl(l)
            dcsdlxp1(k,l,id) = ( dn(k+1)         ) / yscl(l)
   12       continue
      end if
c
c.......................................................................
c     gas temperature; equations ag65 - ag67
c.......................................................................
c
      if (wt(3) .gt. 0.0D0 .and. lady(3) .eq. 1) then
            l = 3
            do 13 k = ngrs, ngre - 1
             cs     (k,l   ) = ( tn(k+1) - tn(k) ) / yscl(l)
            dcsdlx00(k,l,it) = (         - tn(k) ) / yscl(l)
            dcsdlxp1(k,l,it) = ( tn(k+1)         ) / yscl(l)
   13       continue
      end if
c
c.......................................................................
c     radiation energy; equations ag68 - ag70
c.......................................................................
c
      if (wt(4) .gt. 0.0D0 .and. lady(4) .eq. 1) then
            l = 4
            do 14 k = ngrs, ngre - 1
             cs     (k,l   ) = ( ern(k+1) - ern(k) ) / yscl(l)
            dcsdlx00(k,l,ie) = (          - ern(k) ) / yscl(l)
            dcsdlxp1(k,l,ie) = ( ern(k+1)          ) / yscl(l)
   14       continue
      end if
c
c.......................................................................
c     gas pressure; equations ag71 - ag75 
c.......................................................................
c
      if (wt(5) .gt. 0.0D0 .and. lady(5) .eq. 1) then
            l = 5
            do 15 k = ngrs, ngre - 1
             cs     (k,l   ) = ( pgn(k+1) -   pgn   (k  )) / yscl(l)
            dcsdlx00(k,l,id) = - pgn(k  ) * dlpgdldn(k  )  / yscl(l)
            dcsdlxp1(k,l,id) = + pgn(k+1) * dlpgdldn(k+1)  / yscl(l)
            dcsdlx00(k,l,it) = - pgn(k  ) * dlpgdltn(k  )  / yscl(l)
            dcsdlxp1(k,l,it) = + pgn(k+1) * dlpgdltn(k+1)  / yscl(l)
   15       continue
      end if
c
c.......................................................................
c     gas energy density; equations ag76 - ag80
c.......................................................................
c
      if (wt(6) .gt. 0.0D0 .and. lady(6) .eq. 1) then
            l = 6
            do 16 k = ngrs, ngre - 1
             cs     (k,l   ) = ( egn(k+1) -   egn   (k  )) / yscl(l)
            dcsdlx00(k,l,id) = - egn(k  ) * dlegdldn(k  )  / yscl(l)
            dcsdlxp1(k,l,id) = + egn(k+1) * dlegdldn(k+1)  / yscl(l)
            dcsdlx00(k,l,it) = - egn(k  ) * dlegdltn(k  )  / yscl(l)
            dcsdlxp1(k,l,it) = + egn(k+1) * dlegdltn(k+1)  / yscl(l)
   16       continue
      end if
c
c.......................................................................
c     opacity; equations ag81 - ag85
c.......................................................................
c
      if (wt(7) .gt. 0.0D0 .and. lady(7) .eq. 1) then
            l = 7
            do 17 k = ngrs, ngre - 1
             cs     (k,l   ) = ( chifn(k+1) -   chifn (k  )) / yscl(l)
            dcsdlx00(k,l,id) = - chifn(k  ) * dlcfdldn(k  )  / yscl(l)
            dcsdlxp1(k,l,id) = + chifn(k+1) * dlcfdldn(k+1)  / yscl(l)
            dcsdlx00(k,l,it) = - chifn(k  ) * dlcfdltn(k  )  / yscl(l)
            dcsdlxp1(k,l,it) = + chifn(k+1) * dlcfdltn(k+1)  / yscl(l)
   17       continue
      end if
c
c.......................................................................
c     viscosity; equations ag99 - ag109
c.......................................................................
c
      if (wt(8) .gt. 0.0D0 .and. lady(8) .eq. 1) then
            l = 8
            do 18 k = ngrs, ngre - 1
             cs     (k,l   ) = ( qx     (k+1) -  qx     (k  )) / yscl(l)
            dcsdlx00(k,l,ir) =                - dqxdlr00(k  )  / yscl(l)
            dcsdlxp1(k,l,ir) = (dqxdlr00(k+1) - dqxdlrp1(k  )) / yscl(l)
            dcsdlxp2(k,l,ir) =  dqxdlrp1(k+1)                  / yscl(l)
            dcsdlx00(k,l,iu) =                - dqxdlu00(k  )  / yscl(l)
            dcsdlxp1(k,l,iu) = (dqxdlu00(k+1) - dqxdlup1(k  )) / yscl(l)
            dcsdlxp2(k,l,iu) =  dqxdlup1(k+1)                  / yscl(l)
            dcsdlx00(k,l,id) =                - dqxdld00(k  )  / yscl(l)
            dcsdlxp1(k,l,id) =  dqxdld00(k+1)                  / yscl(l)
            dcsdlx00(k,l,it) =                - dqxdlt00(k  )  / yscl(l)
            dcsdlxp1(k,l,it) =  dqxdlt00(k+1)                  / yscl(l)
   18       continue
      end if
c
c.......................................................................
c     radiation flux
c.......................................................................
c
      if (wt(9) .gt. 0.0D0 .and. lady(9) .eq. 1) then
            l = 9
            do 19 k = ngrs, ngre - 1
             cs     (k,l   ) = ( frn(k+1) - frn(k) ) / yscl(l)
            dcsdlx00(k,l,if) = (          - frn(k) ) / yscl(l)
            dcsdlxp1(k,l,if) = ( frn(k+1)          ) / yscl(l)
   19       continue
      end if
c
c.......................................................................
c     luminosity or velocity
c.......................................................................
c
      if (wt(10) .gt. 0.0D0 .and. lady(10) .eq. 1) then
            l = 10
            do 20 k = ngrs, ngre - 1
             cs     (k,l   ) = ( bri(k+1) - bri(k) ) / yscl(l)
            dcsdlx00(k,l,jb) = (          - bri(k) ) / yscl(l)
            dcsdlxp1(k,l,jb) = ( bri(k+1)          ) / yscl(l)
   20       continue
      end if
c 
c
c-----------------------------------------------------------------------
c     logarithmic resolution
c-----------------------------------------------------------------------
c.......................................................................
c     interior mass; equations ag110 - ag112
c.......................................................................
c
      if( wt(1) .gt. 0.0D0 .and. lady(1) .eq. 2 ) then
            l = 1
                                     xmscl = yscl(1)
            do 21 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( xmen(k+1) - xmen(k) )
     .                       /( 2.D0*xmscl - ( xmen(k+1) + xmen(k) ))
            dcsdlx00(k,l,im) =  2.D0*xmen(k)*( xmen(k+1) - xmscl)
     .                       /( 2.D0*xmscl - ( xmen(k+1) + xmen(k) ))**2
            dcsdlxp1(k,l,im) =  2.D0*xmen(k+1) * ( xmscl - xmen(k) )
     .                       /( 2.D0*xmscl - ( xmen(k+1) + xmen(k) ))**2
   21       continue
      end if
c
c.......................................................................
c     exterior mass; equations ag110 - ag112
c.......................................................................
c
      if( wt(11) .gt. 0.0D0 .and. lady(11) .eq. 2 ) then
            l = 11
            do 210 k = ngrs, ngre - 1
             cs     (k,l   ) =             ( xmen(k+1) - xmen(k) )
     .                                   / ( xmen(k+1) + xmen(k) )
            dcsdlx00(k,l,im) =             - xmen(k+1) * xmen(k)  * 2.D0
     .                                   / ( xmen(k+1) + xmen(k) )**2
            dcsdlxp1(k,l,im) =             + xmen(k+1) * xmen(k)  * 2.D0
     .                                   / ( xmen(k+1) + xmen(k) )**2
  210       continue
      end if
c
c.......................................................................
c     density; equations ag113 - ag115
c.......................................................................
c
      if (wt(2) .gt. 0.0D0 .and. lady(2) .eq. 2) then
            l = 2
            do 22 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( dn(k+1) - dn(k) )
     .                                     / ( dn(k+1) + dn(k) )
            dcsdlx00(k,l,id) =               - dn(k+1) * dn(k)  * 2.D0
     .                                     / ( dn(k+1) + dn(k) )**2
            dcsdlxp1(k,l,id) =               + dn(k+1) * dn(k)  * 2.D0
     .                                     / ( dn(k+1) + dn(k) )**2
   22       continue
      end if
c
c.......................................................................
c     gas temperature; equations ag116 - ag118
c.......................................................................
c
      if (wt(3) .gt. 0.0D0 .and. lady(3) .eq. 2) then
            l = 3
            do 23 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( tn(k+1) - tn(k) )
     .                                     / ( tn(k+1) + tn(k) )
            dcsdlx00(k,l,it) =               - tn(k+1) * tn(k)  * 2.D0
     .                                     / ( tn(k+1) + tn(k) )**2
            dcsdlxp1(k,l,it) =               + tn(k+1) * tn(k)  * 2.D0
     .                                     / ( tn(k+1) + tn(k) )**2
   23       continue
      end if
c
c.......................................................................
c     radiation energy; equations ag119 - ag121
c.......................................................................
c
      if (wt(4) .gt. 0.0D0 .and. lady(4) .eq. 2) then
            l = 4
            do 24 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( ern(k+1) - ern(k) )
     .                                     / ( ern(k+1) + ern(k) )
            dcsdlx00(k,l,ie) =               - ern(k+1) * ern(k)  * 2.D0
     .                                     / ( ern(k+1) + ern(k) )**2
            dcsdlxp1(k,l,ie) =               + ern(k+1) * ern(k)  * 2.D0
     .                                     / ( ern(k+1) + ern(k) )**2
   24       continue
      end if
c
c.......................................................................
c     gas pressure; equations ag122 - ag126
c.......................................................................
c
      if (wt(5) .gt. 0.0D0 .and. lady(5) .eq. 2) then
            l = 5
            do 25 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( pgn(k+1) - pgn(k) )
     .                                     / ( pgn(k+1) + pgn(k) ) 
            dcsdlx00(k,l,it) =               - pgn(k+1) * pgn(k)  * 2.D0
     .                     * dlpgdltn(k  ) / ( pgn(k+1) + pgn(k) )**2
            dcsdlxp1(k,l,it) =               + pgn(k+1) * pgn(k)  * 2.D0
     .                     * dlpgdltn(k+1) / ( pgn(k+1) + pgn(k) )**2
            dcsdlx00(k,l,id) =               - pgn(k+1) * pgn(k)  * 2.D0
     .                     * dlpgdldn(k  ) / ( pgn(k+1) + pgn(k) )**2
            dcsdlxp1(k,l,id) =               + pgn(k+1) * pgn(k)  * 2.D0
     .                     * dlpgdldn(k+1) / ( pgn(k+1) + pgn(k) )**2
   25       continue
      end if
c
c.......................................................................
c     gas energy density; equations ag127 - ag131
c.......................................................................
c
      if (wt(6) .gt. 0.0D0 .and. lady(6) .eq. 2) then
            l = 6
            do 26 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( egn(k+1) - egn(k) )
     .                                     / ( egn(k+1) + egn(k) ) 
            dcsdlx00(k,l,it) =               - egn(k+1) * egn(k)  * 2.D0
     .                     * dlegdltn(k  ) / ( egn(k+1) + egn(k) )**2
            dcsdlxp1(k,l,it) =               + egn(k+1) * egn(k)  * 2.D0
     .                     * dlegdltn(k+1) / ( egn(k+1) + egn(k) )**2
            dcsdlx00(k,l,id) =               - egn(k+1) * egn(k)  * 2.D0
     .                     * dlegdldn(k  ) / ( egn(k+1) + egn(k) )**2
            dcsdlxp1(k,l,id) =               + egn(k+1) * egn(k)  * 2.D0
     .                     * dlegdldn(k+1) / ( egn(k+1) + egn(k) )**2
   26       continue
      end if
c
c.......................................................................
c     opacity; equations ag132 - ag136
c.......................................................................
c
      if (wt(7) .gt. 0.0D0 .and. lady(7) .eq. 2) then
            l = 7
            do 27 k = ngrs, ngre - 1
             cs     (k,l   ) =           ( chifn(k+1) - chifn(k) )
     .                                 / ( chifn(k+1) + chifn(k) ) 
            dcsdlx00(k,l,it) =           - chifn(k+1) * chifn(k)  * 2.D0
     .                 * dlcfdltn(k  ) / ( chifn(k+1) + chifn(k) )**2
            dcsdlxp1(k,l,it) =           + chifn(k+1) * chifn(k)  * 2.D0
     .                 * dlcfdltn(k+1) / ( chifn(k+1) + chifn(k) )**2
            dcsdlx00(k,l,id) =           - chifn(k+1) * chifn(k)  * 2.D0
     .                 * dlcfdldn(k  ) / ( chifn(k+1) + chifn(k) )**2
            dcsdlxp1(k,l,id) =           + chifn(k+1) * chifn(k)  * 2.D0
     .                 * dlcfdldn(k+1) / ( chifn(k+1) + chifn(k) )**2
   27       continue
      end if
c
c.......................................................................
c     viscosity; equations ag137 - ag147
c.......................................................................
c
      if (wt(8) .gt. 0.0D0 .and. lady(8) .eq. 2) then
            l = 8
            do 28 k = ngrs, ngre - 1
             cs     (k,l   ) =                      (qx(k+1) - qx(k))
     .                                            / (qx(k+1) + qx(k)) 
            dcsdlx00(k,l,ir) = -2.D0*dqxdlr00(k  ) * qx(k+1)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlxp1(k,l,ir) =  2.D0*dqxdlr00(k+1) *           qx(k)
     .                                            / (qx(k+1) + qx(k))**2
     .                         -2.D0*dqxdlrp1(k  ) * qx(k+1)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlxp2(k,l,ir) =  2.D0*dqxdlrp1(k+1) *           qx(k)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlx00(k,l,iu) =  2.D0*dqxdlu00(k  ) * qx(k+1)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlxp1(k,l,iu) =  2.D0*dqxdlu00(k+1) *           qx(k) 
     .                                            / (qx(k+1) + qx(k))**2
     .                         -2.D0*dqxdlup1(k  ) * qx(k+1)  
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlxp2(k,l,iu) =  2.D0*dqxdlup1(k+1) *           qx(k)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlx00(k,l,id) = -2.D0*dqxdld00(k  ) * qx(k+1)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlxp1(k,l,id) =  2.D0*dqxdld00(k+1) *           qx(k)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlx00(k,l,it) = -2.D0*dqxdlt00(k  ) * qx(k+1)
     .                                            / (qx(k+1) + qx(k))**2
            dcsdlxp1(k,l,it) =  2.D0*dqxdlt00(k+1) *           qx(k) 
     .                                            / (qx(k+1) + qx(k))**2
   28       continue
      end if
c
c.......................................................................
c     radiative flux
c.......................................................................
c
      if (wt(9) .gt. 0.0D0 .and. lady(9) .eq. 2) then
            l = 9
            do 29 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( frn(k+1) - frn(k) )
     .                                     / ( frn(k+1) + frn(k) )
            dcsdlx00(k,l,if) =               - frn(k+1) * frn(k)  * 2.D0
     .                                     / ( frn(k+1) + frn(k) )**2
            dcsdlxp1(k,l,if) =               + frn(k+1) * frn(k)  * 2.D0
     .                                     / ( frn(k+1) + frn(k) )**2
   29       continue
      end if
c
c.......................................................................
c     luminosity or velocity
c.......................................................................
c
      if (wt(10) .gt. 0.0D0 .and. lady(10) .eq. 2) then
            l = 10
            do 30 k = ngrs, ngre - 1
             cs     (k,l   ) =               ( bri(k+1) - bri(k) )
     .                                     / ( bri(k+1) + bri(k) )
            dcsdlx00(k,l,jb) =               - bri(k+1) * bri(k)  * 2.D0
     .                                     / ( bri(k+1) + bri(k) )**2
            dcsdlxp1(k,l,jb) =               + bri(k+1) * bri(k)  * 2.D0
     .                                     / ( bri(k+1) + bri(k) )**2
   30       continue
      end if
c
c-----------------------------------------------------------------------
c     harmonic resolution
c-----------------------------------------------------------------------
c.......................................................................
c     interior mass; equations ag148 - ag150
c.......................................................................
c
      if (wt(1) .gt. 0.0D0 .and. lady(1) .eq. 3) then
            l = 1
                              xmscl = yscl(1)
            do 31 k = ngrs, ngre - 1
            cs     (k,l   ) =  
     .                       (xmscl - xmen(k+1)) / (xmscl - xmen(k  ))
     .                     - (xmscl - xmen(k  )) / (xmscl - xmen(k+1))
           dcsdlx00(k,l,im) =       ( xmen(k  )  / (xmscl - xmen(k  )) )
     .                   * ( (xmscl - xmen(k+1)) / (xmscl - xmen(k  ))
     .                     + (xmscl - xmen(k  )) / (xmscl - xmen(k+1)) )
           dcsdlxp1(k,l,im) =     - ( xmen(k+1)  / (xmscl - xmen(k+1)) )
     .                   * ( (xmscl - xmen(k+1)) / (xmscl - xmen(k  ))
     .                     + (xmscl - xmen(k  )) / (xmscl - xmen(k+1)) )

   31       continue
      end if
c.......................................................................
c     exterior mass; equations ag148 - ag150
c.......................................................................
c
      if (wt(11) .gt. 0.0D0 .and. lady(11) .eq. 3) then
            l = 11
            do 310 k = ngrs, ngre - 1
             cs     (k,l   ) =   xmen(k+1)/xmen(k) - xmen(k)/xmen(k+1)
            dcsdlx00(k,l,im) = - xmen(k+1)/xmen(k) - xmen(k)/xmen(k+1)
            dcsdlxp1(k,l,im) =   xmen(k+1)/xmen(k) + xmen(k)/xmen(k+1)
  310       continue
      end if
c
c
c.......................................................................
c     density; equations ag151 - ag153
c.......................................................................
c
      if (wt(2) .gt. 0.0D0 .and. lady(2) .eq. 3) then
            l = 2
            do 32 k = ngrs, ngre - 1
             cs     (k,l   ) =   dn(k+1) / dn(k) - dn(k) / dn(k+1)
            dcsdlx00(k,l,id) = - dn(k+1) / dn(k) - dn(k) / dn(k+1)
            dcsdlxp1(k,l,id) =   dn(k+1) / dn(k) + dn(k) / dn(k+1)
   32       continue
      end if
c
c.......................................................................
c     temperature; equations ag154 - ag156
c.......................................................................
c
      if (wt(3) .gt. 0.0D0 .and. lady(3) .eq. 3) then
            l = 3
            do 33 k = ngrs, ngre - 1
             cs     (k,l   ) =   tn(k+1) / tn(k) - tn(k) / tn(k+1)
            dcsdlx00(k,l,it) = - tn(k+1) / tn(k) - tn(k) / tn(k+1)
            dcsdlxp1(k,l,it) =   tn(k+1) / tn(k) + tn(k) / tn(k+1)
   33       continue
      end if
c
c.......................................................................
c     radiation energy density; equations ag157 - ag159
c.......................................................................
c
      if (wt(4) .gt. 0.0D0 .and. lady(4) .eq. 3) then
            l = 4
            do 34 k = ngrs, ngre - 1
             cs     (k,l   ) =   ern(k+1) / ern(k) - ern(k) / ern(k+1)
            dcsdlx00(k,l,ie) = - ern(k+1) / ern(k) - ern(k) / ern(k+1)
            dcsdlxp1(k,l,ie) =   ern(k+1) / ern(k) + ern(k) / ern(k+1)
   34       continue
      end if
c
c.......................................................................
c     gas pressure; equations ag160 - ag164
c.......................................................................
c
      if (wt(5) .gt. 0.0D0 .and. lady(5) .eq. 3) then
            l = 5
            do 35 k = ngrs, ngre - 1
             cs     (k,l   ) =     pgn(k+1)/pgn(k) - pgn(k)/pgn(k+1)
c
            dcsdlx00(k,l,id) = 
     .         - dlpgdldn(k  ) * ( pgn(k+1)/pgn(k) + pgn(k)/pgn(k+1) )
c
            dcsdlxp1(k,l,id) =  
     .           dlpgdldn(k+1) * ( pgn(k+1)/pgn(k) + pgn(k)/pgn(k+1) )
c
            dcsdlx00(k,l,it) = 
     .         - dlpgdltn(k  ) * ( pgn(k+1)/pgn(k) + pgn(k)/pgn(k+1) )
c
            dcsdlxp1(k,l,it) =
     .           dlpgdltn(k+1) * ( pgn(k+1)/pgn(k) + pgn(k)/pgn(k+1) )
   35       continue
      end if
c
c.......................................................................
c     gas energy density; equations ag165 - ag169
c.......................................................................
c
      if (wt(6) .gt. 0.0D0 .and. lady(6) .eq. 3) then
            l = 6
            do 36 k = ngrs, ngre - 1
             cs     (k,l   ) =     egn(k+1)/egn(k) - egn(k)/egn(k+1)
c
            dcsdlx00(k,l,id) = 
     .         - dlegdldn(k  ) * ( egn(k+1)/egn(k) + egn(k)/egn(k+1) )
c
            dcsdlxp1(k,l,id) =  
     .           dlegdldn(k+1) * ( egn(k+1)/egn(k) + egn(k)/egn(k+1) )
c
            dcsdlx00(k,l,it) = 
     .         - dlegdltn(k  ) * ( egn(k+1)/egn(k) + egn(k)/egn(k+1) )
c
            dcsdlxp1(k,l,it) =
     .           dlegdltn(k+1) * ( egn(k+1)/egn(k) + egn(k)/egn(k+1) )
   36       continue
      end if
c
c.......................................................................
c     opacity; equations ag170 - ag174
c.......................................................................
c
      if (wt(7) .gt. 0.0D0 .and. lady(7) .eq. 3) then
        l = 7
        do 37 k = ngrs, ngre - 1
        cs     (k,l   ) =   chifn(k+1)/chifn(k) - chifn(k)/chifn(k+1)
c
       dcsdlx00(k,l,id) = 
     .   - dlcfdldn(k  ) * (chifn(k+1)/chifn(k) + chifn(k)/chifn(k+1))
c
       dcsdlxp1(k,l,id) =  
     .     dlcfdldn(k+1) * (chifn(k+1)/chifn(k) + chifn(k)/chifn(k+1))
c
       dcsdlx00(k,l,it) = 
     .   - dlcfdltn(k  ) * (chifn(k+1)/chifn(k) + chifn(k)/chifn(k+1))
c
       dcsdlxp1(k,l,it) =
     .     dlcfdltn(k+1) * (chifn(k+1)/chifn(k) + chifn(k)/chifn(k+1))
   37       continue
      end if
c
c.......................................................................
c     viscosity; equations ag175 - ag185
c.......................................................................
c
      if (wt(8) .gt. 0.0D0 .and. lady(8) .eq. 3) then
           l = 8
          do 38 k = ngrs, ngre - 1
           cs     (k,l   ) =           (qx(k+1)/qx(k) - qx(k)/qx(k+1))
c
          dcsdlx00(k,l,ir) = 
     .     - (dqxdlr00(k  )/qx(k  )) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlxp1(k,l,ir) =  
     .       (dqxdlr00(k+1)/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
     .     - (dqxdlrp1(k  )/qx(k  )) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlxp2(k,l,ir) = 
     .       (dqxdlrp1(k+1)/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlx00(k,l,id) = 
     .    - (dqxdld00(k  )/qx(k  )) *  (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlxp1(k,l,id) =
     .       (dqxdld00(k+1)/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlx00(k,l,iu) = 
     .     - (dqxdlu00(k  )/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlxp1(k,l,iu) = 
     .       (dqxdlu00(k+1)/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
     .     - (dqxdlup1(k  )/qx(k  )) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlxp2(k,l,iu) =  
     .       (dqxdlup1(k+1)/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlx00(k,l,it) = 
     .     - (dqxdlt00(k  )/qx(k  )) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
c
          dcsdlxp1(k,l,it) = 
     .       (dqxdlt00(k+1)/qx(k+1)) * (qx(k+1)/qx(k) + qx(k)/qx(k+1))
   38     continue
      end if
c
c.......................................................................
c     radiative flux
c.......................................................................
c
      if (wt(9) .gt. 0.0D0 .and. lady(9) .eq. 3) then
            l = 9
            do 39 k = ngrs, ngre - 1
             cs     (k,l   ) =   frn(k+1) / frn(k) - frn(k) / frn(k+1)
            dcsdlx00(k,l,if) = - frn(k+1) / frn(k) - frn(k) / frn(k+1)
            dcsdlxp1(k,l,if) =   frn(k+1) / frn(k) + frn(k) / frn(k+1)
   39       continue
      end if
c
c.......................................................................
c     luminosity or velocity
c.......................................................................
c
      if (wt(10) .gt. 0.0D0 .and. lady(10) .eq. 3) then
            l = 10
            do 40 k = ngrs, ngre - 1
             cs     (k,l   ) =   bri(k+1) / bri(k) - bri(k) / bri(k+1)
            dcsdlx00(k,l,jb) = - bri(k+1) / bri(k) - bri(k) / bri(k+1)
            dcsdlxp1(k,l,jb) =   bri(k+1) / bri(k) + bri(k) / bri(k+1)
   40       continue
      end if
c
c=======================================================================
c     sums and derivatives
c=======================================================================
c
c     equation ag42
c
      do 43  l = 1, mad
      if (wt(l) .gt. 0.0D0) then
            do 50 k  = ngrs, ngre - 1
            ss   (k) = ss(k) + wt(l) * cs(k,l)**2
   50       continue
c
c     equations ag186 - ag188
c
            do 42 m = 1, meqn
            do 41 k = ngrs, ngre - 1
            dssdlx00(k,m) = dssdlx00(k,m) + 2.0D0 *  cs     (k,l) 
     .                                    * wt(l) * dcsdlx00(k,l,m)
            dssdlxp1(k,m) = dssdlxp1(k,m) + 2.0D0 *  cs     (k,l)
     .                                    * wt(l) * dcsdlxp1(k,l,m)
            dssdlxp2(k,m) = dssdlxp2(k,m) + 2.0D0 *  cs     (k,l)
     .                                    * wt(l) * dcsdlxp2(k,l,m)
   41       continue
   42       continue
      end if
   43 continue
c
      do 44 k = ngrs + 1, ngre - 1
c
c     equation ag43
c
      rr(k) = sqrt( 1.0D0 + xnu(k)**2 * ss(k) )
c
c     equations  ag44 - ag47
c
      dabdlrm1(k) =    dncdlrm1(k) / rr(k)
      dabdlr00(k) =    dncdlr00(k) / rr(k) 
     .              - (dnudlr00(k) / rr(k)**3) * xnc(k)*xnu(k)*ss(k)
      dabdlrp1(k) =    dncdlrp1(k) / rr(k) 
     .              - (dnudlrp1(k) / rr(k)**3) * xnc(k)*xnu(k)*ss(k)
      dabdlrp2(k) =    dncdlrp2(k) / rr(k)
   44 continue
c
      do 45 m = 1, meqn
      do 46 k = ngrs + 1, ngre - 1
c
c     equations ag189 - ag191
c
      dcddlx00(k,m) = 0.5D0 * dssdlx00(k,m) * xnc(k)*xnu(k)**2/rr(k)**3
      dcddlxp1(k,m) = 0.5D0 * dssdlxp1(k,m) * xnc(k)*xnu(k)**2/rr(k)**3
      dcddlxp2(k,m) = 0.5D0 * dssdlxp2(k,m) * xnc(k)*xnu(k)**2/rr(k)**3
   46 continue
   45 continue
c
c=======================================================================
c     matrix elements
c=======================================================================
c
      k = ngrs + 1
c
c     equations ag25 - ag28
c
      rhs(ir,     k) = xnu(k) - xnu(k-1)
c
      em1(ir, jr, k) = em1(ir, jr, k) - dnudlr00(k-1)
      e00(ir, jr, k) = e00(ir, jr, k) - dnudlrp1(k-1) + dnudlr00(k)
      ep1(ir, jr, k) = ep1(ir, jr, k)                 + dnudlrp1(k)
c
c     equations ag29 & ag48 - ag52
c
      do 47 k = ngrs + 2, ngre - 1
      rhs(ir,     k) = xnc(k) / rr(k) - xnc(k-1) / rr(k-1)
c
      em2(ir, jr, k) = em2(ir, jr, k) - dabdlrm1(k-1)
      em1(ir, jr, k) = em1(ir, jr, k) - dabdlr00(k-1) + dabdlrm1(k)
      e00(ir, jr, k) = e00(ir, jr, k) - dabdlrp1(k-1) + dabdlr00(k) 
      ep1(ir, jr, k) = ep1(ir, jr, k) - dabdlrp2(k-1) + dabdlrp1(k) 
      ep2(ir, jr, k) = ep2(ir, jr, k)                 + dabdlrp2(k)
   47 continue
c
c     equations ag192 - ag215
c
      do 48 m = 1, meqn
      do 49 k = ngrs + 2, ngre - 1
      em1(ir, m, k) = em1(ir, m, k) + dcddlx00(k-1, m)
      e00(ir, m, k) = e00(ir, m, k) + dcddlxp1(k-1, m) - dcddlx00(k, m)
      ep1(ir, m, k) = ep1(ir, m, k) + dcddlxp2(k-1, m) - dcddlxp1(k, m)
      ep2(ir, m, k) = ep2(ir, m, k)                    - dcddlxp2(k, m)
   49 continue
   48 continue
c
      k = ngre 
c
c     equations ag25 - ag28
c
      rhs(ir,     k) = xnu(k) - xnu(k-1)
c
      em1(ir, jr, k) = em1(ir, jr, k) - dnudlr00(k-1)
      e00(ir, jr, k) = e00(ir, jr, k) - dnudlrp1(k-1) + dnudlr00(k)
      ep1(ir, jr, k) = ep1(ir, jr, k)                 + dnudlrp1(k)
c
c=======================================================================
c     adaptive grid boundary conditions
c=======================================================================
c
      k = ngrs 
c
c-----------------------------------------------------------------------
c     eulerian inner boundary; equation bc13
c-----------------------------------------------------------------------
c
      if (leibc .gt. 0) then
            rhs(ir,     k) = rn(k) - ro(k)
            e00(ir, jr, k) = rn(k)
      end if
c
c-----------------------------------------------------------------------
c     lagrangean inner boundary; equation bc64
c-----------------------------------------------------------------------
c
      if (llibc .eq. 1) then
	    rhs(ir,     k) = (rn(k) - ro(k)) / dtime - u   (k)
	    e00(ir, jr, k) =  rn(k)          / dtime
            e00(ir, ju, k) =                         - unom(k) * thet
      end if
c
      k = ngre + 1
c
c-----------------------------------------------------------------------
c     eulerian outer boundary; equation bc26
c-----------------------------------------------------------------------
c
      if (leobc .gt. 0) then
            rhs(ir,     k) = rn(k) - ro(k)
            e00(ir, jr, k) = rn(k)
      end if
c
c-----------------------------------------------------------------------
c     lagrangean outer boundary; equation bc67
c-----------------------------------------------------------------------
c
      if (llobc .gt. 0) then
	    rhs(ir,     k) = (rn(k) - ro(k)) / dtime - u   (k)
	    e00(ir, jr, k) =  rn(k)          / dtime
            e00(ir, ju, k) =                         - unom(k) * thet
      end if
c
c+++++++++ lgrid .eq. 1 ++++++++++++++++++++++++++++++++++++++++++++ end
c
      end if
c
c=======================================================================
c     general boundary conditions  
c=======================================================================
c
      k = ngrs - 1
c
c-----------------------------------------------------------------------
c     eulerian inner boundary
c-----------------------------------------------------------------------
c
c     zero flux; equation pz1 & pz2
c
      if (leibc .eq. 1) then
           rhs(ir,     k) = rn(k) - 2.0D0 * rn(k+1) + rn(k+2)
           ep2(ir, jr, k) =                           rn(k+2)
           ep1(ir, jr, k) =       - 2.0D0 * rn(k+1)
           e00(ir, jr, k) = rn(k)
      end if
c
c     nonzero flux; equation pz11 & pz12
c
      if (leibc .eq. 2) then
	    rhs(ir,     k) = rn(k) - rn(k+1) + delrl
	    ep1(ir, jr, k) =       - rn(k+1)
	    e00(ir, jr, k) = rn(k)
      end if
c
c-----------------------------------------------------------------------
c     lagrangean inner boundary 
c-----------------------------------------------------------------------
c
c     equation pz21 & pz22
c
      if (llibc .eq. 1) then
        rhs(ir,     k) =
     .                 - delml/dn(k) + (rmup1n(k+1) - rmup1n(k)) / xmup1
        e00(ir, jr, k) =                            - rmup1n(k)
        ep1(ir, jr, k) =                rmup1n(k+1)
        e00(ir, jd, k) = delml/dn(k)
      end if
c
      k = ngre + 2
c
c-----------------------------------------------------------------------
c     eulerian outer boundary
c-----------------------------------------------------------------------
c
c     zero flux; equations pz39 & pz40
c
      if (leobc .eq. 1) then
            rhs(ir,     k) = rn(k) - 2.0D0 * rn(k-1) + rn(k-2)
            e00(ir, jr, k) = rn(k)
            em1(ir, jr, k) =       - 2.0D0 * rn(k-1)
            em2(ir, jr, k) =                           rn(k-2)
      end if
c
c     nonzero flux; equations pz49 & pz50
c
      if (leobc .eq. 2) then
	    rhs(ir,     k) = rn(k) - rn(k-1) - delrr
	    e00(ir, jr, k) = rn(k)
	    em1(ir, jr, k) =       - rn(k-1)
      end if
c
c     transmitting; equations pz59 & pz60
c
      if (leobc .eq. 3) then
	    rhs(ir,     k) = rn(k) - rn(k-1) - delrr
	    e00(ir, jr, k) = rn(k)
	    em1(ir, jr, k) =       - rn(k-1)
      end if
c
c-----------------------------------------------------------------------
c     lagrangean outer boundary
c-----------------------------------------------------------------------
c
c     equation pz69 & pz70
c
      if (llobc .gt. 0) then
	rhs(ir,     k) = 
     .                  -delmr/dn(k-1) + (rmup1n(k) - rmup1n(k-1))/xmup1
	em1(ir, jr, k) =                            - rmup1n(k-1)
	e00(ir, jr, k) =                  rmup1n(k)
	em1(ir, jd, k) = delmr/dn(k-1)
      end if
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 60 k = 1, ngrs - 2
      rhs(ir,     k) = 0.0D0
      e00(ir, jr, k) = 1.0D0
      ep2(ir, jr, k) = 0.0D0
      ep1(ir, jr, k) = 0.0D0
      em1(ir, jr, k) = 0.0D0
      em2(ir, jr, k) = 0.0D0
   60 continue
c
      do 61 k = ngre + 3, mgr
      rhs(ir,     k) = 0.0D0
      e00(ir, jr, k) = 1.0D0
      ep2(ir, jr, k) = 0.0D0
      ep1(ir, jr, k) = 0.0D0
      em1(ir, jr, k) = 0.0D0
      em2(ir, jr, k) = 0.0D0
   61 continue
c
c-----------------------------------------------------------------------
c
      return
      end
