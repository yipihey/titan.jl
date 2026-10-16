      subroutine grid3nit (ydl, ytl, ydr, ytr, radc, rw)
c
c***********************************************************************
c     grid initialization for SOD's shock tube
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      dimension dg(mgr), tg(mgr)
      dimension dlddlr00(mgr), dlddlrp1(mgr)
      dimension dltdlr00(mgr), dltdlrp1(mgr)
c
c=======================================================================
c     define functions
c=======================================================================
c
      f0 (x) = ( exp(x) - exp(-x) ) / ( exp(-x) + exp(x) ) 
      f1 (x) = + 2.0D0              / ( exp(-x) + exp(x) )
c
c=======================================================================
c     modify original code parameters 
c=======================================================================
c
      iro = ir
      imo = im
      ido = id
      iuo = iu
      ito = it
      jro = jr
      jmo = jm
      jdo = jd
      juo = ju
      jto = jt
c
      ir = 1
      id = 2
      it = 3
      jr = ir
      jd = id
      jt = it
c
      lady2   = lady(2)
      lady(2) = 2
      lady5   = lady(5)
      lady(5) = 2
      lady6   = lady(6)
      lady(6) = 2
      lady9   = lady(9)
      lady(9) = 0
c
      wt2   = wt(2)
      wt(2) = 1.0D0
      wt5   = wt(5)
      wt(5) = 1.0D0
      wt6   = wt(6)
      wt(6) = 1.0D0
c
      neqo  = neqn
      neqn  = 3
      tauo  = tau
      tau   = 1.D-2
c
      dtime = tau
c
c=======================================================================
c     NEWTON-RAPHSON iteration
c=======================================================================
c
      iter = 1
c
c.......................................................................
c     calculate derivatives 
c.......................................................................
c
        dg    (ngrs-2) = ydl
      dlddlr00(ngrs-2) = 0.0D0
      dlddlrp1(ngrs-2) = 0.0D0
        tg    (ngrs-2) = ytl
      dltdlr00(ngrs-2) = 0.0D0
      dltdlrp1(ngrs-2) = 0.0D0
      do 1 k = ngrs - 1, ngre + 1
c
      rmid   = ( 0.5D0 * ( rn(k) + rn(k+1) ) - radc ) / rw
      rmid   = max( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(0.5D0, 0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        tg    (k) = ytr + 0.5D0*( 1.0D0 - f0 (rmid) )*( ytl - ytr )
      dltdlr00(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
      dltdlrp1(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
c
        dg    (k) = ydr + 0.5D0*( 1.0D0 - f0 (rmid) )*( ydl - ydr )
      dlddlr00(k) = - switch*( ydl - ydr ) * 0.5D0/rw * f1(rmid)**2
      dlddlrp1(k) = - switch*( ydl - ydr ) * 0.5D0/rw * f1(rmid)**2
c
    1 continue
        dg    (ngre+2) = ydr
      dlddlr00(ngre+2) = 0.0D0
      dlddlrp1(ngre+2) = 0.0D0
        tg    (ngre+2) = ytr
      dltdlr00(ngre+2) = 0.0D0
      dltdlrp1(ngre+2) = 0.0D0
c
      do 2 k = ngrs - 2, ngre + 2
      dn  (k) = dg(k)
      tn  (k) = tg(k)
      ro  (k) = rn(k)
    2 continue
c
c.......................................................................
c     initialize all matrices to zero
c.......................................................................
c
  100 continue
      do 30 i = 1, meqn
      do 20 j = 1, meqn
      do 10 k = 1, mgr
      em2(i, j, k) = 0.0D0
      em1(i, j, k) = 0.0D0
      e00(i, j, k) = 0.0D0
      ep1(i, j, k) = 0.0D0
      ep2(i, j, k) = 0.0D0
   10 continue
   20 continue
   30 continue
c
      do 50 i = 1, meqn
      do 40 k = 1, mgr
      rhs(i, k) = 0.0D0
   40 continue
   50 continue
c
      do 70 i = 1, 3*mpd - 2
      do 60 j = 1, meqn*mgr
      bm(i, j) = 0.0D0
   60 continue
   70 continue
c
      do 80 i = 1, meqn*mgr
      br(i) = 0.0D0
   80 continue
c
c.......................................................................
c     calculate matrix elements in continuity and gas energy equation
c.......................................................................
c
      e00(id, jd, ngrs-2) =  1.0D0
      e00(it, jt, ngrs-2) =  1.0D0
c
      do 3 k = ngrs - 1, ngre + 1
c
      rhs(id,     k) =  dn(k) - dg(k)
      e00(id, jd, k) =  dn(k) 
      e00(id, jr, k) =        - dlddlr00(k) * rn(k  )
      ep1(id, jr, k) =        - dlddlrp1(k) * rn(k+1)
c
      rhs(it,     k) =  tn(k) - tg(k)
      e00(it, jt, k) =  tn(k) 
      e00(it, jr, k) =        - dltdlr00(k) * rn(k  )
      ep1(it, jr, k) =        - dltdlrp1(k) * rn(k+1)
    3 continue
c
      e00(id, jd, ngre+2) =  1.0D0
      e00(it, jt, ngre+2) =  1.0D0
c
c.......................................................................
c     calculate matrix elements in grid equation
c.......................................................................
c
      call grid
c
c.......................................................................
c     solve the system
c.......................................................................
c
      call matslv
c
c.......................................................................
c     find maximum changes 
c.......................................................................
c
      kdmx = isamax(ngre - ngrs + 2, rhs(ir,ngrs), meqn)
      dmax = abs( rhs(ir,ngrs + kdmx - 1) ) + 1.0D-30
c
      do 4 k = ngrs, ngre
      cx  (k) = ( rhs(ir, k+1) * rn(k+1) - rhs(ir, k) * rn(k) )
     .        / (                rn(k+1) -              rn(k) )
    4 continue
      kcmx = isamax( ngre - ngrs + 1, cx(ngrs), 1 )
      cmax = abs( cx(ngrs + kcmx - 1) ) + 1.0D-30
c
      sc1 = min( 1.0D0, dtol / dmax )
      sc2 = min( 1.0D0, ctol / cmax )
      sc  = min( sc1, sc2 )
c
      do 5 k = ngrs, ngre + 1
      ro  (k) =   rn(k)
      rn  (k) = ( 1.0D0 + sc * rhs(ir,k) ) * ro(k)
      cx  (k) = ( rn(k) - ro(k) ) / ro(k)
    5 continue
      ksmx = isamax(ngre - ngrs + 2, cx(ngrs), 1)
      smax = abs( cx(ngrs + ksmx - 1) ) + 1.0D-30
c
      iter = iter + 1
c.......................................................................
c     summary output
c.......................................................................
c
      write(itty,998) iter,timen,ksmx,smax,kdmx,dmax
  998 format(" iter ="i3" time ="1pe10.3
     .       " d_r/ro("i3") ="e9.2" rhs("i3") ="e9.2 )
c
c-----------------------------------------------------------------------
c     check for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre + 1
      if ( rn(k+1) .le. rn(k) ) then
           write(itty,'(" grid extrapolation nonmonotonic")')
           write(iout,'(" grid extrapolation nonmonotonic")')
           call clzfil
           stop 'gridinit03'
           end if
    6 continue
c
c-----------------------------------------------------------------------
c     compute new timestep and update quantities
c-----------------------------------------------------------------------
c
      timeo = timen
      timen = timeo + dtime
c
        dg    (ngrs-2) = ydl
      dlddlr00(ngrs-2) = 0.0D0
      dlddlrp1(ngrs-2) = 0.0D0
        tg    (ngrs-2) = ytl
      dltdlr00(ngrs-2) = 0.0D0
      dltdlrp1(ngrs-2) = 0.0D0
      do 7 k = ngrs - 1, ngre + 1
c
      rmid   = ( 0.5D0 * ( rn(k) + rn(k+1) ) - radc ) / rw
      rmid   = max( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(0.5D0, 0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        tg    (k) = ytr + 0.5D0*( 1.0D0 - f0 (rmid) )*( ytl - ytr )
      dltdlr00(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
      dltdlrp1(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
c
        dg    (k) = ydr + 0.5D0*( 1.0D0 - f0 (rmid) )*( ydl - ydr )
      dlddlr00(k) = - switch*( ydl - ydr ) * 0.5D0/rw * f1(rmid)**2
      dlddlrp1(k) = - switch*( ydl - ydr ) * 0.5D0/rw * f1(rmid)**2
c
    7 continue
        dg    (ngre+2) = ydr
      dlddlr00(ngre+2) = 0.0D0
      dlddlrp1(ngre+2) = 0.0D0
        tg    (ngre+2) = ytr
      dltdlr00(ngre+2) = 0.0D0
      dltdlrp1(ngre+2) = 0.0D0
c
      do 8 k = ngrs - 2, ngre + 2
      dn  (k) = dg(k)
      tn  (k) = tg(k)
    8 continue
c
      call eos
c
      if (smax .gt. 1.D-8) go to 100
c
c=======================================================================
c     restore original code parameters 
c=======================================================================
c
      ir = iro
      im = imo
      id = ido
      iu = iuo
      it = ito
      jr = jro
      jm = jmo
      jd = jdo
      ju = juo
      jt = jto
c
      neqn  = neqo
      tau   = tauo
c
      lady(2) = lady2
      lady(5) = lady5
      lady(6) = lady6
      lady(9) = lady9
c
      wt(2) = wt2
      wt(5) = wt5
      wt(6) = wt6
c
c=======================================================================
c
      return
      end
      subroutine grid4nit (ydl1,ydl2,ytl1,ytl2,
     .                     ydr1,ydr2,ytr1,ytr2,x1,x2,xw)
c
c***********************************************************************
c     grid initialization for WOODWARD-COLLELA blast wave
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      dimension dg(mgr), tg(mgr)
      dimension dlddlr00(mgr), dlddlrp1(mgr)
      dimension dltdlr00(mgr), dltdlrp1(mgr)
c
c=======================================================================
c     define functions
c=======================================================================
c
      f0 (x) = ( exp(x) - exp(-x) ) / ( exp(-x) + exp(x) ) 
      f1 (x) = + 2.0D0              / ( exp(-x) + exp(x) )
c
c=======================================================================
c     modify original code parameters 
c=======================================================================
c
      iro = ir
      imo = im
      ido = id
      iuo = iu
      ito = it
      jro = jr
      jmo = jm
      jdo = jd
      juo = ju
      jto = jt
c
      ir = 1
      id = 2
      it = 3
      jr = ir
      jd = id
      jt = it
c
      neqo  = neqn
      neqn  = 3
      tauo  = tau
      tau   = 1.D-2
c
      xw    = 1.3D-3
      x2    = 1.9D0
      x1    = 1.1D0
      dtime = tau
c
c=======================================================================
c     NEWTON-RAPHSON iteration
c=======================================================================
c
      iter = 1
c
c.......................................................................
c     calculate derivatives 
c.......................................................................
c
        dg    (ngrs-2) = 1.0D0
      dlddlr00(ngrs-2) = 0.0D0
      dlddlrp1(ngrs-2) = 0.0D0
        tg    (ngrs-2) = ytl1
      dltdlr00(ngrs-2) = 0.0D0
      dltdlrp1(ngrs-2) = 0.0D0
        dg    (ngrs-1) = 1.0D0
      dlddlr00(ngrs-1) = 0.0D0
      dlddlrp1(ngrs-1) = 0.0D0
        tg    (ngrs-1) = ytl1
      dltdlr00(ngrs-1) = 0.0D0
      dltdlrp1(ngrs-1) = 0.0D0
      do 1 k = ngrs, ngre
c
      rmid1   =      (  0.5D0 * ( rn(k) + rn(k+1) ) - x1 ) / xw
      rmid1   =   max( -3.D+2, min( rmid1, +3.D+2 ) )      ! 6.90775D+2
      switch1 = cvmgt(  0.5D0, 0.0D0, abs(f1 (rmid1)) .lt. 1.D+150)
      rmid2   =      (  0.5D0 * ( rn(k) + rn(k+1) ) - x2 ) / xw
      rmid2   =   max( -3.D+2, min( rmid2, +3.D+2 ) )      ! 6.90775D+2
      switch2 = cvmgt(  0.5D0, 0.0D0, abs(f1 (rmid2)) .lt. 1.D+150)
c
      dg      (k) = 1.0D0
      dlddlr00(k) = 0.0D0
      dlddlrp1(k) = 0.0D0
c
      tg      (k) = ytr1 + 0.5D0*( 1.0D0 - f0 (rmid1) )*( ytl1 - ytr1 )
     .                   - 0.5D0*( 1.0D0 + f0 (rmid2) )*( ytl2 - ytr2 )
      dltdlr00(k) = - switch1*( ytl1 - ytr1 ) * 0.5D0/xw * f1(rmid1)**2
     .              - switch2*( ytl2 - ytr2 ) * 0.5D0/xw * f1(rmid2)**2
      dltdlrp1(k) = - switch1*( ytl1 - ytr1 ) * 0.5D0/xw * f1(rmid1)**2
     .              - switch2*( ytl2 - ytr2 ) * 0.5D0/xw * f1(rmid2)**2
c
    1 continue
        dg    (ngre+1) = 1.0D0
      dlddlr00(ngre+1) = 0.0D0
      dlddlrp1(ngre+1) = 0.0D0
        tg    (ngre+1) = ytr2
      dltdlr00(ngre+1) = 0.0D0
      dltdlrp1(ngre+1) = 0.0D0
        dg    (ngre+2) = 1.0D0
      dlddlr00(ngre+2) = 0.0D0
      dlddlrp1(ngre+2) = 0.0D0
        tg    (ngre+2) = ytr2
      dltdlr00(ngre+2) = 0.0D0
      dltdlrp1(ngre+2) = 0.0D0
c
      do 2 k = ngrs - 2, ngre + 2
      dn  (k) = dg(k)
      tn  (k) = tg(k)
      ro  (k) = rn(k)
    2 continue
c
c.......................................................................
c     initialize all matrices to zero
c.......................................................................
c
  100 continue
      do 30 i = 1, meqn
      do 20 j = 1, meqn
      do 10 k = 1, mgr
      em2(i, j, k) = 0.0D0
      em1(i, j, k) = 0.0D0
      e00(i, j, k) = 0.0D0
      ep1(i, j, k) = 0.0D0
      ep2(i, j, k) = 0.0D0
   10 continue
   20 continue
   30 continue
c
      do 50 i = 1, meqn
      do 40 k = 1, mgr
      rhs(i, k) = 0.0D0
   40 continue
   50 continue
c
      do 70 i = 1, 3*mpd - 2
      do 60 j = 1, meqn*mgr
      bm(i, j) = 0.0D0
   60 continue
   70 continue
c
      do 80 i = 1, meqn*mgr
      br(i) = 0.0D0
   80 continue
c
c.......................................................................
c     calculate matrix elements in continuity and gas energy equation
c.......................................................................
c
      e00(id, jd, ngrs-2) =  1.0D0
      e00(it, jt, ngrs-2) =  1.0D0
c
      do 3 k = ngrs - 1, ngre + 1
c
      rhs(id,     k) =  dn(k) - dg      (k)
      e00(id, jd, k) =  dn(k) 
      e00(id, jr, k) =        - dlddlr00(k) * rn(k  )
      ep1(id, jr, k) =        - dlddlrp1(k) * rn(k+1)
c
      rhs(it,     k) =  tn(k) - tg      (k)
      e00(it, jt, k) =  tn(k) 
      e00(it, jr, k) =        - dltdlr00(k) * rn(k  )
      ep1(it, jr, k) =        - dltdlrp1(k) * rn(k+1)
    3 continue
c
      e00(id, jd, ngre+2) =  1.0D0
      e00(it, jt, ngre+2) =  1.0D0
c
c.......................................................................
c     calculate matrix elements in grid equation
c.......................................................................
c
      call grid
c
c.......................................................................
c     solve the system
c.......................................................................
c
      call matslv
c
c.......................................................................
c     find maximum changes 
c.......................................................................
c
      kdmx = isamax(ngre - ngrs + 2, rhs(ir,ngrs), meqn)
      dmax = abs( rhs(ir,ngrs + kdmx - 1) ) + 1.0D-30
c
      do 4 k = ngrs, ngre
      cx  (k) = ( rhs(ir, k+1) * rn(k+1) - rhs(ir, k) * rn(k) )
     .        / (                rn(k+1) -              rn(k) )
    4 continue
      kcmx = isamax( ngre - ngrs + 1, cx(ngrs), 1 )
      cmax = abs( cx(ngrs + kcmx - 1) ) + 1.0D-30
c
      sc1 = min( 1.0D0, dtol / dmax )
      sc2 = min( 1.0D0, ctol / cmax )
      sc  = min( sc1, sc2 )
c
      do 5 k = ngrs, ngre + 1
      ro  (k) =   rn(k)
      rn  (k) = ( 1.0D0 + sc * rhs(ir,k) ) * ro(k)
      cx  (k) = ( rn(k) - ro(k) ) / ro(k)
    5 continue
      ksmx = isamax(ngre - ngrs + 2, cx(ngrs), 1)
      smax = abs( cx(ngrs + ksmx - 1) ) + 1.0D-30
c
      iter = iter + 1
c.......................................................................
c     summary output
c.......................................................................
c
      write(itty,998) iter,timen,ksmx,smax,kdmx,dmax
  998 format(" iter ="i4" time ="1pe10.3
     .       " d_r/ro("i3") ="e9.2" rhs("i3") ="e9.2 )
c
c-----------------------------------------------------------------------
c     check for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre + 1
      if ( rn(k+1) .le. rn(k) ) then
           write(itty,'(" grid extrapolation nonmonotonic")')
           write(iout,'(" grid extrapolation nonmonotonic")')
           call clzfil
           stop 'gridinit04'
           end if
    6 continue
c
c-----------------------------------------------------------------------
c     compute new timestep
c-----------------------------------------------------------------------
c
      timeo = timen
      timen = timeo + dtime
c
        dg    (ngrs-2) = 1.0D0
      dlddlr00(ngrs-2) = 0.0D0
      dlddlrp1(ngrs-2) = 0.0D0
        tg    (ngrs-2) = ytl1
      dltdlr00(ngrs-2) = 0.0D0
      dltdlrp1(ngrs-2) = 0.0D0
        dg    (ngrs-1) = 1.0D0
      dlddlr00(ngrs-1) = 0.0D0
      dlddlrp1(ngrs-1) = 0.0D0
        tg    (ngrs-1) = ytl1
      dltdlr00(ngrs-1) = 0.0D0
      dltdlrp1(ngrs-1) = 0.0D0
      do 7 k = ngrs, ngre
c
      rmid1   =      (  0.5D0 * ( rn(k) + rn(k+1) ) - x1 ) / xw
      rmid1   =   max( -3.D+2, min( rmid1, +3.D+2 ) )      ! 6.90775D+2
      switch1 = cvmgt(  0.5D0, 0.0D0, abs(f1 (rmid1)) .lt. 1.D+150)
      rmid2   =      (  0.5D0 * ( rn(k) + rn(k+1) ) - x2 ) / xw
      rmid2   =   max( -3.D+2, min( rmid2, +3.D+2 ) )      ! 6.90775D+2
      switch2 = cvmgt(  0.5D0, 0.0D0, abs(f1 (rmid2)) .lt. 1.D+150)
c
      dg      (k) = 1.0D0
      dlddlr00(k) = 0.0D0
      dlddlrp1(k) = 0.0D0
c
      tg      (k) = ytr1 + 0.5D0*( 1.0D0 - f0 (rmid1) )*( ytl1 - ytr1 )
     .                   - 0.5D0*( 1.0D0 + f0 (rmid2) )*( ytl2 - ytr2 )
      dltdlr00(k) = - switch1*( ytl1 - ytr1 ) * 0.5D0/xw * f1(rmid1)**2
     .              - switch2*( ytl2 - ytr2 ) * 0.5D0/xw * f1(rmid2)**2
      dltdlrp1(k) = - switch1*( ytl1 - ytr1 ) * 0.5D0/xw * f1(rmid1)**2
     .              - switch2*( ytl2 - ytr2 ) * 0.5D0/xw * f1(rmid2)**2
c
    7 continue
        dg    (ngre+1) = 1.0D0
      dlddlr00(ngre+1) = 0.0D0
      dlddlrp1(ngre+1) = 0.0D0
        tg    (ngre+1) = ytr2
      dltdlr00(ngre+1) = 0.0D0
      dltdlrp1(ngre+1) = 0.0D0
        dg    (ngre+2) = 1.0D0
      dlddlr00(ngre+2) = 0.0D0
      dlddlrp1(ngre+2) = 0.0D0
        tg    (ngre+2) = ytr2
      dltdlr00(ngre+2) = 0.0D0
      dltdlrp1(ngre+2) = 0.0D0
c
      do 8 k = ngrs - 2, ngre + 2
      dn  (k) = dg(k)
      tn  (k) = tg(k)
    8 continue
c
      call eos
c
      if (smax .gt. 1.D-7) go to 100
c
c=======================================================================
c     restore original code parameters 
c=======================================================================
c
      ir = iro
      im = imo
      id = ido
      iu = iuo
      it = ito
      jr = jro
      jm = jmo
      jd = jdo
      ju = juo
      jt = jto
c
      neqn  = neqo
      tau   = tauo
c
c=======================================================================
c
      return
      end
      subroutine grid5nit (ytl, ytr, radc, rw)
c
c***********************************************************************
c     grid initialization for SEDOV-TAYLOR blast wave
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      dimension tg(mgr)
      dimension dltdlr00(mgr), dltdlrp1(mgr)
c
c=======================================================================
c     define functions
c=======================================================================
c
      f0 (x) = ( exp(x) - exp(-x) ) / ( exp(-x) + exp(x) ) 
      f1 (x) = + 2.0D0              / ( exp(-x) + exp(x) )
c
c=======================================================================
c     modify original code parameters 
c=======================================================================
c
      iro = ir
      imo = im
      ido = id
      iuo = iu
      ito = it
      jro = jr
      jmo = jm
      jdo = jd
      juo = ju
      jto = jt
c
      ir = 1
      it = 2
      jr = ir
      jt = it
c
      lady1   = lady(1)
      lady(1) = 0
      lady2   = lady(2)
      lady(2) = 0
      lady5   = lady(5)
      lady(5) = 2
      lady6   = lady(6)
      lady(6) = 2
      lady9   = lady(9)
      lady(9) = 0
c
      wt2   = wt(2)
      wt(2) = 0.0D0
      wt5   = wt(5)
      wt(5) = 1.0D0
      wt6   = wt(6)
      wt(6) = 1.0D0
c
      neqo  = neqn
      neqn  = 2
      tauo  = tau
      tau   = 1.D-2
c
      dtime = tau
c
c=======================================================================
c     NEWTON-RAPHSON iteration
c=======================================================================
c
      iter = 1
c
c.......................................................................
c     calculate derivatives 
c.......................................................................
c
        tg    (ngrs-2) = ytl
      dltdlr00(ngrs-2) = 0.0D0
      dltdlrp1(ngrs-2) = 0.0D0
        tg    (ngrs-1) = ytl
      dltdlr00(ngrs-1) = 0.0D0
      dltdlrp1(ngrs-1) = 0.0D0
        tg    (ngrs  ) = ytl
      dltdlr00(ngrs  ) = 0.0D0
      dltdlrp1(ngrs  ) = 0.0D0
      do 1 k = ngrs + 1, ngre + 1
c
      rmid   = ( 0.5D0 * ( rn(k) + rn(k+1) ) - radc ) / rw
      rmid   = max( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(0.5D0, 0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        tg    (k) = ytr + 0.5D0*( 1.0D0 - f0 (rmid) )*( ytl - ytr )
      dltdlr00(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
      dltdlrp1(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
c
    1 continue
        tg    (ngre+2) = ytr
      dltdlr00(ngre+2) = 0.0D0
      dltdlrp1(ngre+2) = 0.0D0
c
      do 2 k = ngrs - 2, ngre + 2
      tn  (k) = tg(k)
      ro  (k) = rn(k)
    2 continue
c
c.......................................................................
c     initialize all matrices to zero
c.......................................................................
c
  100 continue
      do 30 i = 1, meqn
      do 20 j = 1, meqn
      do 10 k = 1, mgr
      em2(i, j, k) = 0.0D0
      em1(i, j, k) = 0.0D0
      e00(i, j, k) = 0.0D0
      ep1(i, j, k) = 0.0D0
      ep2(i, j, k) = 0.0D0
   10 continue
   20 continue
   30 continue
c
      do 50 i = 1, meqn
      do 40 k = 1, mgr
      rhs(i, k) = 0.0D0
   40 continue
   50 continue
c
      do 70 i = 1, 3*mpd - 2
      do 60 j = 1, meqn*mgr
      bm(i, j) = 0.0D0
   60 continue
   70 continue
c
      do 80 i = 1, meqn*mgr
      br(i) = 0.0D0
   80 continue
c
c.......................................................................
c     calculate matrix elements in continuity and gas energy equation
c.......................................................................
c
      e00(it, jt, ngrs-2) =  1.0D0
c
      do 3 k = ngrs - 1, ngre + 1
c
      rhs(it,     k) =  tn(k) - tg(k)
      e00(it, jt, k) =  tn(k) 
      e00(it, jr, k) =        - dltdlr00(k) * rn(k  )
      ep1(it, jr, k) =        - dltdlrp1(k) * rn(k+1)
    3 continue
c
      e00(it, jt, ngre+2) =  1.0D0
c
c.......................................................................
c     calculate matrix elements in grid equation
c.......................................................................
c
      call grid
c
c.......................................................................
c     solve the system
c.......................................................................
c
      call matslv
c
c.......................................................................
c     find maximum changes 
c.......................................................................
c
      kdmx = isamax(ngre - ngrs + 2, rhs(ir,ngrs), meqn)
      dmax = abs( rhs(ir,ngrs + kdmx - 1) ) + 1.0D-30
c
      do 4 k = ngrs, ngre
      cx  (k) = ( rhs(ir, k+1) * rn(k+1) - rhs(ir, k) * rn(k) )
     .        / (                rn(k+1) -              rn(k) )
    4 continue
      kcmx = isamax( ngre - ngrs + 1, cx(ngrs), 1 )
      cmax = abs( cx(ngrs + kcmx - 1) ) + 1.0D-30
c
      sc1 = min( 1.0D0, dtol / dmax )
      sc2 = min( 1.0D0, ctol / cmax )
      sc  = min( sc1, sc2 )
c
      do 5 k = ngrs, ngre + 1
      ro  (k) =   rn(k)
      rn  (k) = ( 1.0D0 + sc * rhs(ir,k) ) * ro(k)
      cx  (k) = ( rn(k) - ro(k) ) / ro(k)
    5 continue
      ksmx = isamax(ngre - ngrs + 2, cx(ngrs), 1)
      smax = abs( cx(ngrs + ksmx - 1) ) + 1.0D-30
c
      iter = iter + 1
c.......................................................................
c     summary output
c.......................................................................
c
      write(itty,998) iter,timen,ksmx,smax,kdmx,dmax
  998 format(" iter ="i4" time ="1pe10.3
     .       " d_r/ro("i3") ="e9.2" rhs("i3") ="e9.2 )
c
c-----------------------------------------------------------------------
c     check for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre + 1
      if ( rn(k+1) .le. rn(k) ) then
           write(itty,'(" grid extrapolation nonmonotonic")')
           write(iout,'(" grid extrapolation nonmonotonic")')
           call clzfil
           stop 'gridinit05'
           end if
    6 continue
c
c-----------------------------------------------------------------------
c     compute new timestep and update quantities
c-----------------------------------------------------------------------
c
      timeo = timen
      timen = timeo + dtime
c
        tg    (ngrs-2) = ytl
      dltdlr00(ngrs-2) = 0.0D0
      dltdlrp1(ngrs-2) = 0.0D0
        tg    (ngrs-1) = ytl
      dltdlr00(ngrs-1) = 0.0D0
      dltdlrp1(ngrs-1) = 0.0D0
        tg    (ngrs  ) = ytl
      dltdlr00(ngrs  ) = 0.0D0
      dltdlrp1(ngrs  ) = 0.0D0
      do 7 k = ngrs + 1, ngre + 1
c
      rmid   = ( 0.5D0 * ( rn(k) + rn(k+1) ) - radc ) / rw
      rmid   = max( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(0.5D0, 0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        tg    (k) = ytr + 0.5D0*( 1.0D0 - f0 (rmid) )*( ytl - ytr )
      dltdlr00(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
      dltdlrp1(k) = - switch*( ytl - ytr ) * 0.5D0/rw * f1(rmid)**2
c
    7 continue
        tg    (ngre+2) = ytr
      dltdlr00(ngre+2) = 0.0D0
      dltdlrp1(ngre+2) = 0.0D0
c
      do 8 k = ngrs - 2, ngre + 2
      tn  (k) = tg(k)
    8 continue
c
      call eos
c
      if (smax .gt. 1.D-5) go to 100
c
c=======================================================================
c     restore original code parameters 
c=======================================================================
c
      ir = iro
      im = imo
      id = ido
      iu = iuo
      it = ito
      jr = jro
      jm = jmo
      jd = jdo
      ju = juo
      jt = jto
c
      neqn  = neqo
      tau   = tauo
c
      lady(1) = lady1
      lady(2) = lady2
      lady(5) = lady5
      lady(6) = lady6
      lady(9) = lady9
c
      wt(2) = wt2
      wt(5) = wt5
      wt(6) = wt6
c
c=======================================================================
c
      return
      end
      subroutine grid6nit (yfl, yfr, radc, rw)
c
c***********************************************************************
c     grid initialization for radiative equilibrium
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      dimension  tg     (mgr)
      dimension dltdlr00(mgr), dledlr00(mgr), dledlrp1(mgr)
c
c=======================================================================
c     define functions
c=======================================================================
c
      f0 (x) = ( exp(x) - exp(-x) ) / ( exp(-x) + exp(x) ) 
      f1 (x) = + 2.0D0              / ( exp(-x) + exp(x) )
c
c=======================================================================
c     modify original code parameters 
c=======================================================================
c
      iro = ir
      ito = it
      ieo = ie
      ifo = if
      iuo = iu
      jro = jr
      jto = jt
      jeo = je
      jfo = jf
      juo = ju
      jbo = jb
c
      ir = 1
      if = 2
      ie = 3
      jr = ir
      jf = if
      je = ie
      jb = if
c
      neqo  = neqn
      neqn  = 3
      tauo  = tau
      tau   = 1.D-2
c
      dtime = tau
c
      ern0 = yfr / (4.0D0 * cpi * cc * rn(ngre)**2 ) * sqrt(3.0D0)
      con  = yfr / (4.0D0 * cpi * cc)* chif0 * 3.0D0
c
c=======================================================================
c     NEWTON-RAPHSON iteration
c=======================================================================
c
      iter = 1
c
c.......................................................................
c     calculate derivatives 
c.......................................................................
c
      do 1 k = ngrs - 1, ngre + 1
c
      rmid   =      (  rn(k) - radc ) / rw
      rmid   = max  ( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(  0.5D0,  0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        tg    (k) = yfr + 0.5D0*( 1.0D0 - f0 (rmid) )*( yfl - yfr )
      dltdlr00(k) = - switch*( yfl - yfr ) * 0.5D0/rw * f1(rmid)**2
c
      rmid = ( rn(k) + rn(k+1) )
c
        er    (k) =  ern0 + con * (2.D0 / rmid - 1.D0 / rn(ngre+2))
      dledlr00(k) =       - con *  2.D0 / rmid**2
      dledlrp1(k) =       - con *  2.D0 / rmid**2
c
    1 continue
        tg    (ngre+2) = yfr
      dltdlr00(ngre+2) = 0.0D0
        er    (ngre+2) = ern0
      dledlr00(ngre+2) = 0.0D0
        tg    (ngrs-2) = yfl
      dltdlr00(ngrs-2) = 0.0D0
        er    (ngrs-2) = er(ngrs-1)
      dledlr00(ngrs-2) = 0.0D0
c
      do 2 k = ngrs - 1, ngre + 1
      bri (k) = tg(k) 
      ern (k) = er(k)
      ro  (k) = rn(k)
    2 continue
      bri (ngre+2) = bri(ngre+1)
      ern (ngre+2) = ern0
      bri (ngrs-2) = bri(ngrs-1)
      ern (ngrs-2) = ern(ngrs-1)
c
c.......................................................................
c     initialize all matrices to zero
c.......................................................................
c
  100 continue
      do 30 i = 1, meqn
      do 20 j = 1, meqn
      do 10 k = 1, mgr
      em2(i, j, k) = 0.0D0
      em1(i, j, k) = 0.0D0
      e00(i, j, k) = 0.0D0
      ep1(i, j, k) = 0.0D0
      ep2(i, j, k) = 0.0D0
   10 continue
   20 continue
   30 continue
c
      do 50 i = 1, meqn
      do 40 k = 1, mgr
      rhs(i, k) = 0.0D0
   40 continue
   50 continue
c
      do 70 i = 1, 3*mpd - 2
      do 60 j = 1, meqn*mgr
      bm(i, j) = 0.0D0
   60 continue
   70 continue
c
      do 80 i = 1, meqn*mgr
      br(i) = 0.0D0
   80 continue
c
c.......................................................................
c     calculate matrix elements in continuity and gas energy equation
c.......................................................................
c
      e00(if, jf, ngrs-2) =  1.0D0
      e00(ie, je, ngrs-2) =  1.0D0
c
      do 3 k = ngrs - 1, ngre + 1
c
      rhs(if,     k) =  bri(k) - tg(k)
      e00(if, jf, k) =  bri(k) 
      e00(if, jr, k) =         - dltdlr00(k) * rn(k)
c
      rhs(ie,     k) =  ern(k) - er(k)
      e00(ie, je, k) =  ern(k) 
      e00(ie, jr, k) =         - dledlr00(k) * rn(k  )
      ep1(ie, jr, k) =         - dledlrp1(k) * rn(k+1)
    3 continue
c
      e00(if, jf, ngre+2) =  1.0D0
      e00(ie, je, ngre+2) =  1.0D0
c
c.......................................................................
c     calculate matrix elements in grid equation
c.......................................................................
c
      call grid
c
c.......................................................................
c     solve the system
c.......................................................................
c
      call matslv
c
c.......................................................................
c     find maximum changes 
c.......................................................................
c
      kdmx = isamax(ngre - ngrs + 2, rhs(ir,ngrs), meqn)
      dmax = abs( rhs(ir,ngrs + kdmx - 1) ) + 1.0D-30
c
      do 4 k = ngrs, ngre
      cx  (k) = ( rhs(ir, k+1) * rn(k+1) - rhs(ir, k) * rn(k) )
     .        / (                rn(k+1) -              rn(k) )
    4 continue
      kcmx = isamax( ngre - ngrs + 1, cx(ngrs), 1 )
      cmax = abs( cx(ngrs + kcmx - 1) ) + 1.0D-30
c
      sc1 = min( 1.0D0, dtol / dmax )
      sc2 = min( 1.0D0, ctol / cmax )
      sc  = min( sc1, sc2 )
c
      do 5 k = ngrs, ngre + 1
      ro  (k) =   rn(k)
      rn  (k) = ( 1.0D0 + sc * rhs(ir,k) ) * ro(k)
      cx  (k) = ( rn(k) - ro(k) ) / ro(k)
    5 continue
      ksmx = isamax(ngre - ngrs + 2, cx(ngrs), 1)
      smax = abs( cx(ngrs + ksmx - 1) ) + 1.0D-30
c
      iter = iter + 1
c.......................................................................
c     summary output
c.......................................................................
c
      write(itty,998) iter,timen,ksmx,smax,kdmx,dmax
  998 format(" iter ="i4" time ="1pe10.3
     .       " d_r/ro("i3") ="e9.2" rhs("i3") ="e9.2 )
c
c-----------------------------------------------------------------------
c     check for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre + 1
      if ( rn(k+1) .le. rn(k) ) then
           write(itty,'(" grid extrapolation nonmonotonic")')
           write(iout,'(" grid extrapolation nonmonotonic")')
           call clzfil
           stop 'gridinit06'
           end if
    6 continue
c
c-----------------------------------------------------------------------
c     compute new timestep and update quantities
c-----------------------------------------------------------------------
c
      timeo = timen
      timen = timeo + dtime
c
      do 7 k = ngrs - 1, ngre + 1
c
      rmid   =      (  rn(k) - radc ) / rw
      rmid   = max  ( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(  0.5D0,  0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        tg    (k) = yfr + 0.5D0*( 1.0D0 - f0 (rmid) )*( yfl - yfr )
      dltdlr00(k) = - switch*( yfl - yfr ) * 0.5D0/rw * f1(rmid)**2
c
      rmid = ( rn(k) + rn(k+1) )
c
        er    (k) =  ern0 + con * (2.D0 / rmid - 1.D0 / rn(ngre+2))
      dledlr00(k) =       - con *  2.D0 / rmid**2
      dledlrp1(k) =       - con *  2.D0 / rmid**2
c
    7 continue
        tg    (ngre+2) = yfr
      dltdlr00(ngre+2) = 0.0D0
        er    (ngre+2) = ern0
      dledlr00(ngre+2) = 0.0D0
        tg    (ngrs-2) = yfl
      dltdlr00(ngrs-2) = 0.0D0
        er    (ngrs-2) = er(ngrs-1)
      dledlr00(ngrs-2) = 0.0D0
c
      do 8 k = ngrs - 1, ngre + 1
      bri (k) = tg(k)
      ern (k) = er(k)
    8 continue
      bri (ngre+2) = bri(ngre+1)
      ern (ngre+2) = ern0
      bri (ngrs-2) = bri(ngrs-1)
      ern (ngrs-2) = ern(ngrs-1)
c
      if (smax .gt. 1.D-4) go to 100
c
      do 9 k = ngrs - 1, ngre + 1
      frn (k) = tg(k) / ( 4.D0 * cpi * rn(k)**2 ) 
    9 continue
      frn (ngre+2) = frn(ngre+1)
      frn (ngrs-2) = frn(ngrs-1)
c
c=======================================================================
c     restore original code parameters 
c=======================================================================
c
      ir = iro
      it = ito
      ie = ieo
      if = ifo
      iu = iuo
      jr = jro
      jt = jto
      je = jeo
      jf = jfo
      ju = juo
      jb = jbo
c
      neqn  = neqo
      tau   = tauo
c
c=======================================================================
c
      return
      end
      subroutine grid9nit (radc, rw)
c
c***********************************************************************
c     grid initialization for sub- and supercritical radiation
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      dimension   tg    (mgr),   ug    (mgr)
      dimension dltdlr00(mgr), dludlr00(mgr)
c
c=======================================================================
c     define functions
c=======================================================================
c
      f0 (x) = ( exp(x) - exp(-x) ) / ( exp(-x) + exp(x) )
      f1 (x) = + 2.0D0              / ( exp(-x) + exp(x) )
c
c=======================================================================
c     modify original code parameters 
c=======================================================================
c
      iro = ir
      imo = im
      ido = id
      iuo = iu
      ito = it
      ieo = ie
      ifo = if
      jro = jr
      jmo = jm
      jdo = jd
      juo = ju
      jto = jt
      jeo = je
      jfo = jf
      jbo = jb
c
      ir = 1
      it = 2
      iu = 3
      jr = ir
      jt = it
      ju = iu
      jb = iu
c
      lady2   = lady(2)
      lady(2) = 0
      lady4   = lady(4)
      lady(4) = 0
c
      tauo  = tau
      tau   = 1.D-2
      neqo  = neqn
      neqn  = 3
c
      uextro = uextr
      uextr  = 1.D-10
c
      dtime = tau
c
c=======================================================================
c     NEWTON-RAPHSON iteration
c=======================================================================
c
      iter = 1
c
c.......................................................................
c     calculate derivatives 
c.......................................................................
c
        rn    (ngrs-2) = 1.0D0
        tg    (ngrs-2) = 85.D0
      dltdlr00(ngrs-2) = 1.0D0
        ug    (ngrs-2) = uextl
      dludlr00(ngrs-2) = 1.0D0
        tg    (ngrs-1) = 85.D0
      dltdlr00(ngrs-1) = 1.0D0
        ug    (ngrs-1) = uextl
      dludlr00(ngrs-1) = 1.0D0
      do 1 k = ngrs, ngre
c
        tg    (k) =     10.D0 + 75.D0 * ( rn(k   )-rn(ngre) )
     .                                / ( rn(ngrs)-rn(ngre) )
      dltdlr00(k) =             75.D0 / ( rn(ngrs)-rn(ngre) ) / tg(k)
c
      rmid   =      (  rn(k) - radc ) / rw
      rmid   = max  ( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(  0.5D0, 0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        ug    (k) = uextr + 0.5D0*( 1.0D0 - f0 (rmid) )*(uextl - uextr)
      dludlr00(k) = - switch*( uextl - uextr ) * 0.5D0/rw * f1(rmid)**2
c
    1 continue
        tg    (ngre+1) = 10.D0
      dltdlr00(ngre+1) = 1.0D0
        ug    (ngre+1) = uextr
      dludlr00(ngre+1) = 1.0D0
        tg    (ngre+2) = 10.D0
      dltdlr00(ngre+2) = 1.0D0
        ug    (ngre+2) = uextr
      dludlr00(ngre+2) = 1.0D0
c
      do 2 k = ngrs - 2, ngre + 2
      tn  (k) = tg(k)
      un  (k) = ug(k)
      bri (k) = abs(ug(k))
      ro  (k) = rn(k)
    2 continue
c
c.......................................................................
c     initialize all matrices to zero
c.......................................................................
c
  100 continue
      do 30 i = 1, meqn
      do 20 j = 1, meqn
      do 10 k = 1, mgr
      em2(i, j, k) = 0.0D0
      em1(i, j, k) = 0.0D0
      e00(i, j, k) = 0.0D0
      ep1(i, j, k) = 0.0D0
      ep2(i, j, k) = 0.0D0
   10 continue
   20 continue
   30 continue
c
      do 50 i = 1, meqn
      do 40 k = 1, mgr
      rhs(i, k) = 0.0D0
   40 continue
   50 continue
c
      do 70 i = 1, 3*mpd - 2
      do 60 j = 1, meqn*mgr
      bm(i, j) = 0.0D0
   60 continue
   70 continue
c
      do 80 i = 1, meqn*mgr
      br(i) = 0.0D0
   80 continue
c
c.......................................................................
c     calculate matrix elements in continuity and gas energy equation
c.......................................................................
c
      e00(it, jt, ngrs-2) =  1.0D0
      e00(id, jd, ngrs-2) =  1.0D0
c
      do 3 k = ngrs - 1, ngre + 1
c
      rhs(it,     k) =  tn(k) - tg      (k)
      e00(it, jt, k) =  tn(k) 
      e00(it, jr, k) =        - dltdlr00(k) * rn(k)
c
      rhs(iu,     k) =  un(k) - ug      (k)
      e00(iu, ju, k) =  un(k) 
      e00(iu, jr, k) =        - dludlr00(k) * rn(k)
    3 continue
c
      e00(it, jt, ngre+2) =  1.0D0
      e00(iu, ju, ngre+2) =  1.0D0
c
c.......................................................................
c     calculate matrix elements in grid equation
c.......................................................................
c
      call grid
c
c.......................................................................
c     solve the system
c.......................................................................
c
      call matslv
c
c.......................................................................
c     find maximum changes 
c.......................................................................
c
      kdmx = isamax(ngre - ngrs + 2, rhs(ir,ngrs), meqn)
      dmax = abs( rhs(ir,ngrs + kdmx - 1) ) + 1.0D-30
c
      do 4 k = ngrs, ngre
      cx  (k) = ( rhs(ir, k+1) * rn(k+1) - rhs(ir, k) * rn(k) )
     .        / (                rn(k+1) -              rn(k) )
    4 continue
      kcmx = isamax( ngre - ngrs + 1, cx(ngrs), 1 )
      cmax = abs( cx(ngrs + kcmx - 1) ) + 1.0D-30
c
      sc1 = min( 1.0D0, dtol / dmax )
      sc2 = min( 1.0D0, ctol / cmax )
      sc  = min( sc1, sc2 )
c
      do 5 k = ngrs, ngre + 1
      ro  (k) =   rn(k)
      rn  (k) = ( 1.0D0 + sc * rhs(ir,k) ) * ro(k)
      cx  (k) = ( rn(k) - ro(k) ) / ro(k)
    5 continue
      ksmx = isamax(ngre - ngrs + 2, cx(ngrs), 1)
      smax = abs( cx(ngrs + ksmx - 1) ) + 1.0D-30
c
      iter = iter + 1
c.......................................................................
c     summary output
c.......................................................................
c
      write(itty,998) iter,timen,ksmx,smax,kdmx,dmax
  998 format(" iter ="i3" time ="1pe10.3
     .       " d_r/ro("i3") ="e9.2" rhs("i3") ="e9.2 )
c
c-----------------------------------------------------------------------
c     check for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre + 1
      if ( rn(k+1) .le. rn(k) ) then
           write(itty,'(" grid extrapolation nonmonotonic")')
           write(iout,'(" grid extrapolation nonmonotonic")')
           call clzfil
           stop 'gridinit09'
           end if
    6 continue
c
c-----------------------------------------------------------------------
c     compute new timestep and update quantities
c-----------------------------------------------------------------------
c
      timeo = timen
      timen = timeo + dtime
c
        rn    (ngrs-2) = 1.0D0
        tg    (ngrs-2) = 85.D0
      dltdlr00(ngrs-2) = 1.0D0
        ug    (ngrs-2) = uextl
      dludlr00(ngrs-2) = 1.0D0
        tg    (ngrs-1) = 85.D0
      dltdlr00(ngrs-1) = 1.0D0
        ug    (ngrs-1) = uextl
      dludlr00(ngrs-1) = 1.0D0
      do 7 k = ngrs, ngre
c
        tg    (k) =     10.D0 + 75.D0 * ( rn(k   )-rn(ngre) )
     .                                / ( rn(ngrs)-rn(ngre) )
      dltdlr00(k) =             75.D0 / ( rn(ngrs)-rn(ngre) ) / tg(k)
c
      rmid   =      (  rn(k) - radc ) / rw
      rmid   = max  ( -3.D+2, min( rmid, +3.D+2 ) )
      switch = cvmgt(  0.5D0, 0.0D0, abs(f1 (rmid)) .lt. 1.D+150)
c
        ug    (k) = uextr + 0.5D0*( 1.0D0 - f0 (rmid) )*(uextl - uextr)
      dludlr00(k) = - switch*( uextl - uextr ) * 0.5D0/rw * f1(rmid)**2
c
    7 continue
        tg    (ngre+1) = 10.D0
      dltdlr00(ngre+1) = 1.0D0
        ug    (ngre+1) = uextr
      dludlr00(ngre+1) = 1.0D0
        tg    (ngre+2) = 10.D0
      dltdlr00(ngre+2) = 1.0D0
        ug    (ngre+2) = uextr
      dludlr00(ngre+2) = 1.0D0
c
      do 8 k = ngrs - 2, ngre + 2
      tn  (k) = tg(k)
      un  (k) = ug(k)
      bri (k) = abs(ug(k))
    8 continue
c
      if (smax .gt. 5.D-4) go to 100
c
c=======================================================================
c     restore original code parameters 
c=======================================================================
c
      ir = iro
      im = imo
      id = ido
      iu = iuo
      it = ito
      ie = ieo
      if = ifo
      jr = jro
      jm = jmo
      jd = jdo
      ju = juo
      jt = jto
      je = jeo
      jf = jfo
      jb = jbo
c
      neqn  = neqo
      tau   = tauo
      uextr = uextro
c
      lady(2) = lady2
      lady(4) = lady4
c
c=======================================================================
c
      return
      end
