      subroutine maxdel
c
c***********************************************************************
c     find maximum fractional change in all variables from 
c     NEWTON-RAPHSON iteration; scale them if necessary               
c.......................................................................
c     called by: step
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      ni = ngre - ngrs + 2
      nc = ngre - ngrs + 1
c
      do 1 j = 1, meqn
      kx  (j) = 0
      dx  (j) = 0.0D0
    1 continue
c
      do 2 j = ngrs, ngre
      cx  (j) = 0.0D0
    2 continue
c
c=======================================================================
c     turbulent radiation hydrodynamics, not yet available 
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .ge.1) stop 'maxdel'
c
c=======================================================================
c     radiation hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .eq. 0) then
c
	    kx(ir) = isamax(ni, rhs(ir, ngrs), meqn)
	    kx(im) = isamax(ni, rhs(im, ngrs), meqn)
	    kx(id) = isamax(nc, rhs(id, ngrs), meqn)
	    kx(iu) = isamax(ni, rhs(iu, ngrs), meqn)
	    kx(it) = isamax(nc, rhs(it, ngrs), meqn)
	    kx(ie) = isamax(nc, rhs(ie, ngrs), meqn)
	    kx(if) = isamax(ni, rhs(if, ngrs), meqn)
c
	    dx(ir) = rhs(ir, ngrs + kx(ir) - 1)
	    dx(im) = rhs(im, ngrs + kx(im) - 1)
	    dx(id) = rhs(id, ngrs + kx(id) - 1)
	    dx(iu) = rhs(iu, ngrs + kx(iu) - 1)
	    dx(it) = rhs(it, ngrs + kx(it) - 1)
	    dx(ie) = rhs(ie, ngrs + kx(ie) - 1)
	    dx(if) = rhs(if, ngrs + kx(if) - 1)
c
c           global maximum
c
            dmax = max( abs( dx(ir) ), abs( dx(im) ), abs( dx(id) ), 
     .                  abs( dx(it) ), abs( dx(ie) ) ) 
     .           + 1.0D-30
      end if
c
c=======================================================================
c     turbulent hydrodynamics without radiation, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .ge.1) stop 'maxdel'
c
c=======================================================================
c     pure hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .eq. 0) then
c
	    kx(ir) = isamax(ni, rhs(ir, ngrs), meqn)
	    kx(im) = isamax(ni, rhs(im, ngrs), meqn)
	    kx(id) = isamax(nc, rhs(id, ngrs), meqn)
	    kx(iu) = isamax(ni, rhs(iu, ngrs), meqn)
	    kx(it) = isamax(nc, rhs(it, ngrs), meqn)
c
	    dx(ir) = rhs(ir, ngrs + kx(ir) - 1)
	    dx(im) = rhs(im, ngrs + kx(im) - 1)
	    dx(id) = rhs(id, ngrs + kx(id) - 1)
	    dx(iu) = rhs(iu, ngrs + kx(iu) - 1)
	    dx(it) = rhs(it, ngrs + kx(it) - 1)
c
c           global maximum
c
	    dmax = max( abs( dx(ir) ), abs( dx(im) ), abs( dx(id) ), 
     .                  abs( dx(it) ) )
     .           + 1.0D-30
      end if
c
c=======================================================================
c     time-dependent radiation (static medium)
c=======================================================================
c
      if (lhydr .eq. 0 .and. lrad .eq. 1 .and. ltur. eq. 0) then
c
	    kx(ir) = isamax(ni, rhs(ir, ngrs), meqn)
	    kx(it) = isamax(nc, rhs(it, ngrs), meqn)
	    kx(ie) = isamax(nc, rhs(ie, ngrs), meqn)
	    kx(if) = isamax(ni, rhs(if, ngrs), meqn)
c
	    dx(ir) = rhs(ir, ngrs + kx(ir) - 1)
	    dx(it) = rhs(it, ngrs + kx(it) - 1)
	    dx(ie) = rhs(ie, ngrs + kx(ie) - 1)
	    dx(if) = rhs(if, ngrs + kx(if) - 1)
c
c           global maximum
c
	    dmax = max( abs( dx(it) ), abs( dx(ie) ) )
     .           + 1.0D-30
      end if
c
c=======================================================================
c
      if (lhydr .gt. 1 .or.  lrad .gt. 1                 ) stop 'maxdel'
c
c=======================================================================
c     maximum cell size change
c=======================================================================
c
      cmax = 1.0D-30
c
      if (lgrid .eq. 1) then
c
	    do 3 k = ngrs, ngre
	    cx(k) = ( rhs(ir, k+1) * rn(k+1) - rhs(ir, k) * rn(k) )
     .            / (                rn(k+1) -              rn(k) )
    3       continue
	    kcx = isamax( nc, cx(ngrs), 1 )
	    cmax = abs( cx(ngrs + kcx - 1) ) + 1.0D-30
c
      end if
c
c=======================================================================
c     summary output
c=======================================================================
c
            if (ltalk .gt. 1) write(itty, '(" ")')
            if (ltalk .gt. 0) write(iout, '(" ")')
c
c     normal iteration
c
      if (lboos .eq. 0) then
            if (ltalk .gt. 1) write(itty, 101) jstep, dtime, timen, iter
            if (ltalk .gt. 0) write(iout, 101) jstep, dtime, timen, iter
      end if
c
c     boost iteration
c
      if (lboos .eq. 1) then
            if (ltalk .gt. 1) write(itty, 102) jstep, dtime
            if (ltalk .gt. 0) write(itty, 101) jstep, dtime, timen, iter
            if (ltalk .gt. 0) write(iout, 102) jstep, dtime
            if (ltalk .gt. 0) write(iout, 101) jstep, dtime, timen, iter
      end if
c
      if (ltalk .gt. 1) then
	  write(itty,'(" kmax ="  7i10  )') (kx(i), i = 1, neqn)
          write(itty,'(" dmax ="1p7e10.2)') (dx(i), i = 1, neqn)
          write(itty,'(" cmax ="1p1e10.2)')  cmax
      end if
      if (ltalk .gt. 0) then
          write(iout,'(" kmax ="  7i10  )') (kx(i), i = 1, neqn)
          write(iout,'(" dmax ="1p7e10.2)') (dx(i), i = 1, neqn)
          write(iout,'(" cmax ="1p1e10.2)')  cmax
      end if
c
      if (lboos .eq. 1) return
      if (lboos .gt. 1) stop 'maxdel01'
c
c=======================================================================
c     if necessary scale changes before applying them
c=======================================================================
c
      sc1  = min( 1.0D0, dtol / dmax )
      sc2  = min( 1.0D0, ctol / cmax )
      sc   = min( sc1, sc2 )
c
      if (sc .lt. 1.0D0) then
     	    if (ltalk .gt. 1) 
     .          write(itty,'(" iteration scale factor ="1pe10.2)') sc
     	    if (ltalk .gt. 0) 
     .          write(iout,'(" iteration scale factor ="1pe10.2)') sc
c
            do 4 i = 1, neqn
            do 5 k = ngrs, ngre + 1
            rhs(i, k) = sc * rhs(i, k)
    5	    continue
    4	    continue
      end if
c
c=======================================================================
c
  101 format(' jstep ='i6' dtime ='1pe9.2' time ='e9.2' iter ='i3)
  102 format(' jstep ='i6' dtime ='1pe9.2' boost iteration ')
c
      return
      end
