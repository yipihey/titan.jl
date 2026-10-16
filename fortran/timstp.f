      subroutine timstp
c
c***********************************************************************
c     choose a new timestep
c.......................................................................
c     called by: step
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     calculate maximum fractional change over the last timestep
c=======================================================================
c
      ni = ngre - ngrs + 2
      nc = ngre - ngrs + 1
c
      do 10 j = 1, neqn
      kx   (j) = 0
      dx   (j) = 0.0D0
   10 continue
c
c=======================================================================
c     turbulent radiation hydrodynamics, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .ge. 1) stop'timstp'
c
c=======================================================================
c     radiation hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .eq. 0) then
c
	    do 2 k = ngrs, ngre + 1
	    sx(k, ir) = ( rn  (k) - ro  (k) ) / ro   (k)
	    sx(k, im) = ( xmen(k) - xmeo(k) ) / xmeo (k)
	    sx(k, id) = ( dn  (k) - do  (k) ) / do   (k)
            sx(k, iu) = ( un  (k) - uo  (k) ) / unom (k)
	    sx(k, it) = ( tn  (k) - to  (k) ) / to   (k)
	    sx(k, ie) = ( ern (k) - ero (k) ) / ero  (k)
            sx(k, if) = ( frn (k) - fro (k) ) / frnom(k)
    2       continue
c
	    kx(ir) = isamax(ni, sx(ngrs, ir), 1)
	    kx(im) = isamax(ni, sx(ngrs, im), 1)
	    kx(id) = isamax(nc, sx(ngrs, id), 1)
            kx(iu) = isamax(ni, sx(ngrs, iu), 1)
	    kx(it) = isamax(nc, sx(ngrs, it), 1)
	    kx(ie) = isamax(nc, sx(ngrs, ie), 1)
            kx(if) = isamax(ni, sx(ngrs, if), 1)
c
	    dx(ir) = sx(ngrs + kx(ir) - 1, ir)
	    dx(im) = sx(ngrs + kx(im) - 1, im)
	    dx(id) = sx(ngrs + kx(id) - 1, id)
            dx(iu) = sx(ngrs + kx(iu) - 1, iu)
	    dx(it) = sx(ngrs + kx(it) - 1, it)
	    dx(ie) = sx(ngrs + kx(ie) - 1, ie)
            dx(if) = sx(ngrs + kx(if) - 1, if)
c
c           global maximum
c
            smax = max( abs( dx(ir) ), abs( dx(im) ), abs( dx(id) ), 
     .                  abs( dx(it) ), abs( dx(ie) ) ) + 1.0D-30
c
      end if
c
c=======================================================================
c     turbulent hydrodynamics without radiation, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .ge. 1) stop'timstp'
c
c=======================================================================
c     pure hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .eq. 0) then
c
	    do 4 k = ngrs, ngre + 1
	    sx(k, ir) = ( rn  (k) - ro  (k) ) / ro  (k)
            sx(k, im) = ( xmen(k) - xmeo(k) ) / xmeo(k)
	    sx(k, id) = ( dn  (k) - do  (k) ) / do  (k)
	    sx(k, iu) = ( un  (k) - uo  (k) ) / unom(k)
	    sx(k, it) = ( tn  (k) - to  (k) ) / to  (k)
    4       continue
c
	    kx(ir) = isamax(ni, sx(ngrs, ir), 1)
	    kx(im) = isamax(ni, sx(ngrs, im), 1)
	    kx(id) = isamax(nc, sx(ngrs, id), 1)
	    kx(iu) = isamax(ni, sx(ngrs, iu), 1)
	    kx(it) = isamax(nc, sx(ngrs, it), 1)
c
	    dx(ir) = sx(ngrs + kx(ir) - 1, ir)
	    dx(im) = sx(ngrs + kx(im) - 1, im)
	    dx(id) = sx(ngrs + kx(id) - 1, id)
	    dx(iu) = sx(ngrs + kx(iu) - 1, iu)
	    dx(it) = sx(ngrs + kx(it) - 1, it)
c
c           global maximum
c
            smax = max( abs( dx(ir) ), abs( dx(im) ), abs( dx(id) ), 
     .                  abs( dx(it) ) ) + 1.0D-30
c
      end if
c
c=======================================================================
c     time-dependent radiation (static medium)
c=======================================================================
c
      if (lhydr .eq. 0 .and. lrad .eq. 1 .and. ltur .eq. 0) then
c
	    do 5 k = ngrs, ngre + 1
	    sx(k, it) = ( tn (k) - to (k) ) / to   (k)
	    sx(k, ie) = ( ern(k) - ero(k) ) / ero  (k)
	    sx(k, if) = ( frn(k) - fro(k) ) / frnom(k)
    5       continue
c
	    kx(it) = isamax(nc, sx(ngrs, it), 1)
	    kx(ie) = isamax(nc, sx(ngrs, ie), 1)
	    kx(if) = isamax(ni, sx(ngrs, if), 1)
c
	    dx(it) = sx(ngrs + kx(it) - 1, it)
	    dx(ie) = sx(ngrs + kx(ie) - 1, ie)
	    dx(if) = sx(ngrs + kx(if) - 1, if)
c
c           global maximum
c
	    smax = max( abs( dx(it) ), abs( dx(ie) ) ) + 1.0D-30
c
      end if
c
c=======================================================================
c     choose new timestep
c=======================================================================
c
      x = stol / smax
      phi1 = cvmgt( sqrt(x), x**2, x .gt. 1.0D0 )
c
      phi2 = max( 0.3D0, min( 3.0D0, phi1 ) ) 
      tfac = max( 0.1D0, min( 10.D0, sqrt( tfac * phi2 ) ) )
c
      dtime = tfac * dtime
c
      return
      end
