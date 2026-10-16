      subroutine matgen
c
c***********************************************************************
c     generate all matrix elements
c.......................................................................
c     calling sequence: step > matgen > mass, contin, viscous, grid,
c                                  gasmomrh, gasmomh,
c                                  gasnrgrh, gasnrgh, gasnrgr,
c                                  radmomrh,          radmomr,
c                                  radnrgrh,          radnrgr
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     initialize all matrices to zero
c=======================================================================
c
      do 3 i = 1, meqn
      do 2 j = 1, meqn
      do 1 k = 1, mgr
      em2(i, j, k) = 0.0D0
      em1(i, j, k) = 0.0D0
      e00(i, j, k) = 0.0D0
      ep1(i, j, k) = 0.0D0
      ep2(i, j, k) = 0.0D0
    1 continue
    2 continue
    3 continue
c
      do 5 i = 1, meqn
      do 4 k = 1, mgr
      rhs(i, k) = 0.0D0
    4 continue
    5 continue
c
      do 7 i = 1, 3*mpd - 2
      do 6 j = 1, meqn*mgr
      bm(i, j) = 0.0D0
    6 continue
    7 continue
c
      do 8 i = 1, meqn*mgr
      br(i) = 0.0D0
    8 continue
c
c=======================================================================
c     turbulent-convective radiation hydrodynamics, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .ge. 1) stop'matgen'        
c
c=======================================================================
c     radiation hydrodynamics 
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .eq. 0) then
	    call mass
	    call contin
	    call viscous
	    call gasmomrh
	    call gasnrgrh
	    call radnrgrh
	    call radmomrh
            call grid
	    return
      end if
c
c=======================================================================
c     turbulent hydrodynamics without radiation, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .ge. 1) stop'matgen'
c
c=======================================================================
c     pure hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .eq. 0) then
	    call mass 
	    call contin
	    call viscous
	    call gasmomh
	    call gasnrgh
            call grid
	    return
      end if
c
c=======================================================================
c     time-dependent radiation (static medium)
c=======================================================================
c
      if (lhydr .eq. 0 .and. lrad .eq. 1 .and. ltur .eq. 0) then
            call gasnrgr
            call radnrgr
            call radmomr
            call grid
	    return
      end if
c
c=======================================================================
c
      if (lhydr .gt. 1 .or.  lrad .gt. 1                  ) stop'matgen'
c
      return
      end
