      subroutine update
c
c***********************************************************************
c     update solution at advanced time level after NEWTON-RAPHSON cycle
c.......................................................................
c     calling sequence: step > update > update{r,h,rh}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     turbulent radiation hydrodynamics, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .ge.1) stop 'update'
c
c=======================================================================
c     radiation hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .eq.0) call updaterh
c
c=======================================================================
c     turbulent hydrodynamics without radiation, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .ge.1) stop 'update'
c
c=======================================================================
c     pure hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .eq.0) call updateh
c
c=======================================================================
c     time-dependent radiation (static medium)
c=======================================================================
c
      if (lhydr .eq. 0 .and. lrad .eq. 1 .and. ltur .eq.0) call updater
c
c=======================================================================
c
      if (lhydr .gt. 1 .or.  lrad .gt. 1                 ) stop 'update'
c
      return
      end
