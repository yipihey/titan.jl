      subroutine peq
c
c***********************************************************************
c     prepare for NEWTON-RAPHSON iteration
c.......................................................................
c     calling sequence: step > peq > peq{r,h,rh}
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
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .ge. 1)  stop 'peq'
c
c=======================================================================
c     radiation hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 1 .and. ltur .eq. 0)  call peqrh
c
c=======================================================================
c     turbulent hydrodynamics without radiation, not yet available
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .ge. 1)  stop 'peq'
c
c=======================================================================
c     pure hydrodynamics
c=======================================================================
c
      if (lhydr .eq. 1 .and. lrad .eq. 0 .and. ltur .eq. 0)  call peqh
c
c=======================================================================
c     time-dependent radiation (static medium)
c=======================================================================
c
      if (lhydr .eq. 0 .and. lrad .eq. 1 .and. ltur .eq. 0)  call peqr
c
c=======================================================================
c
      if (lhydr .gt. 1 .or.  lrad .gt. 1                  )  stop 'peq'
c
      return
      end
