      subroutine step
c
c***********************************************************************
c     perform one integration step
c.......................................................................
c     calling sequence: titan > step > peq, matgen, matslv, maxdel, 
c                                      update, clzfil, reset, 
c                                      timstp, totnrg
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     NEWTON-RAPHSON iteration cycle
c=======================================================================
c-----------------------------------------------------------------------
c     extrapolate old solution to new time level
c-----------------------------------------------------------------------
c
      call peq
c
c-----------------------------------------------------------------------
c     calculate matrix elements and solve system
c-----------------------------------------------------------------------
c
      jtry  = 1
      jback = 0
    2 iter  = 1
c
    3 call matgen
c
      call matslv
c
c-----------------------------------------------------------------------
c     find maximum changes
c-----------------------------------------------------------------------
c
      call maxdel
c
c-----------------------------------------------------------------------
c     update solution
c-----------------------------------------------------------------------
c
      call update
c
c-----------------------------------------------------------------------
c     check for convergence
c-----------------------------------------------------------------------
c
      if (dmax .le. conv ) go to 7
          iter = iter + 1
      if (iter .le. niter) go to 3
c
c=======================================================================
c     failure to converge; try again with reduced timestep
c=======================================================================
c
      write( itty,'(/34x " divergence. jtry = " i3)') jtry
      write( iout,'(/34x " divergence. jtry = " i3)') jtry
c
              tfac  = 5.0D-2
      dtime = tfac  * dtime
      timen = timeo + dtime
          jtry = jtry + 1
      if (jtry .gt. ntry) then
          write(iout,'(" failure to converge after"i3" attempts")') ntry
          call writout
          call clzfil
          stop 'step01'
      end if
c
      call reset
      go to 2
c
c=======================================================================
c     convergence: perform boost iteration
c=======================================================================
c
    7 lboos = 1
      call matgen
      call matslv
      call maxdel
      call update
      lboos = 0
c
c=======================================================================
c     choose a new timestep 
c=======================================================================
c
      call timstp
c
c=======================================================================
c     check total change in solution over timestep
c=======================================================================
c
          xs = smax / stol
c                 xchek = cvmgt( 2.0D0, 1.0D3, jstep .ne. jsteps)
 	          xchek =        1.0D6
      if (xs .lt. xchek) then
          call totnrg
          return
      end if
c
c-----------------------------------------------------------------------
c     total change is too large despite convergence; do integration over
c     with reduced timestep
c-----------------------------------------------------------------------
c
      write(itty, 101) xs
      write(iout, 101) xs
c
          jback = jback + 1
      if (jback .gt. nback) then
          write(itty, 102) nback
          write(iout, 102) nback
          call clzfil
          stop 'step02'
      end if
c
      write(itty,103) dtime
      write(iout,103) dtime
c
      call reset
      go to 2
c
c-----------------------------------------------------------------------
c
  101 format(/' changes over timestep'1pe9.2' times larger than limit')
  102 format( ' total changes still too large after'i2' attempts')
  103 format( ' repeat integration step with dtime ='1pe9.2)
c
      end
