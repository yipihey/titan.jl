c-----------------------------------------------------------------------
c     PORT phase 1: callable entry points replacing program titan.
c     The bodies mirror titan.f exactly (same call sequence); the DO loop
c     over jstep becomes one titan_step call per iteration.
c     istat: 0 = continue, 1 = finished (timen >= tmax), 2 = finished
c     (jstep > jstepe).  A Fortran STOP inside the solver still
c     terminates the whole process (as in the original).
c-----------------------------------------------------------------------
      block data titanbd
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.dat'
      end
c
      subroutine titan_init(istat)
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      common /portctl/ ldone
      xindef = 0.0D0
      ldone  = 0
      call opnfil
      call start
      jstep = jsteps
      istat = 0
      if (jstep .gt. jstepe) istat = 2
      return
      end
c
      subroutine titan_step(istat)
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      common /portctl/ ldone
      istat = 0
      if (ldone .ne. 0) then
          istat = ldone
          return
      end if
      call step
      call archive
      call writout
      call writsty
      if (timen .ge. tmax) then
          call clzfil
          ldone = 1
          istat = 1
          return
      end if
      jstep = jstep + 1
      if (jstep .gt. jstepe) then
          call clzfil
          ldone = 2
          istat = 2
      end if
      return
      end
c
      subroutine titan_finalize()
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      common /portctl/ ldone
      if (ldone .eq. 0) call clzfil
      ldone = 3
      if (itty .gt. 0) close(unit = itty)
      return
      end
