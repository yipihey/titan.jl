      subroutine clzfil
c
c***********************************************************************
c     close files
c.......................................................................
c     called by: titan, step, peq{h,rh,th,trh}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      if (iin   .gt. 0) close (unit = iin  )
      if (iout  .gt. 0) close (unit = iout )
      if (idoc  .gt. 0) close (unit = idoc )
      if (ieos  .gt. 0) close (unit = ieos )
      if (iopac .gt. 0) close (unit = iopac)
      if (idump .gt. 0) close (unit = idump)
      if (ihist .gt. 0) close (unit = ihist)
      if (inrg  .gt. 0) close (unit = inrg )
c
      return
      end
