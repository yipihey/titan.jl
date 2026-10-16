      subroutine eddg(nang, ang, xi, sum1)
c
c***********************************************************************
c     calculate flux EDDINGTON factor at boundaries ef48 - ef51
c.......................................................................
c     called by: eddfac
c***********************************************************************
c
      include'titan.imp'
c
      dimension ang(1), xi(1)
c
      sum0 = 0.0D0
      sum1 = 0.0D0
c
      do 1 i = 1, nang - 1
c
      sum0 = sum0 + 0.5D0 * (xi(i) + xi(i+1)) * (ang(i+1) - ang(i))
c
      sum1 = sum1
     .           + (xi(i) * ang(i+1) - xi(i+1) * ang(i))
     .                                     * (ang(i) + ang(i+1)) / 2.0D0
     .
     .           + (xi(i+1) - xi(i))
     .           * (ang(i)**2 + ang(i) * ang(i+1) + ang(i+1)**2) / 3.0D0
    1 continue
c
      sum1 = sum1 / sum0
c
      return
      end
