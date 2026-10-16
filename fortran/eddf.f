      subroutine eddf(nang, ang, xi, sum0, sum2)
c
c***********************************************************************
c     calculate zeroth and second moment of specific intensity
c.......................................................................
c     called by: eddfac
c***********************************************************************
c
      include'titan.imp'
c
      dimension ang(1), xi(1)
c
      do 1 i = 1, nang - 1
c
c     equation ef48
c
      sum0 = sum0 + 0.5D0 * (xi(i) + xi(i+1)) * (ang(i+1) - ang(i))
c
c     equation ef50
c
      sum2 = sum2 
     .            + (xi(i) * ang(i+1) - xi(i+1) * ang(i))
     .            * (ang(i)**2 + ang(i) * ang(i+1) + ang(i+1)**2) / 3.D0
     .
     .            + (xi(i+1) - xi(i)) * (ang(i)    + ang(i+1))
     .                                * (ang(i)**2 + ang(i+1)**2) / 4.D0
    1 continue
c
      return
      end
