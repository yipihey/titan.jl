      subroutine intrc (ks,ke,x,yin,xin,yout,yout1)
c
c***********************************************************************
c     interpolation with monotonized piecewise cubic polynomials
c
c     source: M. STEFFEN (1990) AA, 239, 443
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
c
c.......................................................................
c     input:
      dimension x(mgr), yin(mgr), xin(mgr)
c.......................................................................
c     output:
      dimension yout(mgr), yout1(mgr)
c.......................................................................
c     internal:
      dimension h(mgr), s(mgr), p(mgr), y1(mgr), 
     .          a(mgr), b(mgr), c(mgr),  d(mgr) 
c
c=======================================================================
c-----------------------------------------------------------------------
c     mesh spacing and slope of the secant through the data
c-----------------------------------------------------------------------
c
      do 1 k = ks, ke - 1
      h(k) =     x(k+1) -   x(k)
      s(k) = ( yin(k+1) - yin(k) ) / h(k)
    1 continue
c
c-----------------------------------------------------------------------
c     parabola through the data
c-----------------------------------------------------------------------
c
           k  = ks
      p   (k) =   s(k  ) * ( 1.D0 + h(k) / ( h(k) + h(k+1) ) ) 
     .          - s(k+1) *          h(k) / ( h(k) + h(k+1) ) 
c
      do 2 k  = ks + 1, ke - 1
      p   (k) = ( s(k-1) * h(k) + s(k) * h(k-1) ) / ( h(k-1) + h(k) )
    2 continue
c
           k  = ke
      p   (k) =   s(k-1) * ( 1.D0 + h(k-1) / ( h(k-1) + h(k-2) ) ) 
     .          - s(k-2) *          h(k-1) / ( h(k-1) + h(k-2) ) 
c
c-----------------------------------------------------------------------
c     slope to satisfy monotonicity
c-----------------------------------------------------------------------
c
           k  = ks
      y1  (k) = ( sign(1.D0, s(k  ) ) + sign(1.D0, p(k) ) )
     .          *  min( abs( s(k  ) ), 0.5D0* abs( p(k) ) )
c
      do 3 k  = ks + 1, ke - 1
      y1  (k) = ( sign(1.D0, s(k-1) ) +              sign(1.D0, s(k) ) )
     .          *  min( abs( s(k-1) ), 0.5D0* abs( p(k) ), abs( s(k) ) )
    3 continue
c
           k  = ke
      y1  (k) = ( sign(1.D0, s(k-1) ) + sign(1.D0, p(k) ) )
     .          *  min( abs( s(k-1) ), 0.5D0* abs( p(k) ) )
c
c-----------------------------------------------------------------------
c     coefficients of cubic polynomial
c-----------------------------------------------------------------------
c
      do 4 k  = ks, ke - 1
      a   (k) = + (        y1(k) + y1(k+1) - 2.D0 * s(k) ) / h(k)**2
      b   (k) = - ( 2.D0 * y1(k) + y1(k+1) - 3.D0 * s(k) ) / h(k)
      c   (k) =            y1(k)
      d   (k) =           yin(k)
    4 continue
c
c-----------------------------------------------------------------------
c    find interval [x(k),x(k+1)] containing xin(j) and interpolate
c-----------------------------------------------------------------------
c
      yout (ks) = yin(ks)
      yout1(ks) = 0.D0
      do 5 k  = ks + 1, ke - 1
           do 6 j  = ke - 1, ks, - 1
	   if ( ( xin(j)-x(k) ) .ge. 0.D0 ) go to 7
    6      continue
    7      continue
      yout (k) = ( ( a(k)  *( xin(j)-x(k) ) 
     .             + b(k) )*( xin(j)-x(k) )
     .           +   c(k) )*( xin(j)-x(k) ) + d(k)
      yout1(k) =   ( a(k)  *( xin(j)-x(k) ) 
     .             + b(k) )*( xin(j)-x(k) ) + c(k) 
    5 continue
      yout (ke) = yin(ke)
      yout1(ke) = 0.D0
c
c-----------------------------------------------------------------------
c
      return
      end
