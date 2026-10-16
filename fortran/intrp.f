      subroutine intrp(nk, x, y, f, fx, fy, fxy, xx, yy, e, ex, ey, fl)
c
c***********************************************************************
c     interpolate eos tables using monotonized bicubic hermite         
c     polynomials                                                   
c.......................................................................
c     called by: eos, opac
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
c
c.......................................................................
c     input:
      dimension x(nk), y(nk)
      dimension f(mxe, mye), fx(mxe, mye), fy(mxe, mye), fxy(mxe, mye)
c.......................................................................
c     output:
      dimension xx(nk), yy(nk), e(nk), ex(nk), ey(nk), fl(nk)
c.......................................................................
c     internal:
      dimension xi(mgr), yi(mgr), ii(mgr), jj(mgr) 
c
      dimension f00(mgr), fx00(mgr), fy00(mgr), fxy00(mgr),
     .          f01(mgr), fx01(mgr), fy01(mgr), fxy01(mgr),
     .          f10(mgr), fx10(mgr), fy10(mgr), fxy10(mgr),
     .          f11(mgr), fx11(mgr), fy11(mgr), fxy11(mgr)
c
      dimension h1 (mgr), h2 (mgr), h3 (mgr), h4 (mgr),
     .          hx1(mgr), hx2(mgr), hx3(mgr), hx4(mgr)
c
      dimension g1 (mgr), g2 (mgr), g3 (mgr), g4 (mgr),
     .          gy1(mgr), gy2(mgr), gy3(mgr), gy4(mgr)
c
c-----------------------------------------------------------------------
c     assemble data needed for interpolation
c-----------------------------------------------------------------------
c
      do 1 k = 1, nk
c
c     forbid points off edges of table
c
      xx(k) = max( xmine, min( x(k), xmaxe ) )
      yy(k) = max( ymine, min( y(k), ymaxe ) )
c
c     flag points off edges of table
c
      fl(k) = cvmgt( 1.0D0, 0.0D0, x(k) .ne. xx(k) .or. y(k) .ne. yy(k))
c
c     compute indices of lower left corner of the cell
c
      ii(k) = min( int((xx(k) - xmine) / dxe + 1.00001D0), mxe - 1 )
      jj(k) = min( int((yy(k) - ymine) / dye + 1.00001D0), mye - 1 )
c
c     compute coordinates (xi, yi) relative to lower left corner of
c     cell in units of cell dimensions 
c
      xi(k) = (xx(k) - xmine) / dxe - float( ii(k) - 1 ) 
      yi(k) = (yy(k) - ymine) / dye - float( jj(k) - 1 ) 
c
c     assemble f, fx, fy, fxy at corners of each cell
c
      f  00(k) = f  (ii(k)    , jj(k)    )
      fx 00(k) = fx (ii(k)    , jj(k)    )
      fy 00(k) = fy (ii(k)    , jj(k)    )
      fxy00(k) = fxy(ii(k)    , jj(k)    )
c
      f  01(k) = f  (ii(k)    , jj(k) + 1)
      fx 01(k) = fx (ii(k)    , jj(k) + 1)
      fy 01(k) = fy (ii(k)    , jj(k) + 1)
      fxy01(k) = fxy(ii(k)    , jj(k) + 1)
c
      f  10(k) = f  (ii(k) + 1, jj(k)    )
      fx 10(k) = fx (ii(k) + 1, jj(k)    )
      fy 10(k) = fy (ii(k) + 1, jj(k)    )
      fxy10(k) = fxy(ii(k) + 1, jj(k)    )
c
      f  11(k) = f  (ii(k) + 1, jj(k) + 1)
      fx 11(k) = fx (ii(k) + 1, jj(k) + 1)
      fy 11(k) = fy (ii(k) + 1, jj(k) + 1)
      fxy11(k) = fxy(ii(k) + 1, jj(k) + 1)
c
    1 continue
c
c-----------------------------------------------------------------------
c     evaluate interpolant 
c-----------------------------------------------------------------------
c
      do 2 k = 1, nk
c
c     evaluate basis functions at (xi, yi)
c
      h 2(k) = - (2.0D0 * xi(k) - 3.0D0 ) * xi(k)**2
      g 2(k) = - (2.0D0 * yi(k) - 3.0D0 ) * yi(k)**2
      h 1(k) =    1.0D0 - h2(k)
      g 1(k) =    1.0D0 - g2(k)
      h 3(k) = ( (xi(k) - 2.0D0) * xi(k) + 1.0D0 ) * xi(k) * dxe
      g 3(k) = ( (yi(k) - 2.0D0) * yi(k) + 1.0D0 ) * yi(k) * dye
      h 4(k) =   (xi(k) - 1.0D0) * xi(k)**2 * dxe
      g 4(k) =   (yi(k) - 1.0D0) * yi(k)**2 * dye
      hx2(k) = - 6.0D0 * (xi(k) - 1.0D0) * xi(k) / dxe
      gy2(k) = - 6.0D0 * (yi(k) - 1.0D0) * yi(k) / dye
      hx1(k) = - hx2(k)
      gy1(k) = - gy2(k)
      hx3(k) = (3.0D0 * xi(k) - 4.0D0) * xi(k) + 1.0D0
      gy3(k) = (3.0D0 * yi(k) - 4.0D0) * yi(k) + 1.0D0
      hx4(k) = (3.0D0 * xi(k) - 2.0D0) * xi(k)
      gy4(k) = (3.0D0 * yi(k) - 2.0D0) * yi(k)
c
c     compute e, ex, and ey at (xx, yy)
c
      e (k) = g 1(k) * ( h 1(k) * f  00(k) + h 2(k) * f  10(k)
     .                 + h 3(k) * fx 00(k) + h 4(k) * fx 10(k) )
     .      + g 2(k) * ( h 1(k) * f  01(k) + h 2(k) * f  11(k)
     .                 + h 3(k) * fx 01(k) + h 4(k) * fx 11(k) )
     .      + g 3(k) * ( h 1(k) * fy 00(k) + h 2(k) * fy 10(k)
     .                 + h 3(k) * fxy00(k) + h 4(k) * fxy10(k) )
     .      + g 4(k) * ( h 1(k) * fy 01(k) + h 2(k) * fy 11(k)
     .                 + h 3(k) * fxy01(k) + h 4(k) * fxy11(k) )
c
      ex(k) = g 1(k) * ( hx1(k) * f  00(k) + hx2(k) * f  10(k)
     .                 + hx3(k) * fx 00(k) + hx4(k) * fx 10(k) )
     .      + g 2(k) * ( hx1(k) * f  01(k) + hx2(k) * f  11(k)
     .                 + hx3(k) * fx 01(k) + hx4(k) * fx 11(k) )
     .      + g 3(k) * ( hx1(k) * fy 00(k) + hx2(k) * fy 10(k)
     .                 + hx3(k) * fxy00(k) + hx4(k) * fxy10(k) )
     .      + g 4(k) * ( hx1(k) * fy 01(k) + hx2(k) * fy 11(k)
     .                 + hx3(k) * fxy01(k) + hx4(k) * fxy11(k) )
c
      ey(k) = gy1(k) * ( h 1(k) * f  00(k) + h 2(k) * f  10(k)
     .                 + h 3(k) * fx 00(k) + h 4(k) * fx 10(k) )
     .      + gy2(k) * ( h 1(k) * f  01(k) + h 2(k) * f  11(k)
     .                 + h 3(k) * fx 01(k) + h 4(k) * fx 11(k) )
     .      + gy3(k) * ( h 1(k) * fy 00(k) + h 2(k) * fy 10(k)
     .                 + h 3(k) * fxy00(k) + h 4(k) * fxy10(k) )
     .      + gy4(k) * ( h 1(k) * fy 01(k) + h 2(k) * fy 11(k)
     .                 + h 3(k) * fxy01(k) + h 4(k) * fxy11(k) )
c
    2 continue
c
      return
      end
