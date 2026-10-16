      subroutine mass
c
c***********************************************************************
c     cumulative mass equation m2
c.......................................................................
c     called by: matgen
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     interior zones
c=======================================================================
c
      do 1 k = ngrs + 1, ngre + 1
c
c     equation m14
c
      rhs(im, k) = xmen(k-1) - xmen(k) - dn  (k-1) * dvoln (k-1)
c
c     equations m9 - m13
c
      e00(im, jm, k) = e00(im, jm, k)  - xmen(k  )
      em1(im, jm, k) = em1(im, jm, k)  + xmen(k-1)
c
      e00(im, jr, k) = e00(im, jr, k)  - dn  (k-1) * rmup1n(k  )
      em1(im, jr, k) = em1(im, jr, k)  + dn  (k-1) * rmup1n(k-1)
c
      em1(im, jd, k) = em1(im, jd, k)  - dn  (k-1) * dvoln (k-1)
c
    1 continue
c
c=======================================================================
c     inner boundary condition
c=======================================================================
c
c     eulerian boundary bc14, or lagrangean boundary bc65 (phil0 == 0)
c
      k = ngrs
c
      rhs(im,     k) = xmen(k) - xmeo(k)      + rmu  (k) * phil0 * dtime
      e00(im, jr, k) =    + thet * rn(k) *xmu * rmum1(k) * phil0 * dtime
      e00(im, jm, k) = xmen(k)
c
c=======================================================================
c     phantom zones
c=======================================================================
c
      do 2 k = 1, ngrs - 2
      rhs(im,     k) = 0.0D0
      em2(im, jm, k) = 0.0D0
      em1(im, jm, k) = 0.0D0
      e00(im, jm, k) = 1.0D0
      ep1(im, jm, k) = 0.0D0
      ep2(im, jm, k) = 0.0D0
    2 continue
c
      do 3 k = ngre + 3, mgr
      rhs(im,     k) = 0.0D0
      em2(im, jm, k) = 0.0D0
      em1(im, jm, k) = 0.0D0
      e00(im, jm, k) = 1.0D0
      ep1(im, jm, k) = 0.0D0
      ep2(im, jm, k) = 0.0D0
    3 continue
c
c-----------------------------------------------------------------------
c     eulerian inner boundary
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
c     zero flux; equations pz3 & pz4
c
      if (leibc .eq. 1) then
c
           rhs(im,     k) = xmen(k) - 2.0D0 * xmen(k+1) + xmen(k+2)
c
           ep2(im, jm, k) =                             + xmen(k+2)
           ep1(im, jm, k) =         - 2.0D0 * xmen(k+1)
           e00(im, jm, k) = xmen(k)
c
      end if
c
c     nonzero flux; equations pz13 & pz14
c
      if (leibc .eq. 2) then
c
           rhs(im,     k) = - xmen(k) + xmen(k+1) + dl * dvoln (k  )
c
           ep1(im, jm, k) =           + xmen(k+1)
           e00(im, jm, k) = - xmen(k)
c
           ep1(im, jr, k) =                       + dl * rmup1n(k+1)
           e00(im, jr, k) =                       - dl * rmup1n(k  )
c
      end if
c
c-----------------------------------------------------------------------
c     lagrangean inner boundary
c-----------------------------------------------------------------------
c
c     equations pz23 & pz24
c
      if (llibc .gt. 0) then
c
           rhs(im,     k) = - xmen(k) + xmen(k+1) + delml
c
           ep1(im, jm, k) =           + xmen(k+1)
           e00(im, jm, k) = - xmen(k)
c
      end if
c
c-----------------------------------------------------------------------
c     eulerian outer boundary
c-----------------------------------------------------------------------
c
      k = ngre + 2
c
c     zero flux; equations pz41 & pz42
c
      if (leobc .eq. 1) then
c
           rhs(im,     k) = xmen(k-2) - 2.0D0 * xmen(k-1) + xmen(k)
c
           e00(im, jm, k) =                               + xmen(k)
           em1(im, jm, k) =           - 2.0D0 * xmen(k-1)
           em2(im, jm, k) = xmen(k-2)
c
      end if
c
c     nonzero flux; equations pz51 & pz52
c
      if (leobc .eq. 2) then
c
           rhs(im,     k) = xmen(k-1) - xmen(k) - dr      * dvoln (k-1)
c
           e00(im, jm, k) =           - xmen(k)
           em1(im, jm, k) = xmen(k-1)
c
           e00(im, jr, k) =                     - dr      * rmup1n(k  )
           em1(im, jr, k) =                     + dr      * rmup1n(k-1)
c
      end if
c
c     transmitting; equations pz61 & pz62
c
      if (leobc .eq. 3) then
c
           rhs(im,     k) = xmen(k-1) - xmen(k) - dn(k-1) * dvoln (k-1)
c
           e00(im, jm, k) =           - xmen(k)
           em1(im, jm, k) = xmen(k-1)
c
           e00(im, jr, k) =                     - dn(k-1) * rmup1n(k  )
           em1(im, jr, k) =                     + dn(k-1) * rmup1n(k-1)
c
           em1(im, jd, k) =                     - dn(k-1) * dvoln (k-1)
c
      end if
c
c-----------------------------------------------------------------------
c     lagrangean outer boundary 
c-----------------------------------------------------------------------
c
c     equations pz71 & pz72
c
      if (llobc .gt. 0) then
           rhs(im,     k) = xmen(k-1) - xmen(k) + delmr
c
           e00(im, jm, k) =           - xmen(k)
           em1(im, jm, k) = xmen(k-1)
      end if
c
c-----------------------------------------------------------------------
c
      return
      end
