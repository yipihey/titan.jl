      subroutine advecti
c
c***********************************************************************
c     advection of interface-centered quantities
c.......................................................................
c     called by: gasmom(h,rh}, radmom{r,rh}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq2'
c
      do 1 k = ngrs - 1, ngre + 1
c
c     equation c21
c
      dq(k) = qn(k+1) - qn(k)
    1 continue
c
      do 2 k = ngrs, ngre + 1
c
c     equation c22
c
      sdq = dq(k-1) + dq(k)
c
      denom(k) = sdq + cvmgm( - epsadv, + epsadv, sdq )
c
c     equation c23 & c26 - c28
c
      qr     (k) =     cadv *  dq(k-1) * dq(k)   / denom(k)
c
      dqrdqp1(k) = + ( cadv *  dq(k-1) - qr(k) ) / denom(k)
      dqrdq00(k) = -   cadv *( dq(k-1) - dq(k) ) / denom(k)
      dqrdqm1(k) = - ( cadv *  dq(k  ) - qr(k) ) / denom(k) 
c
c     equation c24 & c29 - c31
c
      qsn    (k) = cvmgp(  qr    (k), 0.0D0, dq(k-1) * dq(k) )
c
      dqsdqp1(k) = cvmgp( dqrdqp1(k), 0.0D0, dq(k-1) * dq(k) )
      dqsdq00(k) = cvmgp( dqrdq00(k), 0.0D0, dq(k-1) * dq(k) )
      dqsdqm1(k) = cvmgp( dqrdqm1(k), 0.0D0, dq(k-1) * dq(k) )
    2 continue
c
      do 3 k = ngrs, ngre
c
c     equation c20
c
      s(k) = cvmgp( 0.5D0, - 0.5D0, flow(k) )
c
c     equation c25
c
      qb(k) = 
     .   (0.5D0 + s(k)) * (        thet  * (qn(k  ) + 0.5D0 * qsn(k  ))
     .                  + (1.0D0 - thet) * (qo(k  ) + 0.5D0 * qso(k  )))
     . + (0.5D0 - s(k)) * (        thet  * (qn(k+1) - 0.5D0 * qsn(k+1))
     .                  + (1.0D0 - thet) * (qo(k+1) - 0.5D0 * qso(k+1)))
c
c     equations c33 - c36
c
      dqbdqp2(k) =
     .           (0.5D0 - s(k)) * thet * (      - 0.5D0 * dqsdqp1(k+1))
      dqbdqp1(k) =
     .           (0.5D0 + s(k)) * thet * (      + 0.5D0 * dqsdqp1(k  ))
     .         + (0.5D0 - s(k)) * thet * (1.0D0 - 0.5D0 * dqsdq00(k+1))
      dqbdq00(k) =
     .           (0.5D0 + s(k)) * thet * (1.0D0 + 0.5D0 * dqsdq00(k  ))
     .         + (0.5D0 - s(k)) * thet * (      - 0.5D0 * dqsdqm1(k+1))
      dqbdqm1(k) = 
     .           (0.5D0 + s(k)) * thet * (      + 0.5D0 * dqsdqm1(k  ))
    3 continue
c
c-----------------------------------------------------------------------
c
      return
      end
