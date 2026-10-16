      subroutine advectc
c
c***********************************************************************
c     advection of cell-centered quantities
c.......................................................................
c     called by: contin, gasnrg{h,rh}, radnrg{r,rh}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq1'
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
      sdq   = dq(k-1) + dq(k)
c
      denom(k) = sdq + cvmgp( + epsadv, - epsadv, sdq )
c
c     equation c23 & c26 - c28
c
      qr      (k) =             cadv *  dq(k-1) * dq(k)  / denom(k)
c
      dqrdlqp1(k) = + qn(k+1) * cadv * (dq(k-1)          / denom(k))**2
      dqrdlq00(k) = - qn(k  ) * cadv * (dq(k-1) - dq(k)) / denom(k)
      dqrdlqm1(k) = - qn(k-1) * cadv * (          dq(k)  / denom(k))**2
c
c     equation c24 & c29 - c31
c
      qsn     (k) = cvmgp(  qr     (k), 0.0D0, dq(k-1) * dq(k) )
c
      dqsdlqp1(k) = cvmgp( dqrdlqp1(k), 0.0D0, dq(k-1) * dq(k) )
      dqsdlq00(k) = cvmgp( dqrdlq00(k), 0.0D0, dq(k-1) * dq(k) )
      dqsdlqm1(k) = cvmgp( dqrdlqm1(k), 0.0D0, dq(k-1) * dq(k) )
    2 continue
c
      do 3 k = ngrs + 1, ngre + 1
c
c     equation c20
c
      s(k) = cvmgp( 0.5D0, - 0.5D0, flow(k) )
c
c     equation c25
c
      qb(k) =
     .   (0.5D0 + s(k)) * (        thet  * (qn(k-1) + 0.5D0 * qsn(k-1))
     .                  + (1.0D0 - thet) * (qo(k-1) + 0.5D0 * qso(k-1)))
     . + (0.5D0 - s(k)) * (        thet  * (qn(k  ) - 0.5D0 * qsn(k  ))
     .                  + (1.0D0 - thet) * (qo(k  ) - 0.5D0 * qso(k  )))
c
c     equations c33 - c36
c
      dqbdlqp1(k) = 
     .        (0.5D0 - s(k)) * thet * (        - 0.5D0 * dqsdlqp1(k  ))
      dqbdlq00(k) =
     .        (0.5D0 + s(k)) * thet * (        + 0.5D0 * dqsdlqp1(k-1))
     .      + (0.5D0 - s(k)) * thet * (qn(k  ) - 0.5D0 * dqsdlq00(k  )) 
      dqbdlqm1(k) =
     .      + (0.5D0 + s(k)) * thet * (qn(k-1) + 0.5D0 * dqsdlq00(k-1))
     .      + (0.5D0 - s(k)) * thet * (        - 0.5D0 * dqsdlqm1(k  ))
      dqbdlqm2(k) =
     .        (0.5D0 + s(k)) * thet * (        + 0.5D0 * dqsdlqm1(k-1))
    3 continue
c
c-----------------------------------------------------------------------
c
      return
      end
