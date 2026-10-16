      subroutine viscous
c
c***********************************************************************
c     viscous energy dissipation and momentum deposition
c.......................................................................
c     called by: matgen
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     define various quantities and their derivatives
c=======================================================================
c-----------------------------------------------------------------------
c     f1 = du/dr - mu/2 <u/r>; equations fe60 - fe64
c-----------------------------------------------------------------------
c
      f1  (xk,xp,uk,up) =  (up-uk)/(xp-xk)    - xmu4*(up/xp   +uk/xk   )
      f1rk(xk,xp,uk,up) = +(up-uk)/(xp-xk)**2 + xmu4*(        +uk/xk**2)
      f1rp(xk,xp,uk,up) = -(up-uk)/(xp-xk)**2 + xmu4*(up/xp**2         )
      f1uk(xk,xp,uk,up) =  (  -en)/(xp-xk)    - xmu4*(        +en/xk   )
      f1up(xk,xp,uk,up) =  (en   )/(xp-xk)    - xmu4*(en/xp            )
c
c=======================================================================
c
      xmu4 = 0.25D0 * xmu
      en   = 1.00D0
c
c-----------------------------------------------------------------------
c     artificial viscosity coefficient
c-----------------------------------------------------------------------
c
      do 1 k = ngrs, ngre
c
c     equations fe65 - fe69
c
       dudr   (k) = f1   (r(k),r(k+1),u(k),u(k+1))
      durdlrp1(k) = f1rp (r(k),r(k+1),u(k),u(k+1)) * thet * rn  (k+1)
      durdlr00(k) = f1rk (r(k),r(k+1),u(k),u(k+1)) * thet * rn  (k  )
      durdlup1(k) = f1up (r(k),r(k+1),u(k),u(k+1)) * thet * unom(k+1)
      durdlu00(k) = f1uk (r(k),r(k+1),u(k),u(k+1)) * thet * unom(k  )
c
c     dissipation length; equations fe70 - fe72
c
       ql  (k) = ql0 + ql1 * 0.5D0 * (r (k) + r (k+1))
      dqldlrp1 =       ql1 * 0.5D0 * (        rn(k+1)) * thet 
      dqldlr00 =       ql1 * 0.5D0 * (rn(k)          ) * thet 
c
c     velocity divergence; equations fe73 - fe77
c
       div (k) =         ( rmu  (k+1) * u(k+1) 
     .                   - rmu  (k  ) * u(k  ) )          / dvol(k)
      ddivdrp1 = + ( xmu * rmum1(k+1) * u(k+1) - rmu(k+1) *  div(k) ) 
     .                                                    / dvol(k)
      ddivdr00 = - ( xmu * rmum1(k  ) * u(k  ) - rmu(k  ) *  div(k) )
     .                                                    / dvol(k)
      ddivdup1 = +         rmu  (k+1) * 1.0D0             / dvol(k) 
      ddivdu00 = -         rmu  (k  ) * 1.0D0             / dvol(k) 
c
c     equations fe78 - fe82
c
       qv     (k) =   min(  div (k)                   , 0.D0)
      dqvdlrp1(k) = cvmgm( ddivdrp1 * thet * rn  (k+1), 0.D0, div(k) )
      dqvdlr00(k) = cvmgm( ddivdr00 * thet * rn  (k  ), 0.D0, div(k) )
      dqvdlup1(k) = cvmgm( ddivdup1 * thet * unom(k+1), 0.D0, div(k) )
      dqvdlu00(k) = cvmgm( ddivdu00 * thet * unom(k  ), 0.D0, div(k) )
c
c     viscosity coefficient; equations fe83 - fe89
c
       qm     (k) =    ql  (k) * ( + cq1 * as(k) 
     .                             - cq2 * ql(k)   *  qv     (k) )
      dqmdlrp1(k) =   dqldlrp1 * ( + cq1 * as(k) 
     .                             - cq2 * ql(k)*2.D0*qv     (k) )
     .                             - cq2 * ql(k)**2* dqvdlrp1(k) 
      dqmdlr00(k) =   dqldlr00 * ( + cq1 * as(k) 
     .                             - cq2 * ql(k)*2.D0*qv     (k) )
     .                             - cq2 * ql(k)**2* dqvdlr00(k) 
      dqmdlup1(k) =                - cq2 * ql(k)**2* dqvdlup1(k)
      dqmdlu00(k) =                - cq2 * ql(k)**2* dqvdlu00(k)
      dqmdld00(k) =    ql     (k)  * cq1 * as(k)   * thet * 0.5D0
     .         * ( - dn(k) / d(k)  + pgn(k) * dlpgdldn(k) / pg(k) ) 
      dqmdlt00(k) =    ql     (k)  * cq1 * as(k)   * thet * 0.5D0
     .         * (                   pgn(k) * dlpgdltn(k) / pg(k) )
c
c     viscous energy dissipaton rate; equation fe90
c
      qe(k) = - cqvis * d(k) * qm(k) * dudr(k)**2 * dvol(k)
c
c     viscous pressure; equations gm48 - gm54
c
       qf     (k) = - cqvis * d (k) *    qm     (k) *  dudr(k)
      dqfdlrp1(k) = - cqvis * d (k) * ( dqmdlrp1(k) *  dudr(k)
     .                                +  qm     (k) * durdlrp1(k) )
      dqfdlr00(k) = - cqvis * d (k) * ( dqmdlr00(k) *  dudr(k)
     .                                +  qm     (k) * durdlr00(k) )
      dqfdlup1(k) = - cqvis * d (k) * ( dqmdlup1(k) *  dudr(k)
     .                                +  qm     (k) * durdlup1(k) )
      dqfdlu00(k) = - cqvis * d (k) * ( dqmdlu00(k) *  dudr(k)
     .                                +  qm     (k) * durdlu00(k) )
      dqfdlt00(k) = - cqvis * d (k) *   dqmdlt00(k) *  dudr(k)
      dqfdld00(k) = - cqvis *(d (k) *   dqmdld00(k) *  dudr(k)
     .               + thet * dn(k) *    qm     (k) *  dudr(k)    )
c
c     "kinetic pressure"; equations ag88 - ag92
c
       pk     (k) =       pg (k) 
     .                             + 0.25D0 * d (k) * (u(k) + u(k+1))**2
      dpkdlt00(k) = thet* pgn(k) * dlpgdltn(k)
      dpkdld00(k) = thet* pgn(k) * dlpgdldn(k)
     .            + thet           * 0.25D0 * dn(k) * (u(k) + u(k+1))**2
      dpkdlu00(k) = thet* unom(k  )* 0.50D0 * d (k) * (u(k) + u(k+1))
      dpkdlup1(k) = thet* unom(k+1)* 0.50D0 * d (k) * (u(k) + u(k+1))
c
c     viscosity indicator for adaptive grid;
c     equations ag87 and ag93 - ag98
c
       qx     (k) =  qf     (k) / pk(k) + q0
      dqxdlrp1(k) = dqfdlrp1(k) / pk(k)
      dqxdlr00(k) = dqfdlr00(k) / pk(k)
      dqxdlup1(k) = dqfdlup1(k) / pk(k) - qf(k) * dpkdlup1(k) / pk(k)**2
      dqxdlu00(k) = dqfdlu00(k) / pk(k) - qf(k) * dpkdlu00(k) / pk(k)**2
      dqxdlt00(k) = dqfdlt00(k) / pk(k) - qf(k) * dpkdlt00(k) / pk(k)**2
      dqxdld00(k) = dqfdld00(k) / pk(k) - qf(k) * dpkdld00(k) / pk(k)**2
c
       qv     (k) = 0.0D0 
    1 continue
c
c=======================================================================
c
      return
      end
