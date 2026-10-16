      subroutine diffuse
c
c***********************************************************************
c     artificial mass and energy diffusion
c.......................................................................
c     called by: contin, gasnrg{h,rh}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.eq3'
c
      if (sig .eq. 0.0D0) return
      xtht = xmu * thet
c
c=======================================================================
c     mass diffusion
c=======================================================================
c
      if (zet .eq. 0.0D0) then
c
            do 1 k = ngrs, ngre + 1
c
c           equation c42
c
             df     (k) =    2.0D0 * sig *  rmu(k) * ( d(k  ) - d(k-1) )
     .                                             / ( r(k+1) - r(k-1) )
c
c           equations c43 - c45
c
            ddfdlrp1(k) = - df(k) * thet * rn(k+1) / ( r(k+1) - r(k-1) )
            ddfdlr00(k) =   df(k) * xtht * rn(k  ) /   r(k  )
            ddfdlrm1(k) = + df(k) * thet * rn(k-1) / ( r(k+1) - r(k-1) )
c
c           equations c46 & c47
c
            ddfdld00(k) = + 2.0D0 * sig * ( rmu(k) / (r(k+1) - r(k-1)) )
     .                            * thet * dn(k  )
            ddfdldm1(k) = - 2.0D0 * sig * ( rmu(k) / (r(k+1) - r(k-1)) )
     .                            * thet * dn(k-1)
    1       continue
            return
      endif
c
c=======================================================================
c     energy diffusion
c=======================================================================
c
      if (zet .eq. 1.0D0) then
c
            do 2 k = ngrs, ngre + 1
c
c           equation c42
c
             df     (k) =           sig * rmu(k) * ( d (k  ) + d (k-1) )
     .                                           * ( eg(k  ) - eg(k-1) ) 
     .                                           / ( r (k+1) - r (k-1) )
c
c           equations c43 - c45
c
            ddfdlrp1(k) = - df(k) * thet * rn(k+1) / ( r(k+1) - r(k-1) )
            ddfdlr00(k) =   df(k) * xtht * rn(k  ) /   r(k  )
            ddfdlrm1(k) = + df(k) * thet * rn(k-1) / ( r(k+1) - r(k-1) )
c
c           equations c46 & c47
c
            ddfdle00    = +  thet * sig * rmu(k  ) / ( r(k+1) - r(k-1) )
     .                                  * egn(k  ) * ( d(k  ) + d(k-1) )
            ddfdlem1    = -  thet * sig * rmu(k  ) / ( r(k+1) - r(k-1) )
     .                                  * egn(k-1) * ( d(k  ) + d(k-1) )
c
c           equations c48 & c49
c
            ddfdld00(k) =   df(k) * thet * dn(k  ) / ( d(k  ) + d(k-1) ) 
     .                  + ddfdle00 * dlegdldn(k  )
            ddfdldm1(k) =   df(k) * thet * dn(k-1) / ( d(k  ) + d(k-1) )
     .                  + ddfdlem1 * dlegdldn(k-1)
c
c           equations fe58 & fe59
c
            ddfdlt00(k) = ddfdle00 * dlegdltn(k  )
            ddfdltm1(k) = ddfdlem1 * dlegdltn(k-1)
    2       continue
            return
      endif
c
c-----------------------------------------------------------------------
c
      stop 'diffuse'
      end
