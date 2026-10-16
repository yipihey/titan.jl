      subroutine eddfac
c
c***********************************************************************
c     calculate EDDINGTON factors ef2
c.......................................................................
c     calling sequence: update{r,rh} > eddfac > eddf, eddg
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     use EDDINGTON approximation for diffusion
c=======================================================================
c
      if (ltran .gt. 1) then
c
	    do 200 k = ngrs - 1, ngre + 2
            fedd(k) = 1.0D0 / 3.0D0
  200       continue
c
	    geddl = 0.5D0
	    geddr = 0.5D0
c
	    return
c
      end if
c
c=======================================================================
c     planar geometry
c=======================================================================
c
c++++++++ mu .eq. 0 ++++++++++++++++++++++++++++++++++++++++++++++ begin
      if (mu .eq. 0) then
c
c-----------------------------------------------------------------------
c     source function and optical depth; equations ef3 & ef5
c-----------------------------------------------------------------------
c
      do 1 k = ngrs, ngre
      sf(k) = 
     . ((chifn(k) * dn(k) - csige * xnen(k))* plfn(k) 
     .                    + csige * xnen(k) * ern (k) * cc / (4.D0*cpi))
     .                                              / (chifn(k) * dn(k))
      dtau(k) = chifn(k) * dn(k) * (rn(k+1) - rn(k)) 
    1 continue
c
c.......................................................................
c     monotonize slopes of source function
c.......................................................................
c
      do 2 k = ngrs, ngre - 1
c
c     equation ef42
c
      dsf(k) = sf(k) - sf(k+1)
    2 continue
c
      do 3 k = ngrs + 1, ngre - 1
c
c     equation ef41
c
      dxm = 0.5D0 * (dtau(k-1) + dtau(k  ))
      dxp = 0.5D0 * (dtau(k  ) + dtau(k+1))
c
c     equation ef44
c
      sum   = dxm * dsf(k) + dxp * dsf(k-1)
      prod  = dsf(k-1) * dsf(k)
      denom = cvmgm( sum - epsedf, sum + epsedf, sum )
c
c     equation ef6
c
      bsf(k) = cvmgp( 2.0D0 * prod / denom, 0.0D0, prod)
      asf(k) = sf(k) - 0.5D0 * dtau(k) * bsf(k)
    3 continue
c
      bsf(ngrs) = bsf(ngrs + 1)
      bsf(ngre) = bsf(ngre - 1)
      asf(ngrs) = sf(ngrs) - 0.5D0 * dtau(ngrs) * bsf(ngrs)
      asf(ngre) = sf(ngre) - 0.5D0 * dtau(ngre) * bsf(ngre)
c
c-----------------------------------------------------------------------
c     for case of two optically reflecting boundaries must derive 
c     globally consistent incident intensity at inner boundary
c-----------------------------------------------------------------------
c
      if (lribc .ne. 2) go to 11
c
c.......................................................................
c     k = ngrs; equation ef13
c.......................................................................
c
      do 4 j = 2, nang
      cjk(j) =   0.0D0
      ejk(j) = - 1.0D0
      sjk(j) =   0.0D0
    4 continue
c
c.......................................................................
c     k = ngrs + 1, ... , ngre + 1; equations ef14, ef17 & ef18
c.......................................................................
c
      do 6 k = ngrs, ngre
      do 5 j = 2, nang
      dt = dtau(k) / ang(j)
      ex = exp(- dt)
      cjk(j) = - ex 
      ejk(j) = - cjk(j) * ejk(j)
      sjk(j) = - cjk(j) * sjk(j)
     .         + asf(k) * (1.0D0 - ex) 
     .         + bsf(k) * ang(j) * (1.0D0 - (1.0D0 + dt) * ex) 
c
    5 continue
    6 continue
c
c.......................................................................
c     k = ngre + 2; equations ef15
c.......................................................................
c
      do 7 j = 2, nang
      cjk(j) = - 1.0D0
      ejk(j) = - cjk(j) * ejk(j)
      sjk(j) = - cjk(j) * sjk(j)
    7 continue
c
c.......................................................................
c     k = ngre + 3,..., 2 * ngre + 3 - ngrs. (take k = 2 * ngre + 3 - kp
c     where kp = ngre, ... , ngrs); equation ef16
c.......................................................................
c
      do 9 kp = ngre, ngrs, - 1
      do 8 j = 2, nang
c
      dt = dtau(k) / ang(j)
      ex = exp(- dt)
      cjk(j) = - ex 
      ejk(j) = - cjk(j) * ejk(j)
      sjk(j) = - cjk(j) * sjk(j)
     .         + asf(k) * (1.0D0 - ex) 
     .         + bsf(k) * ang(j) * (dt - 1.0D0 + ex) 
c
    8 continue
    9 continue
c
c.......................................................................
c     finally solve for incoming intensity at inner boundary; by 
c     reflection equals outgoing intensity needed for boundary condition
c     equations ef19 - ef21
c.......................................................................
c
      do 10 j = 2, nang
      xim(j) = sjk(j) / ( 1.0D0 + ejk(j) )
   10 continue
c
c     handle singularity at mu = 0. 
      xim(1) = xim(2)
c
c     proceed with standard planar solution.
c
c-----------------------------------------------------------------------
c     outgoing rays
c-----------------------------------------------------------------------
c.......................................................................
c     inner boundary condition
c.......................................................................
c
   11 if (lribc .eq. 1) then
c
c           equation ef9
c
	    do 12 j = 1, nang
            xim(j) = xipl
   12       continue
c
      end if
c
c     if (lribc .eq. 2) incident xim(j) was determined from globally
c                       consistent solution above.
c
      if (lribc .eq. 3) then
c
c           equation ef10
c
	    do 13 j = 1, nang
	    xim(j) = ( cc * ern(ngrs) + 1.5D0 * frn(ngrs) * dtau(ngrs) ) 
     .                                                    / (4.D0 * cpi)
     .                        + 3.D0 * ang(j) * frn(ngrs) / (4.D0 * cpi)
   13       continue
c
      end if
c
c.......................................................................
c     outward sweep; equation ef7
c.......................................................................
c
c     integrate first to midpoint then to next interface
c
      do 16 k = ngrs, ngre
c
      do 14 j = 2, nang
      dt = 0.5D0 * dtau(k) / ang(j)
      ex = exp(- dt)
      xi0(j) =  sf(k) * (1.0D0 - ex) 
     .       + bsf(k) * ang(j) * (1.0D0 - (1.0D0 + dt) * ex) 
     .       + xim(j) * ex
c
      dt = 2.0D0 * dt
      ex = exp(- dt)
      xip(j) = asf(k) * (1.0D0 - ex) 
     .       + bsf(k) * ang(j) * (1.0D0 - (1.0D0 + dt) * ex) 
     .       + xim(j) * ex
c
   14 continue
c
c     handle singularity at mu = 0.
c
      xi0(1) = xi0(2)
      xip(1) = xip(2)
c
c    swap intensity on outer boundary of cell k to inner boundary of k+1
c
      do 15 j = 1, nang
      xim(j) = xip(j)
   15 continue
c
c.......................................................................
c     angle moments; equations ef48 & ef50
c.......................................................................
c
c     contribution of outgoing radiation to zeroth and second moments
c
      xmom0(k) = 0.0D0
      xmom2(k) = 0.0D0
      call eddf (nang, ang, xi0, xmom0(k), xmom2(k))
c
c     flux EDDINGTON factor at outer boundary
c
      if (k .eq. ngre) call eddg (nang, ang, xip, geddr)
c
   16 continue
c
c-----------------------------------------------------------------------
c     incoming rays
c-----------------------------------------------------------------------
c.......................................................................
c     outer boundary condition
c.......................................................................
c
      if (lrobc .eq. 1) then
c
c           equation ef11
c
	    do 17 j = 1, nang
	    xip(j) = ximr
   17       continue
      end if
c
c     if (lrobc .eq. 2) xip(j) already is (reflected) outgoing intensity
c
c.......................................................................
c     inward sweep; equation ef8
c.......................................................................
c
c     integrate first to midpoint then to next interface
c
      do 20 k = ngre, ngrs, - 1
c
      do 18 j = 2, nang
      dt = 0.5D0 * dtau(k) / ang(j)
      ex = exp(- dt)
      xi0(j) = asf(k) * (1.0D0 - ex) 
     .       + bsf(k) * ang(j) * (dt - 1.0D0 + ex) 
     .       + xip(j) * ex
c
      dt = 2.0D0 * dt
      ex = exp(- dt)
      xim(j) = asf(k) * (1.0D0 - ex) 
     .       + bsf(k) * ang(j) * (dt - 1.0D0 + ex) 
     .       + xip(j) * ex
c
   18 continue
c
c     handle singularity at mu = 0.
c
      xi0(1) = xi0(2)
      xim(1) = xim(2)
c
c    swap intensity on inner boundary of cell k to outer boundary of k-1
c
      do 19 j = 1, nang
      xip(j) = xim(j)
   19 continue
c
c.......................................................................
c     angle moments; equations ef48 & ef50
c.......................................................................
c
c     contribution of incoming radiation to zeroth and second moments
c
      call eddf (nang, ang, xi0, xmom0(k), xmom2(k))
c
c     compute EDDINGTON factor; equation ef51
c
      fedd(k) = xmom2(k) / xmom0(k)
c
      fedd(k) = max( 0.0D0, min( fedd(k), 1.0D0 ) )
c
c     flux EDDINGTON factor at inner boundary
c
      if (k .eq. ngrs) call eddg (nang, ang, xim, geddl)
c
      geddl = max( 0.5D0, min( geddl, 1.0D0 ) )
c
   20 continue
c
      return
c
c++++++++ mu .eq. 2 ++++++++++++++++++++++++++++++++++++++++++++++++ end
c
      end if
c
c=======================================================================
c     spherical geometry
c=======================================================================
c-----------------------------------------------------------------------
c     source function and optical depth; equations ef3 & ef5
c-----------------------------------------------------------------------
c
c++++++++ mu .eq. 2 ++++++++++++++++++++++++++++++++++++++++++++++ begin
      if (mu .eq. 2) then
c
c     equations ef3 & ef5
c
      do 21 k = ngre, ngrs, -1
      sf(k) = 
     . ((chifn(k) * dn(k) - csige * xnen(k))* plfn(k) 
     .                    + csige * xnen(k) * ern (k) * cc / (4.D0*cpi))
     .                                              / (chifn(k) * dn(k))
      dtau(k) = chifn(k) * dn(k) * (rn(k+1) - rn(k)) 
   21 continue
c
c.......................................................................
c     monotonize slopes of source function
c.......................................................................
c
      do 22 k = ngrs, ngre - 1
c
c     equation ef42
c
      dsf(k) = sf(k) - sf(k+1)
   22 continue
c
      do 23 k = ngrs + 1, ngre - 1
c
c     equation ef43
c
      dxm = 0.5D0 * (dtau(k-1) + dtau(k  ))
      dxp = 0.5D0 * (dtau(k  ) + dtau(k+1))
c
c     equation ef44
c
      sum   = dxm * dsf(k) + dxp * dsf(k-1)
      prod  = dsf(k-1) * dsf(k)
      denom = cvmgm( sum - epsedf, sum + epsedf, sum )
c
c     equation ef6
c
      bsf(k) = cvmgp( 2.0D0 * prod / denom, 0.0D0, prod)
   23 continue
c
      bsf(ngrs) = bsf(ngrs + 1)
      bsf(ngre) = bsf(ngre - 1)
c    
c-----------------------------------------------------------------------
c     define augmented radial mesh, interface values of source function
c-----------------------------------------------------------------------
c
      do 24 k = ngrs, ngre 
c
      l = 2 * (k - ngrs) + 1
      rp(l  ) = rn(k)
      rp(l+1) = 0.5D0 * (rn(k) + rn(k+1))
c
c     equations ef22 - ef25
c
      sl(l  ) = sf(k) + 0.5D0 * dtau(k) * bsf(k)
      sr(l  ) = sf(k)
      sl(l+1) = sf(k)
      sr(l+1) = sf(k) - 0.5D0 * dtau(k) * bsf(k)
c
   24 continue
c    
      lmax = 2 * (ngre - ngrs + 1) + 1
      rp(lmax) = rn(ngre+1)
c    
c-----------------------------------------------------------------------
c     define impact parameters
c-----------------------------------------------------------------------
c
c     equation ef27
c
      do 25 j = 1, ncor
      p(j) = float(j - 1) * rp(1) / float(ncor)
   25 continue
c
c     equation ef29
c
      nray = (ncor + 1) + 2 * (ngre + 1 - ngrs)
c
c     equation ef28
c
      do 26 j = ncor + 1, nray
      p(j) = rp(j - ncor)
   26 continue
c
c-----------------------------------------------------------------------
c     incoming rays
c-----------------------------------------------------------------------
c.......................................................................
c     outer boundary condition; equation ef11
c.......................................................................
c
      if (lrobc .eq. 1) then
	    do 27 j = 1, nray
	    xim(j) = ximr
   27       continue
c
      end if
c
c.......................................................................
c     inward sweep
c.......................................................................
c
      do 30 l = lmax - 1, 1, - 1
      k  = ngrs + (l - 1) / 2
      nj = ncor + l
c
      do 28 j = 1, nj
c
c     equation ef30
c
      xjl   = sqrt( rp(l  )**2 - p(j)**2 )
c
c     equation ef31
c
      xjlp1 = sqrt( rp(l+1)**2 - p(j)**2 )
c
c     equation ef32
c
      dtjl  = dn(k) * chifn(k) * (xjlp1 - xjl)
c
c     equations ef33 & ef34
c
      ex = exp(- dtjl)
      xip (j) = (1.0D0 - ex) * sl(l) 
     .        + (1.0D0 - (1.0D0 + dtjl) * ex) * (sr(l) - sl(l)) / dtjl
     .        + ex * xim(j)
c
c     equation ef35
c
      xang(j) = sqrt(1.0D0 - p(j)**2 / rp(l)**2)
   28 continue
c
c     swap intensity on inner boundary of current cell to outer boundary
c     of next cell
c
      do 29 j = 1, nj
      xim(j) = xip(j)
   29 continue
c
c.......................................................................
c     angle moments; equations ef48 & ef50
c.......................................................................
c
c     contribution of incoming radiation to zeroth and second moments
c
      if (mod(l,2) .eq. 0) then
	    k = ngrs + (l - 1) / 2
	    xmom0(k) = 0.0D0
	    xmom2(k) = 0.0D0
	    call eddf (nj, xang, xip, xmom0(k), xmom2(k))
      end if
c
c     flux EDDINGTON factor at inner boundary
c
      if (l .eq. 1) call eddg (nj, xang, xip, geddl)
c
   30 continue
c
c-----------------------------------------------------------------------
c     outgoing rays
c-----------------------------------------------------------------------
c.......................................................................
c     inner boundary condition
c.......................................................................
c
      if (lribc .eq. 1) then
c
c     equation ef39
c
	    l = 1
	    do 31 j = 1, ncor + 1
	    xang(j) = sqrt( 1.0D0 - p(j)**2 / rp(l)**2)
            xim (j) = xipl
   31       continue
c
      end if
c
      if (lribc .eq. 3) then
c
c     equation ef40
c
	    l = 1
	    do 32 j = 1, ncor + 1
	    xang(j) = sqrt( 1.0D0 - p(j)**2 / rp(l)**2)
	    xim(j) = ( cc * ern(ngrs) + 1.5D0 * frn(ngrs) * dtau(ngrs) ) 
     .                                                    / (4.D0 * cpi)
     .                        + 3.D0 * ang(j) * frn(ngrs) / (4.D0 * cpi)
   32       continue
c
      end if
c
c.......................................................................
c     outward sweep
c.......................................................................
c
      do 35 l = 2, lmax
      k  = ngrs + (l - 2) / 2
      nj = ncor +  l - 1
c
      do 33 j = 1, nj
c
c     equation ef30
c
      xjl   = sqrt( rp(l  )**2 - p(j)**2 )
c
c     equation ef31
c
      xjlm1 = sqrt( rp(l-1)**2 - p(j)**2 )
c
c     equation ef32
c
      dtjl  = dn(k) * chifn(k) * (xjl - xjlm1)
c
c     equations ef37 & ef38
c
      ex   = exp(- dtjl)
      xip (j) = (1.D0 - ex) * sr(l-1) 
     .        + (1.D0 - (1.D0 + dtjl) * ex) * (sl(l-1) - sr(l-1)) / dtjl
     .        + ex * xim(j)
c
c     equation ef35
c
      xang(j) = sqrt( 1.0D0 - p(j)**2 / rp(l)**2)
c
   33 continue
c
c     swap intensity on outer boundary of cell (k,k+1) to inner boundary
c     of next cell
c
      do 34 j = 1, nj
      xim(j) = xip(j)
   34 continue
c
c.......................................................................
c     angle moments; equations ef48 & ef50
c.......................................................................
c
c     contribution of outgoing radiation to zeroth and second moments
c
      if (mod(l,2) .eq. 0) then
	    k = ngrs + (l - 1) / 2
	    xang(nj+1) = 0.0D0
	    call eddf (nj + 1, xang, xim, xmom0(k), xmom2(k))
c
c           compute EDDINGTON factor; equation ef51
c
	    fedd(k) = xmom2(k) / xmom0(k)
c
            fedd(k) = max( 0.0D0, min( fedd(k), 1.0D0 ) )
      end if
c
c     flux EDDINGTON factor at outer boundary
c
      if (l .eq. lmax) call eddg (nj, xang, xip, geddr)
c
      geddr = max( 0.5D0, min( geddr, 1.0D0 ) )
c
   35 continue
c
      return
c
c++++++++ mu .eq. 2 ++++++++++++++++++++++++++++++++++++++++++++++++ end
c
      end if
c
      end
