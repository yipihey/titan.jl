      subroutine updaterh
c
c***********************************************************************
c     update solution at advanced time level after NEWTON-RAPHSON cycle
c     for radiation hydrodynamics, no turbulence
c.......................................................................
c     calling sequence: update > updaterh > eos, opac, eddfac
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     apply the computed changes
c=======================================================================
c
      do 11 k = ngrs, ngre + 1
      rn  (k) = (1.0D0 + rhs(ir, k)) * rn   (k)
      xmen(k) = (1.0D0 + rhs(im, k)) * xmen (k)
      un  (k) = un (k) + rhs(iu, k)  *  unom(k)
      frn (k) = frn(k) + rhs(if, k)  * frnom(k)
   11 continue
c
      do 12 k = ngrs, ngre
      dn  (k) = (1.0D0 + rhs(id, k)) * dn (k)
      tn  (k) = (1.0D0 + rhs(it, k)) * tn (k)
      ern (k) = (1.0D0 + rhs(ie, k)) * ern(k)
   12 continue
c
c=======================================================================
c     phantom zones
c=======================================================================
c-----------------------------------------------------------------------
c     inner boundary
c-----------------------------------------------------------------------
c
      k = ngrs - 1
c
c     zero flux eulerian; equations pz1, pz3, pz5, pz7, pz9
c
      if (leibc .eq. 1) then
          tn  (k) =   tn  (k+1)
          dn  (k) =   dn  (k+1)
	  rn  (k) =   rn  (k+1) * 2.0D0 - rn  (k+2)
	  xmen(k) =   xmen(k+1) * 2.0D0 - xmen(k+2)
          un  (k) = - un  (k+2)
      end if
c
c     nonzero flux eulerian; equations pz11, pz13, pz15, pz17, pz19
c
      if (leibc .eq. 2) then
	  tn  (k) =   tl
	  dn  (k) =   dl
	  rn  (k) =   rn  (k+1) - delrl 
	  xmen(k) =   xmen(k+1) + dl   *rxm1*(rn(k+1)**mup1-rn(k)**mup1)
          un  (k) =   ul
      end if
c
c     lagrangean; equations pz21, pz23, pz25, pz27, pz29
c
      if (llibc .eq. 1) then
          tn  (k) =   tn  (k+1)
	  dn  (k) =   dn  (k+1)
          rn  (k) = ( rn  (k+1)**mup1 - xmup1 * delml / dn(k) )**rxm1
	  xmen(k) =   xmen(k+1)               + delml
          un  (k) =   un  (k+1)
      end if
c
c.......................................................................
c     radiation
c.......................................................................
c
c                             equation pz31
c
	                ern(k) =   ern(k+1)
c
c     optically transmitting; equation pz33 
c
      if (lribc .eq. 1) frn(k) = + frn(k+1)
c
c     optically reflecting; equation pz35
c
      if (lribc .eq. 2) frn(k) = - frn(k+2)
c
c     imposed flux; equation pz37
c
      if (lribc .eq. 3) frn(k) =   frn(k+1) * (rn(k+1) / rn(k))**mu
c
c-----------------------------------------------------------------------
c     outer boundary
c-----------------------------------------------------------------------
c
      k = ngre + 1
c
c     zero flux eulerian; equations pz39, pz41, pz43, pz45, pz47
c
      if (leobc .eq. 1) then
	  tn  (k  ) =   tn  (k-1)
	  tn  (k+1) =   tn  (k  )
	  dn  (k  ) =   dn  (k-1)
	  dn  (k+1) =   dn  (k  )
	  rn  (k+1) =   rn  (k  ) * 2.0D0 - rn  (k-1)
          xmen(k+1) =   xmen(k  ) * 2.0D0 - xmen(k-1)
          un  (k+1) = - un  (k-1)
      end if
c
c     nonzero flux eulerian; equations pz49, pz51, pz53, pz55, pz57
c
      if (leobc .eq. 2) then
	  tn  (k  ) = tr
	  tn  (k+1) = tr
	  dn  (k  ) = dr
	  dn  (k+1) = dr
	  rn  (k+1) = rn  (k  ) + delrr
	  xmen(k+1) = xmen(k  ) - dr   *rxm1*(rn(k+1)**mup1-rn(k)**mup1)
          un  (k+1) = ur
      end if
c
c     transmitting eulerian; equations pz59, pz61, pz63, pz65, pz67
c
      if (leobc .eq. 3) then
	  tn  (k  ) = tn  (k-1)
	  tn  (k+1) = tn  (k  )
	  dn  (k  ) = dn  (k-1)
	  dn  (k+1) = dn  (k  )
	  rn  (k+1) = rn  (k  ) + delrr
	  xmen(k+1) = xmen(k  ) - dn(k)*rxm1*(rn(k+1)**mup1-rn(k)**mup1)
          un  (k+1) = un  (k  )
      end if
c
c     lagrangean; equations pz69, pz71, pz73, pz75, pz77
c
      if (llobc .gt. 0) then
	  tn  (k  ) =   tn  (k-1)
	  tn  (k+1) =   tn  (k  )
	  dn  (k  ) =   dn  (k-1)
	  dn  (k+1) =   dn  (k  )
	  rn  (k+1) = ( rn  (k  )**mup1 + xmup1 * delmr / dn(k) )**rxm1
	  xmen(k+1) =   xmen(k  )               - delmr
          un  (k+1) =   un  (k  )
      end if
c
c.......................................................................
c     radiation
c.......................................................................
c
c     optically transmitting; equations pz79 & pz81
c
      if (lrobc .eq. 1) then
          ern(k  ) =   ern(k-1)
          ern(k+1) =   ern(k  )
          frn(k+1) = + frn(k  )
      end if
c
c     optically reflecting; equations pz83 & pz85
c
      if (lrobc .eq. 2) then
          ern(k  ) =   ern(k-1)
          ern(k+1) =   ern(k-2)
          frn(k+1) = - frn(k-1)
      end if
c
c     imposed net flux; equations pz87 & pz89
c
      if (lrobc .eq. 3) then
          ern(k  ) =   ern(k-1)
          ern(k+1) =   ern(k  )
          frn(k+1) =   frn(k  ) * (rn(k) / rn(k+1))**mu
      end if
c
c=======================================================================
c     update radius and mass related quantities
c=======================================================================
c
      do 13 k = ngrs - 1, ngre + 2
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
      xmn   (k) = xmtot - xmen(k)
   13 continue
c
      do 14 k = ngrs - 1, ngre + 1
      dvoln(k) = ( rmup1n(k+1) - rmup1n(k) ) / xmup1
   14 continue
c
c=======================================================================
c     update eos, opacity, EDDINGTON factor
c=======================================================================
c
      call eos
      call opac
c
      do 15 k  = ngrs - 1, ngre + 2
      egrn (k) =   egn(k) + ( ern(k) / dn(k) )
      plfn (k) =   csigr  * tn(k)**4 / cpi
      dplfdltn(k) = 4.0D0 * plfn(k)
   15 continue
c
      do 16 k = ngrs, ngre + 1
      avchin(k) = 2.0D0 / ( 1.0D0 / chifn(k-1) + 1.0D0 / chifn(k) )
   16 continue
c
      if (ltran .eq. 0) call eddfac
c
c=======================================================================
c     calculate time-centered quantities
c=======================================================================
c
      do 17 k = ngrs - 1, ngre + 2
      r    (k) = thet * rn   (k) + (1.D0 - thet) * ro   (k)
      u    (k) = thet * un   (k) + (1.D0 - thet) * uo   (k)
      xm   (k) = thet * xmn  (k) + (1.D0 - thet) * xmo  (k)
      xme  (k) = thet * xmen (k) + (1.D0 - thet) * xmeo (k)
      fr   (k) = thet * frn  (k) + (1.D0 - thet) * fro  (k)
      dmdt (k) = - ( xmen(k) - xmeo(k) ) / dtime
      drdt (k) =   ( rn  (k) - ro  (k) ) / dtime
      urel (k) =     u   (k) - drdt(k)
      rmu  (k) = r(k)**mu
      rmup1(k) = r(k)**mup1
      rmum1(k) = r(k)**mum1
      frnom(k) = (1.0D0 - 0.5D0 * xmu) * csigr * teff**4
     .                  + 0.5D0 * xmu  * xlum  / (4.0D0 * cpi * rmun(k))
c                                                        ! equation re10
      frnom(k) = abs(frn(k)) + 1.0D-30 ! that's the choice for all tests
       unom(k) = abs(un (k)) + 1.0D-30
   17 continue
c
      do 18 k = ngrs - 1, ngre + 1
      d    (k) = thet * dn   (k) + (1.D0 - thet) * do   (k)
      t    (k) = thet * tn   (k) + (1.D0 - thet) * to   (k)
      eg   (k) = thet * egn  (k) + (1.D0 - thet) * ego  (k)
      pg   (k) = thet * pgn  (k) + (1.D0 - thet) * pgo  (k)
      as   (k) = thet * asn  (k) + (1.D0 - thet) * aso  (k)
      er   (k) = thet * ern  (k) + (1.D0 - thet) * ero  (k)
      egr  (k) = thet * egrn (k) + (1.D0 - thet) * egro (k)
      chif (k) = thet * chifn(k) + (1.D0 - thet) * chifo(k)
      xke  (k) = thet * xken (k) + (1.D0 - thet) * xkeo (k)
      xkp  (k) = thet * xkpn (k) + (1.D0 - thet) * xkpo (k)
      plf  (k) = thet * plfn (k) + (1.D0 - thet) * plfo (k)
      dvol (k) =  ( rmup1(k+1) - rmup1(k) ) / xmup1
   18 continue
c
      do 19 k = ngrs, ngre + 1
      avchi(k) = 2.0D0 / ( 1.0D0 / chif(k-1) + 1.0D0 / chif(k) )
   19 continue
c
c=======================================================================
c     zero out turbulence-related quantities
c=======================================================================
c
      do 20 k = ngrs - 1, ngre + 2
      eto   (k) = 0.0D0
      et    (k) = 0.0D0
      etn   (k) = 0.0D0
      etso  (k) = 0.0D0
      etb   (k) = 0.0D0
      etsn  (k) = 0.0D0
      egso  (k) = 0.0D0
      egb   (k) = 0.0D0
      egsn  (k) = 0.0D0
      egrto (k) = egro(k)
      egrt  (k) = egr (k)
      egrtn (k) = egrn(k)
      egrtso(k) = 0.0D0
      egrtb (k) = 0.0D0
      egrtsn(k) = 0.0D0
      bri   (k) = unom(k)
   20 continue
c
      return
      end
