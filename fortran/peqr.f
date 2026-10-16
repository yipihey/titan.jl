      subroutine peqr
c
c***********************************************************************
c     prepare for NEWTON-RAPHSON iteration; 
c     time-dependent radiation in a static medium
c.......................................................................
c     calling sequence: peq > peqr > writout, clzfil, eos, opac, eddfac
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     replace solution at old time with solution at new time
c=======================================================================
c
      timeo = timen
      timen = timeo + dtime
c
      teso  = tes
      tewo  = tew
      telo  = tel
      teqo  = teq
c
      do 1 k = ngrs - 1, ngre + 2
      to    (k) = tn    (k)
      do    (k) = dn    (k)
      ro    (k) = rn    (k)
      rmuo  (k) = rmun  (k)
      rmup1o(k) = rmup1n(k)
      rmum1o(k) = rmum1n(k)
      dvolo (k) = dvoln (k)
      ego   (k) = egn   (k)
      ero   (k) = ern   (k)
      fro   (k) = frn   (k)
      drdto (k) = drdt  (k)
      dmdto (k) = dmdt  (k)
      xnto  (k) = xnt   (k)
    1 continue
c
      do 2 k = ngrs - 1, ngre + 1
      dso  (k) = dsn  (k)
      pgo  (k) = pgn  (k)
      aso  (k) = asn  (k)
      egso (k) = egsn (k)
      erso (k) = ersn (k)
      egro (k) = egrn (k)
      egrso(k) = egrsn(k)
      frso (k) = frsn (k)
      chifo(k) = chifn(k)
      xkeo (k) = xken (k)
      xkpo (k) = xkpn (k)
      xneo (k) = xnen (k)
      plfo (k) = plfn (k)
    2 continue
c
      do 3 k = ngrs, ngre + 1
      avchio(k) = avchin(k)
    3 continue
c
c=======================================================================
c     extrapolate radii and masses to new time level
c=======================================================================
c
      jext = 1
c
    4 if (leibc .gt. 0) drdt(ngrs  ) =   0.0D0
      if (leobc .gt. 0) drdt(ngre+1) =   0.0D0
      if (llibc .gt. 0) drdt(ngrs  ) =   0.0D0
      if (llobc .gt. 0) drdt(ngre+1) =   0.0D0
c
      do 5 k = ngrs, ngre + 1
      rn  (k) = ro (k) + drdt(k) * dtime
    5 continue
c
c-----------------------------------------------------------------------
c     check radius and mass for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre
      if ( rn(k+1) .le. rn(k) ) then
c
c          nonmonotonic, cut timestep and try again
c
           write(itty, 100) jext, dtime
           write(iout, 100) jext, dtime
	   write(iout, 101) k+1,rn(k+1),k,rn(k)
           tfac  = 0.5D0
           dtime = tfac  * dtime
           timen = timeo + dtime
c
               jext = jext + 1
           if (jext .le. next) go to 4
c
c          failure
c
           write(itty, 102) next
           write(iout, 102) next
           call writout
           call clzfil
           stop 'peq'
      end if
    6 continue
  100 format(' grid extrapolation nonmonotonic. jext ='i3
     .                                          ' dtime =' 1pe10.2)
  101 format( ' rn('i3') ='1pe22.14 ' rn('i3') ='e22.14)
  102 format(' grid extrapolation nonmonotonic after'i3' attempts')
c
c=======================================================================
c     boundary conditions
c=======================================================================
c-----------------------------------------------------------------------
c     radiation flux
c-----------------------------------------------------------------------
c
c     inner boundary; equations bc22 - bc25
c
      k = ngrs
c
      if (lribc .eq. 1) frn(k) = + (2.0D0 * geddl + 1.0D0) * cpi * xipl
     .                           -     cc * geddl * ern(k)
      if (lribc .eq. 2) frn(k) =    0.0D0
      if (lribc .eq. 3) frn(k) =    csigr * teff**4
     .                           * (1.0D0 - 0.5D0 * xmu)
     .                           +  0.5D0 *  xlum * xmu 
     .                           / (4.0D0*cpi*rn(k)**mu)
c
c     outer boundary; equations bc35 & bc36
c
      k = ngre + 1
c
      if (lrobc .eq. 1) frn(k) = - (2.0D0 * geddr + 1.0D0) * cpi * ximr 
     .                           +     cc * geddr * ern(k-1)
      if (lrobc .eq. 2) frn(k) =    0.0D0
      if (lrobc .eq. 3) frn(k) = frn(k-1) * (rn(k-1) / rn(k))**mu
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
c     zero flux eulerian; equations pz1, pz3, pz5, pz9 
c
      if (leibc .eq. 1) then
          tn (k) =   tn (k+1)
          dn (k) =   dn (k+1)
	  rn (k) =   rn (k+1) * 2.0D0 - rn (k+2)
      end if
c
c     nonzero flux eulerian; equations pz11, pz13, pz15, pz19 
c
      if (leibc .eq. 2) then
	  tn (k) =   tl
	  dn (k) =   dl
	  rn (k) =   rn (k+1) - delrl 
      end if
c
c     lagrangean; equations pz21, pz23, pz25, pz29
c
      if (llibc .eq. 1) then
          tn (k) =   tn (k+1)
	  dn (k) =   dn (k+1)
	  rn (k) = ( rn (k+1)**mup1 - xmup1 * delml / dn(k) )**rxm1
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
c     zero flux eulerian; equations pz39, pz41, pz43, pz47 
c
      if (leobc .eq. 1) then
	  tn (k  ) =   tn (k-1)
	  tn (k+1) =   tn (k  )
	  dn (k  ) =   dn (k-1)
	  dn (k+1) =   dn (k  )
	  rn (k+1) =   rn (k  ) * 2.0D0 - rn (k-1)
      end if
c
c     nonzero flux eulerian; equations pz49, pz51, pz53, pz57 
c
      if (leobc .eq. 2) then
	  tn (k  ) =   tr
	  tn (k+1) =   tr
	  dn (k  ) =   dr
	  dn (k+1) =   dr
	  rn (k+1) =   rn (k  ) + delrr
      end if
c
c     transmitting eulerian; equations pz59, pz61, pz63, pz67 
c
      if (leobc .eq. 3) then
	  tn (k  ) =   tn (k-1)
	  tn (k+1) =   tn (k  )
	  dn (k  ) =   dn (k-1)
	  dn (k+1) =   dn (k  )
	  rn (k+1) =   rn (k  ) + delrr
      end if
c
c     lagrangean; equations pz69, pz71, pz73, pz77
c
      if (llobc .gt. 0) then
	  tn (k  ) =   tn (k-1)
	  tn (k+1) =   tn (k  )
	  dn (k  ) =   dn (k-1)
	  dn (k+1) =   dn (k  )
	  rn (k+1) = ( rn (k  )**mup1 + xmup1 * delmr / dn(k) )**rxm1
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
c     compute new radius-related quantities
c=======================================================================
c
      do 7 k = ngrs - 1, ngre + 2
      rmun  (k) =  rn(k)**mu
      rmup1n(k) =  rn(k)**mup1
      rmum1n(k) =  rn(k)**mum1
    7 continue
c
      do 8 k = ngrs - 1, ngre + 1
      dvoln (k) = ( rmup1n(k+1) - rmup1n(k) ) / xmup1
    8 continue
c
c=======================================================================
c     eos, opacity, EDDINGTON factor
c=======================================================================
c
      call eos
      call opac
c
      do 9 k = ngrs - 1, ngre + 2
      egrn  (k) = egn(k) +  ( ern (k) / dn(k) )
      plfn  (k) = csigr/cpi * tn  (k)**4
      dplfdltn(k) =   4.0D0 * plfn(k)
    9 continue
c
      do 10 k = ngrs, ngre + 1
      avchin(k) = 2.0D0 / ( 1.0D0 / chifn(k-1) + 1.0D0 / chifn(k) )
   10 continue
c
      call eddfac
c
c=======================================================================
c     calculate time-centered quantities
c=======================================================================
c
      do 11 k = ngrs - 1, ngre + 2
      r    (k) = thet * rn   (k) + (1.D0 - thet) * ro   (k)
      fr   (k) = thet * frn  (k) + (1.D0 - thet) * fro  (k)
      drdt (k) =   ( rn(k) - ro(k) ) / dtime
      rmu  (k) = r(k)**mu
      rmup1(k) = r(k)**mup1
      rmum1(k) = r(k)**mum1
      frnom(k) = (1.0D0 - 0.5D0 * xmu) * csigr * teff**4
     .                  + 0.5D0 * xmu * xlum / (4.0D0 * cpi * rmun(k))
      frnom(k) = frn(k)              ! that's the choice for all tests
   11 continue
c
      do 12 k = ngrs - 1, ngre + 1
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
   12 continue
c
      do 13 k = ngrs, ngre + 1
      avchi (k) = 2.0D0 / ( 1.0D0 / chif (k-1) + 1.0D0 / chif (k) )
   13 continue
c
c=======================================================================
c     zero out turbulence- and hydro-related quantities
c=======================================================================
c
      do 14  k = ngrs - 1, ngre + 2
      dso   (k) = 0.0D0
      db    (k) = 0.0D0
      dsn   (k) = 0.0D0
      eto   (k) = 0.0D0
      et    (k) = 0.0D0
      etn   (k) = 0.0D0
      etso  (k) = 0.0D0
      etb   (k) = 0.0D0
      etsn  (k) = 0.0D0
      egso  (k) = 0.0D0
      egb   (k) = 0.0D0
      egsn  (k) = 0.0D0
      egto  (k) = 0.0D0
      egt   (k) = 0.0D0
      egtn  (k) = 0.0D0
      egtso (k) = 0.0D0
      egtb  (k) = 0.0D0
      egtsn (k) = 0.0D0
      egrto (k) = egro(k)
      egrt  (k) = egr (k)
      egrtn (k) = egrn(k)
      egrtso(k) = 0.0D0
      egrtb (k) = 0.0D0
      egrtsn(k) = 0.0D0
      uo    (k) = 0.0D0
      u     (k) = 0.0D0
      un    (k) = 0.0D0
      uso   (k) = 0.0D0
      ub    (k) = 0.0D0
      usn   (k) = 0.0D0
      unom  (k) = 0.0D0
      urel  (k) = 0.0D0
      bri   (k) = 4.0D0 * cpi * rmun(k) * frn(k)
   14 continue
c
      return
      end
