      subroutine peqh
c
c***********************************************************************
c     prepare for NEWTON-RAPHSON iteration; 
c     hydrodynamics only, no radiation or turbulence
c.......................................................................
c     calling sequence: peq > peqh > writout, clzfil, eos
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
      xmeo  (k) = xmen  (k)
      xmo   (k) = xmn   (k)
      uo    (k) = un    (k)
      ego   (k) = egn   (k)
      drdto (k) = drdt  (k)
      dmdto (k) = dmdt  (k)
      urelo (k) = urel  (k)
      xnto  (k) = xnt   (k)
    1 continue
c
      do 2 k = ngrs - 1, ngre + 1
      dso (k) = dsn (k)
      uso (k) = usn (k)
      pgo (k) = pgn (k)
      aso (k) = asn (k)
      egso(k) = egsn(k)
    2 continue
c
c=======================================================================
c     extrapolate radii and masses to new time level
c=======================================================================
c
      jext = 1
c
    4 if (leibc .gt. 0) drdt(ngrs  ) =   0.0D0
      if (leibc .gt. 0) dmdt(ngrs  ) = - rn  (ngrs  )**mu * phil0
      if (leobc .gt. 0) drdt(ngre+1) =   0.0D0
      if (leobc .lt. 3) dmdt(ngre+1) = - rn  (ngre+1)**mu * phir0
      if (leobc .eq. 3) dmdt(ngre+1) =   dmdt(ngre+1)
      if (llibc .gt. 0) drdt(ngrs  ) =   un  (ngrs  )
      if (llibc .gt. 0) dmdt(ngrs  ) =   0.0D0
      if (llobc .gt. 0) drdt(ngre+1) =   un  (ngre+1)
      if (llobc .gt. 0) dmdt(ngre+1) =   0.0D0
c
      do 5 k = ngrs, ngre + 1
      rn  (k) = ro  (k) + drdt(k) * dtime
      xmen(k) = xmeo(k) - dmdt(k) * dtime
    5 continue
c
c-----------------------------------------------------------------------
c     check radius and mass for monotonicity
c-----------------------------------------------------------------------
c
      do 6 k = ngrs, ngre
      if ( rn(k+1) .le. rn(k) .or. xmen(k+1) .gt. xmen(k) ) then
c
c          nonmonotonic, cut timestep and try again
c
           write(itty, 100) jext, dtime
           write(iout, 100) jext, dtime
           write(iout, 101) k+1,rn(k+1),k,rn(k), k+1,xmen(k+1),k,xmen(k)
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
           jstep = jstepe
           call writout
           call clzfil
           stop 'peq'
      end if
    6 continue
  100 format(' grid extrapolation nonmonotonic. jext ='i3
     .       ' dtime =' 1pe10.2)
  101 format(  ' rn('i3') ='1pe22.14  ' rn('i3') ='e22.14
     .       ' xmen('i3') ='  e22.14' xmen('i3') ='e22.14)
  102 format(' grid extrapolation nonmonotonic after'i3' attempts')
c
c=======================================================================
c     boundary conditions
c=======================================================================
c-----------------------------------------------------------------------
c     velocity
c-----------------------------------------------------------------------
c
c     inner boundary; equations bc16 & bc66
c
      k = ngrs
c
      if (leibc .le. 2) un(k) = ul
      if (llibc .eq. 1) un(k) = uextl
c
c     outer boundary; equations bc32 & bc69
c
      k = ngre + 1
c
      if (leobc .le. 2) un(k) = ur
      if (llobc .eq. 1) un(k) = uextr
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
          tn  (k) =  tn  (k+1)
	  dn  (k) =  dn  (k+1)
	  rn  (k) = (rn  (k+1)**mup1 - xmup1 * delml / dn(k))**rxm1
	  xmen(k) =  xmen(k+1)               + delml
          un  (k) =  un  (k+1)
      end if
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
	  tn  (k  ) =  tn  (k-1)
	  tn  (k+1) =  tn  (k  )
	  dn  (k  ) =  dn  (k-1)
	  dn  (k+1) =  dn  (k  )
	  rn  (k+1) = (rn  (k  )**mup1 + xmup1 * delmr / dn(k) )**rxm1
	  xmen(k+1) =  xmen(k  )               - delmr
	  un  (k+1) =  un  (k  )
      end if
c
c=======================================================================
c     compute new radius and mass related quantities
c=======================================================================
c
      do 7 k = ngrs - 1, ngre + 2
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
      xmn   (k) = xmtot - xmen(k)
    7 continue
c
      do 8 k = ngrs - 1, ngre + 1
      dvoln(k) = (rmup1n(k+1) - rmup1n(k)) / xmup1
    8 continue
c
c=======================================================================
c     eos
c=======================================================================
c
      call eos
c
c=======================================================================
c     calculate time-centered quantities
c=======================================================================
c
      do 11 k = ngrs - 1, ngre + 2
      r    (k) = thet * rn  (k) + (1.D0 - thet) * ro  (k)
      u    (k) = thet * un  (k) + (1.D0 - thet) * uo  (k)
      xm   (k) = thet * xmn (k) + (1.D0 - thet) * xmo (k)
      xme  (k) = thet * xmen(k) + (1.D0 - thet) * xmeo(k)
      dmdt (k) = - ( xmen(k) - xmeo(k) ) / dtime
      drdt (k) =   ( rn  (k) - ro  (k) ) / dtime
      urel (k) =     u   (k) - drdt(k)
      rmu  (k) = r(k)**mu
      rmup1(k) = r(k)**mup1
      rmum1(k) = r(k)**mum1
      unom (k) = max(abs(un(k)),1.D0)
   11 continue
c
      do 12 k = ngrs - 1, ngre + 1
      d    (k) = thet * dn  (k) + (1.D0 - thet) * do  (k)
      t    (k) = thet * tn  (k) + (1.D0 - thet) * to  (k)
      pg   (k) = thet * pgn (k) + (1.D0 - thet) * pgo (k)
      eg   (k) = thet * egn (k) + (1.D0 - thet) * ego (k)
      as   (k) = thet * asn (k) + (1.D0 - thet) * aso (k)
      dvol (k) =  ( rmup1(k+1) - rmup1(k) ) / xmup1
   12 continue
c
c=======================================================================
c     zero out turbulence- and radiation-related quantities
c=======================================================================
c
      do 14 k = ngrs - 1, ngre + 2
      eto   (k) = 0.0D0
      et    (k) = 0.0D0
      etn   (k) = 0.0D0
      etso  (k) = 0.0D0
      etb   (k) = 0.0D0
      etsn  (k) = 0.0D0
      ero   (k) = 0.0D0
      er    (k) = 0.0D0
      ern   (k) = 0.0D0
      erso  (k) = 0.0D0
      erb   (k) = 0.0D0
      ersn  (k) = 0.0D0
      egto  (k) = 0.0D0
      egt   (k) = 0.0D0
      egtn  (k) = 0.0D0
      egtso (k) = 0.0D0
      egtb  (k) = 0.0D0
      egtsn (k) = 0.0D0
      egro  (k) = ego(k)
      egr   (k) = eg (k)
      egrn  (k) = egn(k)
      egrso (k) = 0.0D0
      egrb  (k) = 0.0D0
      egrsn (k) = 0.0D0
      egrto (k) = egro(k)
      egrt  (k) = egr (k)
      egrtn (k) = egrn(k)
      egrtso(k) = 0.0D0
      egrtb (k) = 0.0D0
      egrtsn(k) = 0.0D0
      fro   (k) = 0.0D0
      fr    (k) = 0.0D0
      frn   (k) = 0.0D0
      frso  (k) = 0.0D0
      frb   (k) = 0.0D0
      frsn  (k) = 0.0D0
      chifo (k) = 0.0D0
      chif  (k) = 0.0D0
      chifn (k) = 0.0D0
      avchio(k) = 0.0D0
      avchi (k) = 0.0D0
      avchin(k) = 0.0D0
      xkeo  (k) = 0.0D0
      xke   (k) = 0.0D0
      xken  (k) = 0.0D0
      xkpo  (k) = 0.0D0
      xkp   (k) = 0.0D0
      xkpn  (k) = 0.0D0
      xneo  (k) = 0.0D0
      xne   (k) = 0.0D0
      xnen  (k) = 0.0D0
      plfo  (k) = 0.0D0
      plf   (k) = 0.0D0
      plfn  (k) = 0.0D0
      fedd  (k) = 0.0D0
      frnom (k) = 0.0D0
      bri   (k) = unom(k)
   14 continue
c
      return
      end
