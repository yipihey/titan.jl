      subroutine totnrg
c
c***********************************************************************
c     check total energy conservation te2
c.......................................................................
c     called by: step
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     calculate variables at actual midpoint of timestep
c=======================================================================
c
      do 1 k = ngrs, ngre + 1
      r    (k) =         (rn  (k) + ro  (k)) / 2.0D0
      xm   (k) = xmtot - (xmen(k) + xmeo(k)) / 2.0D0
      u    (k) =         (un  (k) + uo  (k)) / 2.0D0
      fr   (k) =         (frn (k) + fro (k)) / 2.0D0
      fc   (k) = 0.0D0 
      ft   (k) = 0.0D0 
      rmu  (k) = r(k)**mu
      rmup1(k) = r(k)**mup1
      rmum1(k) = r(k)**mum1
    1 continue
c
      do 2 k = ngrs, ngre
      d   (k) = 0.5D0 * (dn   (k) + do   (k))
      eg  (k) = 0.5D0 * (egn  (k) + ego  (k))
      pg  (k) = 0.5D0 * (pgn  (k) + pgo  (k))
      er  (k) = 0.5D0 * (ern  (k) + ero  (k))
      et  (k) = 0.5D0 * (etn  (k) + eto  (k))
      egrt(k) = 0.5D0 * (egrtn(k) + egrto(k))
      dvol(k) = (rmup1(k+1) - rmup1(k)) / xmup1
    2 continue
c
c-----------------------------------------------------------------------
c     internal, radiation, turbulent, kinetic, potential energy;
c     equation te3
c-----------------------------------------------------------------------
c
      do 3 k = ngrs, ngre
c
      tescr(k) = dn(k) * dvoln(k) * 
     .                            (  egrtn(k)
     .
     .         + 0.25D0 *            (  un(k)**2    +  un(k+1)**2      )
     .
     .+ (0.5D0 - 0.25D0 * xmu) * g * (  rn(k)       +  rn(k+1)         )
     .
     .         - xmu * cpi * cgrav * ( xmn(k)/rn(k) + xmn(k+1)/rn(k+1) )
     .                             )
     .                             - tescr0(k)    !that was set in start
c
    3 continue
c
      tee = ssum( ngre - ngrs + 1, tescr(ngrs), 1 )
c
c-----------------------------------------------------------------------
c     surface losses; equation te4
c-----------------------------------------------------------------------
c
      urel_ne = u(ngre+1) - ( rn(ngre+1)-ro(ngre+1) )/dtime
c
                                         xmdot = phir0
      if (leobc .eq. 3)                  xmdot = d(ngre) * urel_ne
c
      tes = teso + dtime * rmu(ngre+1) * xmdot * 
     .                              ( egrt(ngre  )
     .
     . +      0.5D0               * ( u   (ngre+1)**2 - u0(ngre+1)**2)
     .
     . + g * (0.5D0 - 0.25D0*xmu) * ( r   (ngre+1)    - r0(ngre+1)   )
     .
     . - cgrav * cpi * 2.0D0*xmu  * ( xm  (ngre+1)/r (ngre+1) 
     .                              - xm0 (ngre+1)/r0(ngre+1)        )
     .                              )
c
      tes = tes  - dtime * rmu(ngrs  ) * phil0 * 
     .                              ( egrt(ngrs  ) 
     .
     . +      0.5D0               * ( u   (ngrs  )**2 - u0(ngrs  )**2)
     .
     . + g * (0.5D0 - 0.25D0*xmu) * ( r   (ngrs  )    - r0(ngrs  )   )
     .
     . - cgrav * cpi * 2.0D0*xmu  * ( xm  (ngrs  )/r (ngrs  )
     .                              - xm0 (ngrs  )/r0(ngrs  )        )
     .                              )
c
c-----------------------------------------------------------------------
c     work; equation te5
c-----------------------------------------------------------------------
c
      tew = tewo + dtime *
     .
     .   ( rmu(ngre+1) * u(ngre+1) * (pg(ngre) + fedd(ngre) * er(ngre)
     .                       + (2.0D0 / 3.0D0) *    d(ngre) * et(ngre))
     .
     .   - rmu(ngrs  ) * u(ngrs  ) * (pg(ngrs) + fedd(ngrs) * er(ngrs)
     .                       + (2.0D0 / 3.0D0) *    d(ngrs) * et(ngrs))
     .   )
c
c-----------------------------------------------------------------------
c     luminosity; equation te6
c-----------------------------------------------------------------------
c
      tel = telo + dtime *
     .
     .    ( + rmu(ngre+1) * ( fr(ngre+1) + fc(ngre+1) + ft(ngre+1))
     .      - rmu(ngrs  ) * ( fr(ngrs  ) + fc(ngrs  ) + ft(ngrs  ))
     .    )
c
c-----------------------------------------------------------------------
c     dissipation; equation te11
c-----------------------------------------------------------------------
c
      teq = teqo + dtime *
     .
     .    ( + rmu(ngre+1) * u(ngre+1) * ( qf(ngre+1) + qv(ngre+1))
     .      - rmu(ngrs  ) * u(ngrs  ) * ( qf(ngrs  ) + qv(ngrs  ))
     .    )
c
c=======================================================================
c     sum the total, write it out; equation te2
c=======================================================================
c
      tetot = tee + tes + tew + tel - teq
c
      if (ltalk .eq. 1 .and. (mod(jstep,10) .eq. 0))
     .  write(itty,'(/" total energy check:"14x"tetot ="1pe16.8)') tetot
c
      if (ltalk .gt. 1) then
          write(itty,'(/" total energy check")')
          write(itty, 6) tetot, tee, tes, tew, tel, teq
      end if
      if (ltalk .gt. 0) then
          write(iout,'(/" total energy check")')
          write(iout, 6) tetot, tee, tes, tew, tel, teq
      endif
c
      if (inrg .gt. 0) then
          if (jstep .eq. jsteps) write(inrg, 7) header 
          if (jstep .eq. jsteps) write(inrg, 8) 
	   write(inrg, 9) jstep, timen, dtime, tetot,tee,tes,tew,tel,teq
      end if
c
c=======================================================================
c
    6 format(' total ='1pe16.8' tee ='e16.8' tes ='e16.8/
     .       '   tew ='1pe16.8' tel ='e16.8' teq ='e16.8)
    7 format(/a80)
    8 format(/'jstep       timen       dtime           total'         ,
     . '         tee         tes         tew         tel         teq' )
    9 format(i5,1p2e12.4,e16.8,5e12.4)
c
      return
      end
