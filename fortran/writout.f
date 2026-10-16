      subroutine writout
c
c***********************************************************************
c     output physical variables
c.......................................................................
c     called by: titan, start{01,02,03,...,10}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      real*8       odepth(mgr),  flux(mgr), trad(mgr), 
     .              shock(mgr), radia(mgr), resu(mgr)
      character*2      ch(mgr)
c
           jmod = mod( jstep - jsteps + 1, nout )
      if ( jmod .eq. 0 .or. jstep .eq. jstepe 
     .                 .or. timen .ge. tmax   ) then
c
c=======================================================================
c     auxiliary variables
c=======================================================================
c
           fedd  (ngre+1) = geddr
           fedd  (ngrs-1) = geddl
	   odepth(ngre+2) = 0.0D0
	   odepth(ngre+1) = 0.0D0
           xnu   (ngre+2) = 0.0D0
           xnu   (ngre+1) = timen
	   resu  (ngre+2) = 0.0D0
           resu  (ngre+1) = timen
           resu  (ngrs-1) = 0.0D0
           qx    (ngre+2) = 0.0D0
           qx    (ngre+1) = dtime
           ch    (ngre+2) = '  '
           ch    (ngre+1) = '  '
           ch    (ngrs-1) = '  '
           trad  (ngre+1) =    sqrt(sqrt(abs(ern(ngre+1)/car)))
     .                    *  ern(ngre+1)/abs(ern(ngre+1)+1.D-30)
           trad  (ngrs-1) =    sqrt(sqrt(abs(ern(ngrs-1)/car)))
     .                    *  ern(ngrs-1)/abs(ern(ngrs-1)+1.D-30)
	   flux  (ngre+1) =    frn(ngre+1) 
     .                      +   un(ngre+1) * ern(ngre+1) * (1.D0+geddl)
	   flux  (ngrs-1) =    frn(ngrs-1) 
     .                      +   un(ngrs-1) * ern(ngrs-1) * (1.D0+geddl)
           shock (ngre+1) =              abs( un(ngre+1)/asn(ngre+1))
           shock (ngrs-1) =              abs( un(ngrs-1)/asn(ngrs-1))
c
	   do 20 k = ngre, ngrs, -1
           odepth(k) = odepth(k+1) + dn(k) * chifn(k) * (rn(k+1)-rn(k))
           flux  (k) = frn(k) +  ern(k) * un(k) * (1.D0 + fedd(k))
           radia (k) = frn(k) * rmun(k)
           trad  (k) = sqrt(sqrt(abs(ern(k)/car)))
     .                 *  ern(k)/abs(ern(k)+1.D-30)
           shock (k) =           abs( un(k)/asn(k))
           resu  (k) = ( rn(k+1)+rn(k) )/( rn(k+1)-rn(k) )/2.D0
           ch    (k) = '  '
   20      continue
c
          imach   = isamax(ngre-ngrs, shock(ngrs),1) + ngrs - 1
      if (imach  .ge. ngrs .and. imach  .le. ngre) ch (imach ) = ' M'
          ishock  = isamax(ngre-ngrs, qf   (ngrs),1) + ngrs - 1
      if (ishock .ge. ngrs .and. ishock .le. ngre) ch (ishock) = ' S'
          iradia  = isamax(ngre-ngrs, radia(ngrs),1) + ngrs - 1
      if (iradia .ge. ngrs .and. iradia .le. ngre) ch (iradia) = 'R '
      if (iradia .eq. imach  )                     ch (iradia) = 'RM'
      if (iradia .eq. ishock )                     ch (iradia) = 'RS'
c
c=======================================================================
c     problem  1: first FAVOURITE PROBLEM
c=======================================================================
c
      if (lprob .eq.  1) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
      end if
c
c=======================================================================
c     problem  2: second FAVOURITE PROBLEM
c=======================================================================
c
      if (lprob .eq.  2) then
	   write(iout,  1) header 
	   write(iout,  3) jstep, timen, dtime
      end if
c
c=======================================================================
c     problem  3: SOD SHOCK TUBE
c=======================================================================
c
      if (lprob .eq.  3) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout,  4)
           write(iout,  5) (k,  rn(k), xnu(k), dn(k), qx(k), egn(k)
     .                       , pgn(k),  un(k),     4.D0*cpi*xmen(k)
     .                     ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c     problem  4: WOODWARD-COLLELA BLAST WAVE
c=======================================================================
c
      if (lprob .eq.  4) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout,  4)
           write(iout,  5) (k,  rn(k), xnu(k), dn(k), qx(k), egn(k)
     .                       , pgn(k),  un(k),     4.D0*cpi*xmen(k)
     .                     ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c     problem  5: SEDOV-TAYLOR BLAST WAVE
c=======================================================================
c
      if (lprob .eq.  5) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout,  6)
           write(iout,  7) (k,  rn(k), xnu(k), dn(k), qx(k), egn(k)
     .                       , pgn(k),  un(k), 4.D0*cpi*xmn(k)/cmsol
     .                     ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c     problem  6: HEATING TOWARDS RADIATIVE EQUILIBRIUM
c=======================================================================
c
      if (lprob .eq.  6) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout,  8)
           write(iout,  7)  (k,                  rn(k)/crsol
     .                        ,       4.D0*cpi*xmen(k)/cmsol
     .                        ,  ern(k), frn (k)
     .                        ,   tn(k), fedd(k)
     .                        , 4.D0*cpi*rmun(k)*frn(k)
     .                        , odepth(k) 
     .                      ,k=ngre+1,ngrs-1,-1)
           write(iout,  9)   timen, 4.0D0*cpi*rmun(ngre)*frn(ngre), 
     .                       geddl,geddr
      end if
c
c=======================================================================
c     problem  7: COOLING FROM RADIATIVE EQUILIBRIUM
c=======================================================================
c
      if (lprob .eq.  7) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout,  8)
           write(iout,  7)  (k,                  rn(k)/crsol
     .                        ,       4.D0*cpi*xmen(k)/cmsol
     .                        ,  ern(k), frn (k)
     .                        ,   tn(k), fedd(k)
     .                        , 4.D0*cpi*rmun(k)*frn(k)
     .                        , odepth(k) 
     .                      ,k=ngre+1,ngrs-1,-1)
           write(iout,  9)   timen, 4.0D0*cpi*rmun(ngre)*frn(ngre), 
     .                       geddl,geddr
      end if
c
c=======================================================================
c     problem  8: SUBCRITICAL SHOCK
c=======================================================================
c
      if (lprob .eq.  8) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout, 17)
           write(iout, 18)  (ch(k), k,       rn(k)
     .                        ,   4.D0*cpi*xmen(k)
     .                        ,    dn(k),odepth(k)
     .                        ,    tn(k),rmun(k)*4.D0*cpi*frn(k)/clsol
     .                        ,  trad(k),  un(k)
     .                        ,   pgn(k),fedd(k)
     .                        ,   xnu(k)
     .                      ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c     problem  9: SUPERCRITICAL SHOCK
c=======================================================================
c
      if (lprob .eq.  9) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout, 12)
           write(iout,  7)  (k,   rn(k), trad(k)
     .                        ,   dn(k), tn(k), odepth(k)
     .                        , fedd(k),rmun(k)*4.D0*cpi*frn(k)
     .                        ,   un(k)
     .                      ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c     problem 10: RADIATIVE BLAST WAVE
c=======================================================================
c
      if (lprob .eq. 10) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout, 17)
           write(iout, 18)  (ch(k), k,       rn(k)
     .                        ,   4.D0*cpi*xmen(k)
     .                        ,    dn(k),odepth(k)
     .                        ,    tn(k),rmun(k)*4.D0*cpi*frn(k)/clsol
     .                        ,  trad(k),  un(k)
     .                        ,   pgn(k),fedd(k)
     .                        ,   xnu(k)
     .                      ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c     problem 11: DeGREGORIA WHITE DWARF ACCRETION
c=======================================================================
c
      if (lprob .eq. 11) then
           write(iout,  1) header 
           write(iout,  3) jstep, timen, dtime
c
           write(iout, 17)
           write(iout, 18)  (ch(k), k,       rn(k)
     .                        ,   4.D0*cpi*xmen(k)
     .                        ,    dn(k),odepth(k)
     .                        ,    tn(k),rmun(k)*4.D0*cpi*frn(k)/clsol
     .                        ,  trad(k),  un(k)
     .                        ,   pgn(k),fedd(k)
     .                        ,   xnu(k)
     .                      ,k=ngre+1,ngrs-1,-1)
      end if
c
c=======================================================================
c
      return
      end if
c
    1 format(/a80)
    2 format(/' mass ='1pe10.2' lum ='e10.2' rad ='e10.2' teff ='0pf7.0)
    3 format(/' jstep ='i6' timen ='1pe12.4' dtime ='e12.4)
    4 format(/'  k    radius       resolution  density     art_stress',
     .                            ' gas energy gas pressure  velocity')
    5 format(i5, 1pe13.5, 7e12.4)
    6 format(/'  k    radius      resolution  density     art_stress',
     .              ' gas energy gas pressure  velocity    mass'     )
    7 format(i5, 1p8e12.4)
    8 format(/'  k    radius      mass        E_rad          F_rad'  ,
     .           '    T_gas       f_edd       luminosity  tau'       )
    9 format(/'  time ='1pe11.4' luminosity ='e11.4
     .         '  geddl ='e11.4    '  geddr ='e11.4                  )
   10 format(/'  k    radius      ext. M      d           gas T  '   ,
     .          '     opt.depth   fedd        luminosity  rad T  '   )
   11 format(/'  k    radius      mass        d           T      '   ,
     .          '     E           F           u           P      '   )
   12 format(/'  k    radius      M_ext       d           T      '   ,
     .          '     P_g           P_vis     luminosity  velocity'  )
   13 format(/'  k    radius      mass        d           T_gas  '   ,
     .          '     T_rad         f_edd     luminosity  velocity'  )
   14 format(/'  k    radius      ext. M      d           gas T  '   ,
     .          '     opt.depth   fedd        flux        velocity'  )
   15 format(/'  k    radius      ext. M      d           gas T  '   ,
     .          '     tau         luminosity  rad T       u      '   ,
     .          '     P_gas       P_vis'                             )
   16 format(a2, i3, 1p10e12.4)
   17 format(/'    k  radius      ext. M      d           tau    '   ,
     .          '     gas T       luminosity  rad T       u      '   ,
     .          '     gas P       f           resolution'            )
   18 format(a2, i3, 1p2e12.5, 9e12.4)
c
      end
