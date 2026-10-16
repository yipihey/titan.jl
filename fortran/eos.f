      subroutine eos
c
c***********************************************************************
c     equation of state
c.......................................................................
c     calling sequence: peqrh, peqh, peqr, start...,
c                       updaterh, updateh, updater > eos > intrp
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      dimension dlg(mgr), tlg(mgr), xx(mgr), yy(mgr), fl(mgr)
c
      k1 = ngrs - 1
      k2 = ngre + 2
      kk = k2 - k1 + 1
c
c=======================================================================
c     table lookups (MHD equation of state)
c     source: D. MIHALAS, W. DAPPEN, D. HUMMER, and B. MIHALAS, in
c     Ap. J., 331, 794, 1988; Ap. J., 331, 815, 1988;
c     Ap. J., 332, 261, 1988; Ap. J., 350, 300, 1990.
c=======================================================================
c
c+++++++++++++++++++++++++++ leos .eq. 1 +++++++++++++++++++++++++ begin
      if (leos .eq. 1) then
c
      do 11 k = k1, k2
      tlg(k) = log(tn(k)) 
      dlg(k) = log(dn(k))
   11 continue
c
c-----------------------------------------------------------------------
c     gas pressure
c-----------------------------------------------------------------------
c
      call intrp( kk, dlg(k1), tlg(k1), fpg, fpgx, fpgy, fpgxy, xx, yy,
     .                         pgn(k1), dlpgdldn(k1), dlpgdltn(k1), fl)
c
      do 20 k = k1, k2
      pgn(k) = pgn(k) + ( dlg(k) - xx(k) ) * dlpgdldn(k)
     .                + ( tlg(k) - yy(k) ) * dlpgdltn(k)
      pgn(k) = exp( pgn(k) )
   20 continue
c
          noff = int( ssum(kk, fl, 1) )
      if (noff. gt. 0) then
	    write(iout, 100) noff
	    write(iout,'(25f3.0)') fl
      end if
  100 format( ' went of edges of eos table' i4 ' times' )
c
c-----------------------------------------------------------------------
c     gas energy density
c-----------------------------------------------------------------------
c
      call intrp( kk, dlg(k1), tlg(k1), feg, fegx, fegy, fegxy, xx, yy,
     .                         egn(k1), dlegdldn(k1), dlegdltn(k1), fl)
c
      do 30 k = k1, k2
      egn(k) = egn(k) + ( dlg(k) - xx(k) ) * dlegdldn(k)
     .                + ( tlg(k) - yy(k) ) * dlegdltn(k)
      egn(k) = exp( egn(k) )
   30 continue
c
c-----------------------------------------------------------------------
c     electron number density
c-----------------------------------------------------------------------
c
      call intrp( kk, dlg(k1), tlg(k1), fpe, fpex, fpey, fpexy, xx, yy,
     .                        xnen(k1), dlpedldn(k1), dlpedltn(k1), fl)
c
c
      do 40 k = k1, k2
      xnen(k) = xnen(k) + ( dlg(k) - xx(k) ) * dlpedldn(k)
     .                  + ( tlg(k) - yy(k) ) * dlpedltn(k)
      xnen(k) = exp( xnen(k) ) / ( ck * tn(k) )
   40 continue
c
c-----------------------------------------------------------------------
c     sound speed
c-----------------------------------------------------------------------
c
      do 50 k = k1, k2
Cmhg  gamman(k) = dlpgdldn(k) * cpn(k) / cvn(k) [cvn = tn/egn*dlegdltn]
      gamman(k) = gam
      asn   (k) = sqrt( abs( gamman(k) * pgn(k) / dn(k) ) )
   50 continue
c
      return
c+++++++++++++++++++++++++++ leos .eq. 1 +++++++++++++++++++++++++++ end
      end if
c
c=======================================================================
c     perfect gas
c=======================================================================
c
c+++++++++++++++++++++++++++ leos .eq. 2 +++++++++++++++++++++++++ begin
      if (leos .eq. 2) then
c
                       do 60 k = k1, k2
                       pgn  (k) =   ck  / (cm0 * gmu) * dn(k) * tn(k)
                       egn  (k) =   ck  / (cm0 * gmu)         * tn(k)
     .                                  / (gam - 1.D0)
                       xnen (k) = 0.5D0 / (cm0 * gmu) * dn(k)
 		       delsn(k) =         (gam - 1.D0) / gam
 		       cpn  (k) =   ck  / (cm0 * gmu)  / delsn(k)
 		       betan(k) = 1.0D0 / ( tn(k) * dn(k) )
 		      gamman(k) =       gam
                       asn  (k) = sqrt( gam * pgn(k) / dn(k) )
c
                       dlpgdltn(k) = + 1.0D0
                       dlpgdldn(k) = + 1.0D0
                       dlegdltn(k) = + 1.0D0
                       dlegdldn(k) =   0.0D0
                       dlpedltn(k) = + 1.0D0
                       dlpedldn(k) = + 1.0D0
                       dldsdltn(k) =   0.0D0
                       dldsdldn(k) =   0.0D0
                       dlcpdltn(k) =   0.0D0
                       dlcpdldn(k) =   0.0D0
                       dlbedltn(k) = - 1.0D0
                       dlbedldn(k) = - 1.0D0
   60                  continue
                       return
      end if
c+++++++++++++++++++++++++++ leos .eq. 2 +++++++++++++++++++++++++++ end
c
c=======================================================================
c     eos fit according to R.F. STELLINGWERF (1984) ApJ, 262, 336 
c=======================================================================
c
c+++++++++++++++++++++++++++ leos .eq. 3 +++++++++++++++++++++++++ begin
      if (leos .eq. 3) then
		             epseos = 1.D-9
c
                 do 70 k = k1, k2
c
		 call eos3 ( epseos,tn(k),dn(k),xabun,yabun,zabun
     .                     , pe , pe v, pe t, pe vv, pe tt, pe vt
     .                     , pn , pn v, pn t, pn vv, pn tt, pn vt
     .                     , en , en v, en t, en vv, en tt, en vt )
c
		 egn  (k) =   en
		 pgn  (k) =   pn
                 xnen (k) = + pe   / ( ck * tn(k) )
 		 betan(k) = - pn t /  pn v
 		 cpn  (k) =   en t + betan(k) * pn t
 		 delsn(k) =   pn   * betan(k) / cpn(k)
c
                 dlegdltn(k) = + (en t /  egn(k)) * tn(k)
                 dlegdldn(k) = -  en v / (egn(k)  * dn(k))
                 dlpgdltn(k) = + (pn t /  pgn(k)) * tn(k)
                 dlpgdldn(k) = -  pn v / (pgn(k)  * dn(k))
                 dlpedltn(k) = + (pe t /  pe    ) * tn(k)
                 dlpedldn(k) = -  pe v / (pe      * dn(k))
c
                 dbe dtn = - (pn tt - pn t/pn v * pn vt) / pn v
                 dbe dvn = - (pn vt - pn t/pn v * pn vv) / pn v
                 dcp dtn = en tt +  pn t  * betan(k)
     .                   + tn(k) * (pn tt * betan(k) + dbedtn * pn t)
                 dcp dvn = en vt
     .                   + tn(k) * (pn vt * betan(k) + dbedvn * pn t)
                 dds dtn = (betan(k) * (pn t - pn * dcpdtn / cpn(k))
     .                                       + pn * dbedtn)/ cpn(k)
                 dds dvn = (betan(k) * (pn v - pn * dcpdvn / cpn(k))
     .                                       + pn * dbedvn)/ cpn(k)
c
                 dlbedltn(k) = + dbedtn *  tn(k) / betan(k)
                 dlbedldn(k) = - dbedvn / (dn(k) * betan(k))
                 dlcpdltn(k) = + dcpdtn *  tn(k) /   cpn(k)
                 dlcpdldn(k) = - dcpdvn / (dn(k) *   cpn(k))
                 dldsdltn(k) = + ddsdtn *  tn(k) / delsn(k)
                 dldsdldn(k) = - ddsdvn / (dn(k) * delsn(k))
c
		 gamman(k) = - (pn v/en t) * cpn(k) / (pgn(k) * dn(k))
                 asn   (k) =   sqrt( abs( gamman(k) *  pgn(k) / dn(k)))
Cmhg------- need asn derivatives in turscl or set asn^2=gam p/rho ------
   70            continue
                 return
      end if
c+++++++++++++++++++++++++++ leos .eq. 3 +++++++++++++++++++++++++++ end
c
c=======================================================================
c
      stop 'eos01'
      end
      subroutine  eos3       ( epseos, tin,din, x,y,z
     1                       , pe , pe v, pe t, pe vv, pe tt, pe vt
     2                       ,  p ,  p v,  p t,  p vv,  p tt,  p vt
     3                       ,  e ,  e v,  e t,  e vv,  e tt,  e vt )
c
c     equation of state: first and second order derivatives
c     arguments........: tgn (K), dgn (g/cm**3)
c     metals...........: Na,Al always ionized,treated as single element
c                        Mg,Si,Fe            included as single element
c                                                    all others ignored
c     source...........: R.F.STELLINGWERF ApJ 262,p.336-338
c
c     thermodynamic relations :
c
c     e(T,rho) = fe(T,rho) * T
c                fe = feg + fei  gas+ionization
c     p(T,rho) = fp(T,rho) * T * rho
c                fp = fpg        gas
c
c          x = hydrogen mass fraction
c          y = helium   mass fraction
c          z = metallic mass fraction
c
c=======================================================================
c
      include 'titan.imp'
      include 'titan.par'
c
      data  cgas           /8.31451D7/
      data  xm,xh,xhe,xh2  /7.6D0,13.595D0,24.581D0,54.403D0/
c
      t =         tin
      v = 1.0D0 / din
c
      zh   = x /  1.00797D0
      zhe  = y /  4.00260D0
      zme  = z / 17.84128D0
      zh2  = zhe *  2.0D0
      zm   = zme *  5.30325D-2
      zna  = zme *  1.78330D-3
c
      tk   = 8.6167D-5 * t
      dtk  = 1.0D0 / tk
      ts   = 3.3334622D-1/cgas * v * t * sqrt(t)
      ts t = 1.5D0 * ts / t
      ts v = 1.0D0 * ts / v
                     tm  = 6.614536000D-34
      if (t.gt.5.D2) tm  = exp(-xm  *dtk)
                     th  = 5.347896000D-35
      if (t.gt.1.D3) th  = exp(-xh  *dtk)
                     the = 1.665795163D-25
      if (t.gt.5.D3) the = exp(-xhe *dtk)
                     th2 = 3.802592017D-28
      if (t.gt.1.D4) th2 = exp(-xh2 *dtk)
c
      fim      = 1.4D0 * ts * tm
      fih      =         ts * th
      fihe     = 4.0D0 * ts * the
      fih2     =         ts * th2
c
      fim   t  = + fim   * ( 1.5D0 + d tk*xm  ) / t
      fih   t  = + fih   * ( 1.5D0 + d tk*xh  ) / t
      fihe  t  = + fihe  * ( 1.5D0 + d tk*xhe ) / t
      fih2  t  = + fih2  * ( 1.5D0 + d tk*xh2 ) / t
c
      fim   v  = + fim   / v
      fih   v  = + fih   / v
      fihe  v  = + fihe  / v
      fih2  v  = + fih2  / v
c
      fim   tt = + fim   * (.75D0 + d tk*xm  *( 1.D0+d tk*xm  ) ) / t**2
      fih   tt = + fih   * (.75D0 + d tk*xh  *( 1.D0+d tk*xh  ) ) / t**2
      fihe  tt = + fihe  * (.75D0 + d tk*xhe *( 1.D0+d tk*xhe ) ) / t**2
      fih2  tt = + fih2  * (.75D0 + d tk*xh2 *( 1.D0+d tk*xh2 ) ) / t**2
c
      fim   vt = + fim   t / v
      fih   vt = + fih   t / v
      fihe  vt = + fihe  t / v
      fih2  vt = + fih2  t / v
c
c --- converge on electronic molecular weight
c
          wmy i = .5D0* ( sqrt( fih*(fih+4.D0*x) ) - fih )
      if  ( x .eq. 1.D0)  go to 60
          wmy i = .5D0* ( sqrt( fim*(fim+4.D0*z) ) - fim )
      if  ( z .eq. 1.D0)  go to 60
c
          wmy e =     ( x+.5D0*y ) / ( 1.D0+.25D0*y/fihe )
      if  ( wmy e.lt. x )  
     .    wmy e = .5D0* ( sqrt( fih*(fih+4.D0*x) ) - fih )
            delta = abs ( ( wmy i-wmy e )/wmy e )
      if  ( delta.lt.0.95D0)  wmy e = wmy i
c
      do  1000  iter = 1,99
c
      xim     =  fim  / ( wmy e+fim )
      xih     =  fih  / ( wmy e+fih )
      xihe    =  fihe / ( wmy e+fihe*( 1.D0+fih2 /wmy e ) )
      xih2    =  fih2 /   wmy e*xihe
      xim   w = - xim  **2        /fim
      xih   w = - xih  **2        /fih
      xihe  w = - xihe **2 * (1.D0/fihe    +fih2 /wmy e**2 )
      xih2  w = + fih2 /  wmy e*xihe w-xih2 /wmy e
c
        wmy a = + zna + zm*xim   + zh*xih   + zhe*xihe   + zh2*xih2
      d wmy a =       + zm*xim w + zh*xih w + zhe*xihe w + zh2*xih2 w
        wmyle =   log ( wmy e )
        wmyla =   log ( wmy a )
      d wmyla =       d wmy a * wmy e / wmy a
c
        wmyla =  wmyl e -   wmyl a
      d wmyla =  1.0D0  - d wmyl a
c
        delta = - wmyl a / d wmyl a
        wmyli =   wmyl e + delta
        wmy i =   exp      ( wmyl i )
c
        wmylo =   wmyle
        wmy o =   wmy e
        wmy e =   wmy i
c
            prec =  abs ( (wmyli-wmylo)/wmylo ) - epseos
      if  ( prec .lt. 0.D0)  go to 60
1000  continue
60    continue
c
c --- eos and derivatives
c
      xim   =  fim  / ( wmy e+fim )
      xih   =  fih  / ( wmy e+fih )
      xihe  =  fihe / ( wmy e+fihe*( 1.D0+fih2 /wmy e ) )
      xih2  =  fih2 /   wmy e*xihe
c
      xim   w = - xim  **2        /fim
      xih   w = - xih  **2        /fih
      xihe  w = - xihe **2 * (1.D0/fihe  +fih2 /wmy e**2 )
      xih2  w = + fih2 / wmy e*xihe w-xih2 /wmy e
c
      xim   m   = + ( xim  / fim  )**2 * wmy e
      xih   h   = + ( xih  / fih  )**2 * wmy e
      xihe  he  = + ( xihe / fihe )**2 * wmy e
      xihe  h2  = -   xihe         **2 / wmy e
      xih2  he  = +   xihe he  * fih2  / wmy e
      xih2  h2  = +   xihe h2  * fih2  / wmy e + xih2  / fih2
c
      xim   m  m   = - xim  m  *        xim  / fim
      xih   h  h   = - xih  h  *        xih  / fih
      xihe  he he  = + xihe he * ( 2.D0*xihe / fihe * wmy e - 1.D0)
      xihe  he h2  = + xihe he *   2.D0*xihe        / wmy e
      xihe  h2 h2  = - xihe h2 *   2.D0*xihe        / wmy e
      xih2  he he  = + xihe he he  * fih2                   / wmy e
      xih2  he h2  = ( xihe he h2  * fih2  +      xihe he ) / wmy e
      xih2  h2 h2  = ( xihe h2 h2  * fih2  + 2.D0*xihe h2 ) / wmy e
c
      xim   ww = - 2.D0* xim  w*xim  /fim
      xih   ww = - 2.D0* xih  w*xih  /fih
      xihe  ww = - 2.D0* xihe w*xihe /fihe
     1           + 2.D0*(xihe w*xihe +xihe h2 )*fih2 /wmy e**2
      xih2  ww = - 2.D0*(xihe w*fih2 -xih2    )      /wmy e**2
     1           +       xihe ww               *fih2 /wmy e
c
      xim   wm   = ( fim-wmy e )/( fim+wmy e)**3
      xih   wh   = ( fih-wmy e )/( fih+wmy e)**3
      xihe  whe  = xihe * ( xihe      +2.D0*xihe w*wmy e ) / fihe**2
      xihe  wh2  = xihe * ( xihe/wmy e-2.D0*xihe w       ) / wmy e
      xih2  whe  =   ( xihe whe /wmy e         -xihe he  ) * fih2
      xih2  wh2  =   ( xihe wh2 *fih2 + xihe w -xih2 h2  ) / wmy e
c
      xw m   t =  xim  m   * fim  t
      xw h   t =  xih  h   * fih  t
      xw he  t =  xihe he  * fihe t + xihe h2 * fih2  t
      xw h2  t =  xih2 h2  * fih2 t + xih2 he * fihe  t
c
      xw m   v =  xim  m   * fim  v
      xw h   v =  xih  h   * fih  v
      xw he  v =  xihe he  * fihe v + xihe h2 * fih2  v
      xw h2  v =  xih2 h2  * fih2 v + xih2 he * fihe  v
c
      xw m  tt =  xim  m  m  * fim  tt
      xw h  tt =  xih  h  h  * fih  tt
      xw he tt =  xihe he he * fihe tt + xihe h2 h2 * fih2 tt
     1                 + 2.D0* fihe t  * xihe he h2 * fih2 t
      xw h2 tt =  xih2 he he * fihe tt + xih2 h2 h2 * fih2 tt
     1                 + 2.D0* fihe t  * xih2 he h2 * fih2 t
c
      xw m  vt =  xim  m  m  * fim  vt
      xw h  vt =  xih  h  h  * fih  vt
      xw he vt =  xihe he he * fihe vt + xihe h2 h2 * fih2 vt
     1                       + fihe  t * xihe he h2 * fih2 v
     2                       + fihe v  * xihe he h2 * fih2  t
      xw h2 vt =  xih2 he he * fihe vt + xih2 h2 h2 * fih2 vt
     1                       + fihe  t * xih2 he h2 * fih2 v
     2                       + fihe v  * xih2 he h2 * fih2  t
c
      xw m  wt =  xim  w m   * fim  t
      xw h  wt =  xih  w h   * fih  t
      xw he wt =  xihe w he  * fihe t  + xihe w h2 * fih2  t
      xw h2 wt =  xih2 w he  * fihe t  + xih2 w h2 * fih2  t
c
      xw m  wv =  xim  w m   * fim  v
      xw h  wv =  xih  w h   * fih  v
      xw he wv =  xihe w he  * fihe v  + xihe w h2 * fih2  v
      xw h2 wv =  xih2 w he  * fihe v  + xih2 w h2 * fih2  v
c
c --- electron molecular weight and derivatives
c
      wmy e  = zna + zm*xim    + zh*xih    + zhe*xihe    + zh2 *xih2
      wmy w  = 1.D0- zm*xim w  - zh*xih w  - zhe*xihe w  - zh2 *xih2 w
      wmy ww = 1.D0- zm*xim ww - zh*xih ww - zhe*xihe ww - zh2 *xih2 ww
c
      wmy    =  x/1.00797D0 + y/4.0026D0 + z/17.84128D0 + wmy e
c
      wmy t  =  (zm*xwm t + zh*xwh t + zhe*xwhe t + zh2*xwh2  t )/wmy w
      wmy v  =  (zm*xwm v + zh*xwh v + zhe*xwhe v + zh2*xwh2  v )/wmy w
      wmy tt = ((zm*xwmtt + zh*xwhtt + zhe*xwhett + zh2*xwh2 tt )
     1         +(zm*xwmwt + zh*xwhwt + zhe*xwhewt + zh2*xwh2 wt )*wmy t
     2                                                    *2.D0 )/wmyww
      wmy vt = ((zm*xwmvt + zh*xwhvt + zhe*xwhevt + zh2*xwh2 vt )
     1         +(zm*xwmwt + zh*xwhwt + zhe*xwhewt + zh2*xwh2 wt )*wmy v
     2         +(zm*xwmwv + zh*xwhwv + zhe*xwhewv + zh2*xwh2 wv )*wmy t
     3                                                          )/wmyww
      wmy vv =   0.0D0
c
      xim   t  = xim  m  * fim  t + xim  w * wmy t
      xih   t  = xih  h  * fih  t + xih  w * wmy t
      xihe  t  = xihe he * fihe t + xihe w * wmy t + xihe h2 * fih2 t
      xih2  t  = xih2 he * fihe t + xih2 w * wmy t + xih2 h2 * fih2 t
c
      xim   v  = xim  m  * fim  v + xim  w * wmy v
      xih   v  = xih  h  * fih  v + xih  w * wmy v
      xihe  v  = xihe he * fihe v + xihe w * wmy v + xihe h2 * fih2 v
      xih2  v  = xih2 he * fihe v + xih2 w * wmy v + xih2 h2 * fih2 v
c
      xim   tt = xim  m m *fim  tt+xim  w w  *wmy tt
     1              + 2.D0*fim  t *xim  w m  *wmy t
      xih   tt = xih  h h *fih  tt+xih  w w  *wmy tt
     1              + 2.D0*fih  t *xih  w h  *wmy t
      xihe  tt = xihe hehe*fihe tt+xihe w w  *wmy tt+xihe h2 h2 *fih2 tt
     1              + 2.D0*fihe t *xihe w he *wmy t
     2              + 2.D0*fih2 t *xihe w h2 *wmy t
     3              + 2.D0*fihe t *xihe heh2 *fih2  t
      xih2  tt = xih2 hehe*fih2 tt+xih2 w w  *wmy tt+xih2 h2 h2 *fih2 tt
     1              + 2.D0*fihe t *xih2 w he *wmy t
     2              + 2.D0*fih2 t *xih2 w h2 *wmy t
     3              + 2.D0*fihe t *xih2 heh2 *fih2  t
c
      xim   tt = xim  m  m         *fim  tt + xim  w  w          *wmy tt
     1         + xim  w  m  *wmy  t*fim   t * 2.D0
      xih   tt = xih  h  h         *fih  tt + xih  w  w          *wmy tt
     1         + xih  w  h  *wmy  t*fih   t * 2.D0
      xihe  tt = xihe he he        *fihe tt + xihe w  w          *wmy tt
     1         + xihe he h2        *fih2 tt + 2.D0 *
     2         ( xihe w  he *wmy  t*fihe  t + xihe w  h2 *wmy  t*fih2  t
     3         + xihe he h2 *fihe t*fih2  t )
      xih2  tt = xih2 he he        *fihe tt + xih2 w  w          *wmy tt
     1         + xih2 he h2        *fih2 tt + 2.D0 *
     2         ( xih2 w  he *wmy  t*fihe  t + xih2 w  h2 *wmy  t*fih2  t
     3         + xih2 he h2 *fihe t*fih2  t )
c
      xim   vt = xim  m  m         *fim  vt + xim  w  w          *wmy vt
     1         + xim  w  m  *wmy  v*fim   t + xim  w  m  *fim   v*wmy  t
      xih   vt = xih  h  h         *fih  vt + xih  w  w          *wmy vt
     1         + xih  w  h  *wmy  v*fih   t + xih  w  h  *fih   v*wmy  t
      xihe  vt = xihe he he        *fihe vt + xihe w  w          *wmy vt
     1         + xihe w  he *wmy  v*fihe  t + xihe w  he *fihe  v*wmy  t
     2         + xihe he h2        *fih2 vt
     3         + xihe w  h2 *wmy  v*fih2  t + xihe w  h2 *fih2  v*wmy  t
     4         + xihe he h2 *fihe v*fih2  t + xihe he h2 *fih2  v*fihe t
      xih2  vt = xih2 he he        *fihe vt + xih2 w  w          *wmy vt
     1         + xih2 w  he *wmy  v*fihe  t + xih2 w  he *fihe  v*wmy  t
     2         + xih2 he h2        *fih2 vt
     3         + xih2 w  h2 *wmy  v*fih2  t + xih2 w  h2 *fih2  v*wmy  t
     4         + xih2 he h2 *fihe v*fih2  t + xih2 he h2 *fih2  v*fihe t
c
c --- pressure and derivatives
c
      pe    =                                     wmy e  * cgas* t/v
      pe t  = pe * ( wmy t /wmy e                        + 1.D0/ t   )
      pe v  = pe * ( wmy v /wmy e                        - 1.D0/   v )
      pe tt = pe * ( wmy tt/wmy e + 2.D0*wmy t/(t*wmy e)             )
      pe vv = pe * ( wmy vv/wmy e - 2.D0*wmy v/(v*wmy e) + 2.D0/(v*v))
      pe vt = pe * ( wmy vt/wmy e -      wmy v/(t*wmy e)
     1                            -      wmy t/(v*wmy e) + 1.D0/(t*v))
c
      p     =                                     wmy    * cgas* t/v
      p  t  = p  * ( wmy t /wmy                          + 1.D0/ t   )
      p  v  = p  * ( wmy v /wmy                          - 1.D0/   v )
      p  tt = p  * ( wmy tt/wmy   + 2.D0*wmy t/(t*wmy  )             )
      p  vv = p  * ( wmy vv/wmy   - 2.D0*wmy v/(v*wmy  ) + 2.D0/(v*v))
      p  vt = p  * ( wmy vt/wmy   -      wmy v/(t*wmy  )
     1                            -      wmy t/(v*wmy  ) + 1.D0/(t*v))
c
c --- energy and derivatives
c
      eg    =      1.5D0 *  wmy * cgas* t
      eg t  = eg * ( wmy t /wmy + 1.D0/ t )
      eg v  = eg *   wmy v /wmy
      eg tt = eg *   wmy tt/wmy
      eg vt = eg *   wmy vt/wmy
c
              rbk = cgas / 8.6167D-5 
c
      ei    = rbk * ( zh * xh *xih    + zhe *  xhe      * xihe
     1              + zm * xm *xim    + zhe * (xhe+xh2) * xih2
     2              + zm * 0.185735D0 + zh  * 0.7558696D0        )
      ei t  = rbk * ( zh * xh *xih t  + zhe *  xhe      * xihe t
     1              + zm * xm *xim t  + zhe * (xhe+xh2) * xih2 t  )
      ei v  = rbk * ( zh * xh *xih v  + zhe *  xhe      * xihe v
     1              + zm * xm *xim v  + zhe * (xhe+xh2) * xih2 v  )
      ei tt = rbk * ( zh * xh *xih tt + zhe *  xhe      * xihe tt
     1              + zm * xm *xim tt + zhe * (xhe+xh2) * xih2 tt )
      ei vt = rbk * ( zh * xh *xih vt + zhe *  xhe      * xihe vt
     1              + zm * xm *xim vt + zhe * (xhe+xh2) * xih2 vt )
      ei vv = 0.0D0
c
      e    =  eg    + ei
      e t  =  eg t  + ei t
      e v  =  eg v  + ei v
      e tt =  eg tt + ei tt
      e vt =  eg vt + ei vt
      e vv =  0.0D0
c
c-----------------------------------------------------------------------
c
      return
      end
