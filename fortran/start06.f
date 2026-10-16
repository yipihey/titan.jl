      subroutine start06
c
c***********************************************************************
c     setup for HEATING of grey sphere to RADIATIVE EQUILIBRIUM
c.......................................................................
c     calling sequence: start > start06 > pardoc, eddfac, eos, opac,
c                                         gridinit, writout
c***********************************************************************
c     for equilibrium diffusion  comparison results exist in the litera-
c     ture for times  t =  2.1e5,  1.1e6,  2.4e6,  4.9e6,  8.5e6, 1.1e7,
c     1.9e7, 1.8e10 seconds.
c
c     for full radiation transport  comparison results exist in the lit-
c     erature for times  t =  1.6e5,  1.9e6, 5.6e6, 1.2e7, 2.0e7, 3.9E10
c     seconds.
c***********************************************************************
c     notes:
c
c     Information about correct settings for the logical switches is im-
c     bedded in comments preceding each group of switches. 
c
c     The switches  "lcray"  and  "lband"  MUST be set by the  user.  If
c     an  adaptive grid  is used,  "ladx"  and  "lady(m)"  MUST  also be
c     set. These switches are given default values which are illegal, so
c     failure  to set them will cause the  code to abort.  Likewise, re-
c     quired  input parameters which are not  preassigned default values
c     are set to indefinites ( "xindef" ) so that the code will abort if
c     they are not set by the user. In addition  "ncell" , the number of
c     cells in the domain MUST be set.
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      character*80 yheader
c
c=======================================================================
c     general setup
c=======================================================================
c-----------------------------------------------------------------------
c     read in header and parameters
c-----------------------------------------------------------------------
c     xmass = total mass (solar masses)
c     xlum  = imposed luminosity (erg/s)
c     rho0  = density (gm/cm**3)
c     chif0 = opacity (cm**-1)
c-----------------------------------------------------------------------
c
      read (iin,'(a80)') header
      read (iin,'(2e10.4,e11.4,2e10.3)') xmass, xlum, rho0, chif0, ratio
c
      xmass = xmass * cmsol
      rmin  = 1.0D6
c
c-----------------------------------------------------------------------
c     geometry: lgeom = 2  => spherical geometry
c-----------------------------------------------------------------------
c
      lgeom = 2
       mu   = lgeom
       mup1 = mu + 1
       mum1 = mu - 1
      xmu   = float(mu  )
      xmup1 = float(mup1)
      xmum1 = float(mum1)
       rxm1 = 1.0D0 / xmup1
c
c-----------------------------------------------------------------------
c     gravity:
c
c           planar    geometry --
c                                 off: set g = 0.0D0
c                                 on : set g = desired value
c           spherical geometry -- 
c                                 off: set cgrav = 0.0D0
c                                 on : no action needed
c-----------------------------------------------------------------------
c
       g    = 0.0D0              ! that's ok too
      cgrav = 0.0D0              ! geometry is spherical; gravity is off
c     xmext = 0.0D0              ! set in exterior mass definition
c
c-----------------------------------------------------------------------
c     radiation transport: lrad  = 1  => radiation
c
c                          ltran = 0  => full transport solution and
c                                        fedd updated at each iteration
c                          ltran = 1  => full transport solution
c                          ltran = 2  => nonequilibrium diffusion
c                          ltran = 3  =>    equilibrium diffusion
c                                lam   = 0  => no flux limiting 
c                                lam   = 1  =>    flux limiting
c
c     ncor  = number of rays intersecting core
c-----------------------------------------------------------------------
c
      lrad  = 1
      ltran = 1
      lam   = 0
      ncor = 20
c
c-----------------------------------------------------------------------
c     radiation boundary conditions:
c     lribc = 3 => imposed luminosity at  inner boundary
c     lrobc = 1 => optically transmitting outer boundary
c-----------------------------------------------------------------------
c
      lribc = 3
      lrobc = 1
c
      xipl = 0.0D0
      ximl = 0.0D0
      xipr = 0.0D0
      ximr = 0.0D0
c
c-----------------------------------------------------------------------
c     turbulent transport: ltur  = 0  => no turbulence
c-----------------------------------------------------------------------
c
      ltur = 0
c
c-----------------------------------------------------------------------
c     declare equations to be solved (radiation only in static medium)
c     ir > 0 => grid         eqn;  it > 0 => gas energy eqn 
c     if > 0 => rad momentum eqn;  ie > 0 => rad energy eqn 
c
c     neqn  = number of equations at each grid-point
c     ngrs  = index of first (left /inner-most) cell
c     ngre  = index of last  (right/outer-most) cell
c     ncell = number of cells in the domain
c-----------------------------------------------------------------------
c
      ir = 1
      it = 2
      ie = 3
      if = 4
      iu = 5
      jr = ir
      jt = it
      je = ie
      jf = if
      ju = iu
      jb = jf
c
      neqn  = 4
      ngrs  = 5
      ncell = 100 + ngrs                     ! ENSMAN's computation: 290
      ngre  = ncell + 6
c
c-----------------------------------------------------------------------
c     numerics:
c     lcray = 0  => nonCRAY computer; use linpack   subroutines
c     lcray = 1  =>    CRAY computer; use cal-coded subroutines
c
c     lband = 0  => solve block pentadiagonal system
c     lband = 1  => solve banded system
c-----------------------------------------------------------------------
c
      lcray = 0
      lband = 0
      lboos = 0
c
c-----------------------------------------------------------------------
c     hydrodynamics: lhydr = 0  => no hydro
c     thet = time centering parameter
c     cadv   = order of advection. 0 => donor cell
c                                  2 => van leer
c     epsadv = overflow protection in advection switch
c-----------------------------------------------------------------------
c
      lhydr = 0
      thet  = 0.55D0 
c
        cadv = 2.0D0
      epsadv = 1.0D-50
c
c-----------------------------------------------------------------------
c     hydrodynamic boundary conditions (static medium)
c     leibc = 1  => eulerian inner bdy with zero flux
c     leobc = 1  => eulerian outer bdy with zero flux
c-----------------------------------------------------------------------
c     phil0 = mass     flux (cgs) through inner (left ) eulerian bdy
c     phil1 = momentum flux (cgs) through inner (left ) eulerian bdy
c     phil2 = energy   flux (cgs) through inner (left ) eulerian bdy
c     phir0 = mass     flux (cgs) through outer (right) eulerian bdy
c     phir1 = momentum flux (cgs) through outer (right) eulerian bdy
c     phir2 = energy   flux (cgs) through outer (right) eulerian bdy
c-----------------------------------------------------------------------
c
      llibc = 0
      llobc = 0
      leibc = 1
      leobc = 1
c
      uextl = 0.0D0
      uextr = 0.0D0
      phil0 = 0.0D0
      phil1 = 0.0D0
      phil2 = 0.0D0
      phir0 = 0.0D0
      phir1 = 0.0D0
      phir2 = 0.0D0
c
c-----------------------------------------------------------------------
c     grid specification: lgrid = 1  => adaptive   grid
c                         lgrid = 2  => eulerian   grid
c                         lgrid = 3  => lagrangean grid
c-----------------------------------------------------------------------
c     for an adaptive grid:
c     ladx    = 1  => linear      resolution
c     lady(m) = 1  => linear      resolution
c     lady(m) = 2  => logarithmic resolution
c     lady(m) = 3  => harmonic    resolution
c
c     where:
c     m = 1  => interior mass	m = 2  => density
c     m = 3  => temperature 	m = 4  => radiation energy density
c     m = 5  => pressure        m = 6  => internal  energy density
c     m = 7  => opacity         m = 8  => viscosity
c     m = 10 => velocity (bri)  m = 11 => exterior mass
c
c     "lgrid"  MUST always be set. For an adaptive grid  "ladx"  MUST be
c     set. For this problem use linear resolution in x (ladx .eq. 1) and
c     choose a scale factor "xscale".  Set  "lady(m)"  to the desired
c     option for variables actually used in defining the grid;  for all
c     other variables set lady(m) to zero. If linear resolution has been
c     chosen for a variable (lady(m) .eq. 1),   the scale factor yscl(m)
c     MUST  be set.  Choose nonzero weights  "wt(m)"  only for variables
c     used to define grid; take unity as default.Set all others to zero.
c     "yscl(1)" MUST be specified when "lady(1)" is turned on.
c-----------------------------------------------------------------------
c
           lsum  = 0
           lgrid = 1
c
      if (lgrid .eq. 1) then
c
	    alph  = 3.0D0
	    ibet  = 1
	    tau   = 1.0D3
c
	    ladx     = 1
	    lady( 1) = 0
	    lady( 2) = 0
	    lady( 3) = 0
	    lady( 4) = 2
	    lady( 5) = 0
	    lady( 6) = 0
	    lady( 7) = 0
	    lady( 8) = 0
	    lady( 9) = 0
	    lady(10) = 2
	    lady(11) = 0
c
	    if (ladx     .eq. 1) xscale   = 1.0D0
	    if (lady( 1) .gt. 0) yscl( 1) = xindef
	    if (lady( 2) .eq. 1) yscl( 2) = xindef
	    if (lady( 3) .eq. 1) yscl( 3) = 1.0D0
	    if (lady( 4) .eq. 1) yscl( 4) = 1.0D0
	    if (lady( 5) .eq. 1) yscl( 5) = xindef
	    if (lady( 6) .eq. 1) yscl( 6) = 1.0D0
	    if (lady( 7) .eq. 1) yscl( 7) = xindef 
	    if (lady( 8) .eq. 1) yscl( 8) = xindef
	    if (lady( 9) .eq. 1) yscl( 9) = xindef
	    if (lady(10) .eq. 1) yscl(10) = xlum / (4.0D0 * cpi)
	    if (lady(11) .eq. 1) yscl(11) = xindef
c
	    do 66  l = 1, mad
	    lsum = lsum + lady(l)
   66       continue
c
	    wt( 1) = 0.0D0
	    wt( 2) = 0.0D0
	    wt( 3) = 0.0D0
	    wt( 4) = 1.0D0
	    wt( 5) = 0.0D0
	    wt( 6) = 0.0D0
	    wt( 7) = 0.0D0
	    wt( 8) = 0.0D0
	    wt( 9) = 0.0D0
	    wt(10) = 2.0D0
	    wt(11) = 0.0D0
c
      end if
c
c-----------------------------------------------------------------------
c     check the switches
c-----------------------------------------------------------------------
c
      if (lgeom + lcray + lband + lgrid .lt. 0) stop 'start00'
      if (lgeom .ne. 0 .and. lgeom      .ne. 2) stop 'start01'
      if (lgrid .eq. 1 .and. lsum       .le. 0) stop 'start02'
c
c-----------------------------------------------------------------------
c     iteration control; integration control
c-----------------------------------------------------------------------
c     conv  = convergence criterion for NEWTON-RAPHSON iteration
c     ctol  = maximum fractional change in cell size in iteration
c     dtol  = maximum fractional change in physical variables in iter.
c     niter = maximum number of NEWTON-RAPHSON iteration allowed
c     ntry  = maximum number of tries for conv. using reduced timestep
c.......................................................................
c     stol  = maximum fractional change allowed over a timestep
c     nback = maximum number of reintegrations with reduced timestep if
c             maximum change > 2 * stol
c     next  = maximum number of tries to preserve monotonic mesh in
c             extrapolation to new time level
c     jdump = number of last dump file created
c     jstep = number of timestep  at dump number eq jdump
c     nstep = number of timesteps in current run
c-----------------------------------------------------------------------
c
      conv  = 1.0D-5
      ctol  = 0.70D0
      dtol  = 0.10D0
      niter = 20
      ntry  =  6
c
      stol  = 0.10D0
      nback =  6
      next  =  6
      nstep = 10
      jdump =  1
      jstep =  1
c
c-----------------------------------------------------------------------
c     output control
c          ndump = number of timesteps between dumps     of models
c          nout  = number of timesteps between printouts of models
c-----------------------------------------------------------------------
c
      ndump = 50
      nout  = 05
c
c-----------------------------------------------------------------------
c     document all the parameter settings
c-----------------------------------------------------------------------
c
      if (idoc .gt. 0) call pardoc
c
c-----------------------------------------------------------------------
c     eos: leos = 2  => perfect gas
c     set gas gamma and mean molecular weight
c-----------------------------------------------------------------------
c
      leos = 2
      gam  = 5.0D0 / 3.0D0
      gmu  = 0.5D0
c
c-----------------------------------------------------------------------
c     opacity: lopac = 3  => constant opacity
c-----------------------------------------------------------------------
c
      lopac = 3
c     ratio = 1.0D-10 is being read in from titan.in
c
c=======================================================================
c     generate initial solution
c=======================================================================
c
c++++++++ linit .eq. 1 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (linit .eq. 1) then
c
c-----------------------------------------------------------------------
c     initialize variables
c-----------------------------------------------------------------------
c.......................................................................
c     radial grid and mass distribution
c.......................................................................
c     
      dxm = xmass / ( 4.0D0 * cpi * float(ncell) )
      xmn(ngrs-1) = 1.0D0 
      rn (ngrs-1) = 1.0D6
c
      do 1 k = ngrs, ncell
      xmn(k) = xmn(ngrs-1) + dxm * float(k - ngrs + 1)
      rn (k) = ( rn(k-1)**mup1 + xmup1 * dxm / rho0 )**rxm1
    1 continue
c
      do 2 k = ncell, ncell + 7
      dxm    = 0.5D0 * dxm
      xmn(k) =  xmn(k-1) + dxm
      rn (k) = ( rn(k-1)**mup1 + xmup1 * dxm / rho0 )**rxm1
    2 continue
c
      k = ncell + 8
      ngre = k - 2
      xmn(k) =  xmn(k-1) + dxm
      rn (k) = ( rn(k-1)**mup1 + xmup1 * dxm / rho0 )**rxm1
c
      if (ladx .eq. 1) xscale  = rn(ngre+1) - rn(ngrs)
c
      do 3 k = ngrs - 1, ngre + 2
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
    3 continue
c
      do 4 k = ngrs - 1, ngre + 1
      dvoln(k) = (rmup1n(k+1) - rmup1n(k)) / xmup1
    4 continue
c
c.......................................................................
c     density 
c.......................................................................
c
      do 5 k = ngrs - 1, ngre + 2
      dn(k) = rho0
    5 continue
c
c.......................................................................
c     radiation field and temperature
c.......................................................................
c   
c     initial luminosity, surface flux, and surface radiation energy
c     density
c
      xlum0 = xlum  /  1.0D1**1.5D0
      frn0  = xlum0 / (4.0D0 * cpi * rmun(ngre+2))
      ern0  = sqrt(3.0D0) * frn0 / cc
      teff  = (ern0 / car)**0.25D0
c
c     compute initial radiation energy density assuming constant 
c     luminosity and using diffusion approximation
c
      con = 3.0D0 * xlum0 * chif0 / (4.0D0 * cpi * cc)
c
      do 6 k = ngrs - 1, ngre + 1
c
      ern (k) =  ern0 + con * (2.D0/(rn(k) + rn(k+1)) - 1.D0/rn(ngre+2))
      tn  (k) = (ern(k) / car)**0.25D0
      plfn(k) =  ern(k) * cc / ( 4.0D0 * cpi ) 
      dplfdltn(k) = 4.0D0 * plfn(k)
    6 continue
           k = ngre + 2
      ern (k) =  ern0 
      tn  (k) = (ern0   / car)**0.25D0
c
c     radiation flux. apply full luminosity at inner boundary
c     and use initial flux everywhere else
c
           k = ngrs - 1
      frn (k) = xlum  / (4.0D0 * cpi * rmun(k))
           k = ngrs
      frn (k) = xlum  / (4.0D0 * cpi * rmun(k))
      do 7 k = ngrs + 1, ngre + 2
      frn (k) = xlum0 / (4.0D0 * cpi * rmun(k))
    7 continue
c
c.......................................................................
c     opacity
c.......................................................................
c
      call opac
c
c.......................................................................
c     EDDINGTON factors
c.......................................................................
c
      call eddfac
c
c.......................................................................
c     for full transport, recalculate energy density using EDDINGTON
c     factors derived above.
c.......................................................................
c
      if (ltran .eq. 1) then
c
          ern(ngre  ) = frn(ngre) / (cc * geddr)
          ern(ngre+1) = frn(ngre) / (cc * geddr)
          ern(ngre+2) = frn(ngre) / (cc * geddr)
c
          do 8 k = ngre - 1, ngrs - 1, -1
          ern(k) = (chifn(k+1) * dn(k+1) * frn(k+1) / cc
     .              +  2.0D0 * fedd(k+1) * ern(k+1) / (rn(k+2) - rn(k))
     .      + 0.5D0 * (3.0D0 * fedd(k+1) - 1.0D0) * ern(k+1) / rn(k+1) )
     .              / (2.0D0 * fedd(k  )            / (rn(k+2) - rn(k)) 
     .      - 0.5D0 * (3.0D0 * fedd(k  ) - 1.0D0)            / rn(k+1) )
    8     continue
c
          do 9 k = ngrs - 1, ngre + 2
          tn  (k) = (ern(k) / car)**0.25D0
    9     continue
      end if
c
c.......................................................................
c     eos
c.......................................................................
c
      call eos
c
c.......................................................................
c     phantom zones
c.......................................................................
c
      delml = dn(ngrs) * dvoln(ngrs)
      delmr = dn(ngre) * dvoln(ngre)
c
      k = ngrs - 1
c
      ern(k) =   ern(k+1)
c
      k = ngre + 1
c
      ern(k  ) =   ern(k-1)
      ern(k+1) =   ern(k  )
c
c.......................................................................
c     zero unneeded variables
c.......................................................................
c
      do 10 k = ngrs - 1, ngre + 1
      dsn   (k) = 0.0D0
      egsn  (k) = 0.0D0
      etn   (k) = 0.0D0
      etsn  (k) = 0.0D0
      ersn  (k) = 0.0D0
      egtn  (k) = 0.0D0
      egtsn (k) = 0.0D0
      egrsn (k) = 0.0D0
      egrtn (k) = 0.0D0
      egrtsn(k) = 0.0D0
      frsn  (k) = 0.0D0
   10 continue
c
      do 11 k = ngrs - 1, ngre + 2
      r0   (k) = rn (k)
      xm0  (k) = xmn(k)
      un   (k) = 0.0D0
      usn  (k) = 0.0D0
      urel (k) = 0.0D0
      dmdt (k) = 0.0D0
      drdt (k) = 0.0D0
      frnom(k) = frn(k)
       unom(k) = 0.0D0
   11 continue
c
c.......................................................................
c     zero initial errors 
c.......................................................................
c
      cmax = 0.0D0
      dmax = 0.0D0
      smax = 0.0D0
c
c.......................................................................
c     initialize old time solution
c.......................................................................
c
      do 12 k = ngrs - 1, ngre + 2
      ro   (k) = rn(k)
      to   (k) = tn(k)
      do   (k) = dn(k)
   12 continue
c
c.......................................................................
c     initialize time and timestep
c.......................................................................
c
      timen = 0.0D0
      timeo = 0.0D0
      tfac  = 1.0D0
c
c     need to choose a reasonable value for first timestep. 
c     diffusion time = 1.6e+7 sec.
c
      dtime = 2.0D-10
c
c.......................................................................
c     initialize grid equation
c.......................................................................
c
      rw         = 5.00D0 * crsol
      radc       = rn(ngrs) + 2.D0 * rw
      rn(ngrs-1) = rn(ngrs) - 2.D0 * rw
c
      if (lgrid .eq. 1) call grid6nit (xlum, xlum0, radc, rw)
c
c-----------------------------------------------------------------------
c     radial grid and mass distribution
c-----------------------------------------------------------------------
c
      do 13 k = ngrs - 1, ngre + 2
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
   13 continue
c
      do 14 k = ngrs - 1, ngre + 1
      dvoln(k) = (rmup1n(k+1) - rmup1n(k)) / xmup1
   14 continue
c
      xmn(ngrs-1) = 1.0D0 
      do 15 k = ngrs - 1, ngre + 1
      xmn (k) = xmn(k-1) + rho0 * dvoln(k-1)
   15 continue
c
c.......................................................................
c     exterior mass
c.......................................................................
c
      xmass = xmn(ngre+2)
      xmext = 1.D-10
      xmtot = xmass + xmext
c
      do 155 k = ngrs - 1, ngre + 2
      xmen(k) = xmtot - xmn(k)
  155 continue
c
c-----------------------------------------------------------------------
c     radiation field and temperature
c-----------------------------------------------------------------------
c
      do 16 k = ngrs - 1, ngre + 1
c
      ern (k) =  ern0 + con * (2.D0/(rn(k) + rn(k+1)) - 1.D0/rn(ngre+2))
      tn  (k) = (ern(k) / car)**0.25D0
      plfn(k) =  ern(k) * cc / ( 4.D0 * cpi ) 
      dplfdltn(k) = 4.0D0 * plfn(k)
   16 continue
           k = ngre + 2
      ern (k) =  ern0 
      tn  (k) = (ern0   / car)**0.25D0
c
c-----------------------------------------------------------------------
c     opacity
c-----------------------------------------------------------------------
c
      call opac
c
c-----------------------------------------------------------------------
c     EDDINGTON factors
c-----------------------------------------------------------------------
c
      call eddfac
c
c-----------------------------------------------------------------------
c     for full transport, recalculate energy density using EDDINGTON
c     factors derived above.
c-----------------------------------------------------------------------
c
      if (ltran .eq. 1) then
c
          ern(ngre  ) = frn(ngre) / (cc * geddr)
          ern(ngre+1) = frn(ngre) / (cc * geddr)
          ern(ngre+2) = frn(ngre) / (cc * geddr)
c
          do 17 k = ngre - 1, ngrs - 1, -1
          ern(k) = (chifn(k+1) * dn(k+1) * frn(k+1) / cc
     .              +  2.0D0 * fedd(k+1) * ern(k+1) / (rn(k+2) - rn(k))
     .      + 0.5D0 * (3.0D0 * fedd(k+1) - 1.0D0) * ern(k+1) / rn(k+1) )
     .              / (2.0D0 * fedd(k  )            / (rn(k+2) - rn(k)) 
     .      - 0.5D0 * (3.0D0 * fedd(k  ) - 1.0D0)            / rn(k+1) )
   17     continue
c
          do 18 k = ngrs - 1, ngre + 2
          tn  (k) = (ern(k) / car)**0.25D0
   18     continue
      end if
c
c-----------------------------------------------------------------------
c     EDDINGTON factors
c-----------------------------------------------------------------------
c
      call eddfac
c
c-----------------------------------------------------------------------
c     eos
c-----------------------------------------------------------------------
c
      call eos
c
      do 19 k = ngrs - 1, ngre + 2
      egrn (k) =  ern(k) / dn(k) + egn(k)
   19 continue
c
c-----------------------------------------------------------------------
c     phantom zones
c-----------------------------------------------------------------------
c
      delml = dn(ngrs) * dvoln(ngrs)
      delmr = dn(ngre) * dvoln(ngre)
c
      k = ngrs - 1
c
      ern(k) =   ern(k+1)
      frn(k) = ( rn (k+1) / rn(k) )**mu * frn(k+1)
c
      k = ngre + 1
c
      ern(k  ) =   ern(k-1)
      ern(k+1) =   ern(k  )
      frn(k+1) = ( rn (k) / rn(k+1) )**mu * frn(k)
c 
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 30 k = ngrs, ngre + 1
      xm0   (k) = xmn(k)
      r0    (k) = rn (k)
      u0    (k) = un (k)
      tescr (k) = dn(k) * dvoln(k) * egrn(k)
      tescr0(k) = 0.0D0
   30 continue
c
      tee   = 0.0D0
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
c
c-----------------------------------------------------------------------
c     zero unneeded variables
c-----------------------------------------------------------------------
c
      do 31  k = ngrs - 1, ngre + 1
      dsn   (k) = 0.0D0
      egsn  (k) = 0.0D0
      etn   (k) = 0.0D0
      etsn  (k) = 0.0D0
      ersn  (k) = 0.0D0
      egtn  (k) = 0.0D0
      egtsn (k) = 0.0D0
      egrsn (k) = 0.0D0
      egrtn (k) = 0.0D0
      egrtsn(k) = 0.0D0
      frsn  (k) = 0.0D0
   31 continue
c
      do 32 k = ngrs - 1, ngre + 2
      un   (k) = 0.0D0
      usn  (k) = 0.0D0
      urel (k) = 0.0D0
      dmdt (k) = 0.0D0
      drdt (k) = 0.0D0
      frnom(k) = frn(k)
       unom(k) = 0.0D0
   32 continue
c
c-----------------------------------------------------------------------
c     initialize time and timestep
c-----------------------------------------------------------------------
c
      timen = 0.0D0
      timeo = 0.0D0
      tfac  = 1.0D0
c
      dtime = 2.0D-5
c
c-----------------------------------------------------------------------
c     zero initial errors 
c-----------------------------------------------------------------------
c
      cmax = 0.0D0
      dmax = 0.0D0
      smax = 0.0D0
c
c=======================================================================
c     write initial model on dump file
c=======================================================================
c
      if (idump .gt. 0) then
                               irec = 1
          write ( idump, rec = irec, iostat = ier) 
     .                            header, xmass, xlum, rho0, chif0, ximr
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start05' 
          end if
      end if
c
      read (iin,'(3i5,e10.3)') jdump, jstep, nstep, tmax
c
      jsteps = jstep + 1
      jstepe = jsteps + nstep - 1
      if (jdump .ne. 1 .or. jstep .ne. 1) stop 'start06'
c
      write (itty, 20) jdump, jstep
      write (iout, 20) jdump, jstep
c
      if (idump .gt. 0) then
                               irec = 2
          write ( idump, rec = irec, iostat = ier) jdump, jstep
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start07' 
          end if
                               irec = jdump + 2
          write ( idump, rec = irec, iostat = ier)
     .          jdump , jstep , cmax  , dmax  , smax  , timen , dtime , 
     .          tfac  , tetot , tee   , tes   , tew   , tel   , teq   ,
     .          r0    , xm0   , d0    , u0    , egrt0 , dvol0 , tescr0,
     .          geddl , geddr , xmass , xmext , xmtot ,
     .          rn    , rmun  , rmup1n, rmum1n, avchin, urel  ,
     .          dvoln , xmen  , dn    , dsn   , chifn , drdt  ,
     .          un    , usn   , tn    , asn   , xken  , dmdt  ,
     .          egn   , egsn  , pgn   , xnen  , xkpn  , xnt   ,
     .          ern   , ersn  , egrn  , egrsn , plfn  , 
     .          etn   , etsn  , egtn  , egtsn , fedd  ,
     .          frn   , frsn  , egrtn , egrtsn
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start08' 
          end if
      end if
c
      call writout
c
      return
      end if
c
c+++++++++++++++++++++++++++++++ linit .eq. 1 +++++++++++++++++++++++end
c
c=======================================================================
c     recover model from dump file
c=======================================================================
c
c++++++++ linit .ne. 1 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (linit .ne. 1) then
c
c-----------------------------------------------------------------------
c     compare model header
c-----------------------------------------------------------------------
c
      if (idump .gt. 0) then
                               irec = 1
          read  ( idump, rec = irec, iostat = ier) 
     .                         yheader, ymass, ylum, yrho0, ychif0, ximr
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start09' 
          end if
c
          if (yheader .ne. header) then
              write(itty, '(a80/a80)') header, yheader
              write(iout, '(a80/a80)') header, yheader
              stop 'start10' 
          end if
c
          if (ymass  .ne. xmass .or. ylum  .ne. xlum .or.
     .        ychif0 .ne. chif0 .or. yrho0 .ne. rho0 ) then
              write(itty,'(1p8e10.3)') 
     .              xmass, ymass, xlum, ylum, rho0, yrho0, chif0, ychif0
              write(iout,'(1p8e10.3)') 
     .              xmass, ymass, xlum, ylum, rho0, yrho0, chif0, ychif0
              stop 'start11' 
          end if
      end if
c
c-----------------------------------------------------------------------
c     read dump number and timestep number of initial model; check that
c     specified starting model exists
c-----------------------------------------------------------------------
c
      read (iin,'(3i5,e10.3)') jdumpi, jstepi, nstep, tmax
      jdump  = jdumpi 
      jsteps = jstepi + 1
      jstepe = jsteps + nstep - 1
c
      if (idump .gt. 0) then
                               irec = 2
          read  ( idump, rec = irec, iostat = ier) jdumpo, jstepo
c
          if (ier .gt. 0) then
    	      write (itty, 21) irec, ier, ier, ier
    	      write (iout, 21) irec, ier, ier, ier
	      stop 'start12' 
          end if
c
          if (jdumpo .lt. jdumpi .or. jstepo .lt. jstepi) then
    	      write (itty, 22) jdumpo, jdumpi, jstepo, jstepi
	      stop 'start13'
          end if
c
c-----------------------------------------------------------------------
c     read model
c-----------------------------------------------------------------------
c
                               irec = jdumpi + 2
          read  ( idump, rec = irec, iostat = ier)
     .          jdumpd, jstepd, cmax  , dmax  , smax  , timen , dtime ,
     .          tfac  , tetot , tee   , tes   , tew   , tel   , teq   ,
     .          r0    , xm0   , d0    , u0    , egrt0 , dvol0 , tescr0,
     .          geddl , geddr , xmass , xmext , xmtot ,
     .          rn    , rmun  , rmup1n, rmum1n, avchin, urel  ,
     .          dvoln , xmen  , dn    , dsn   , chifn , drdt  ,
     .          un    , usn   , tn    , asn   , xken  , dmdt  ,
     .          egn   , egsn  , pgn   , xnen  , xkpn  , xnt   ,
     .          ern   , ersn  , egrn  , egrsn , plfn  , 
     .          etn   , etsn  , egtn  , egtsn , fedd  ,
     .          frn   , frsn  , egrtn , egrtsn
c
          xlum0 = xlum  /  1.0D1**1.5D0
          frn0  = xlum0 / (4.0D0 * cpi * rmun(ngre+2))
          ern0  = sqrt(3.0D0) * frn0 / cc
          teff  = (ern0 / car)**0.25D0
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start14' 
          end if
c
c-----------------------------------------------------------------------
c     check that we got the right one
c-----------------------------------------------------------------------
c
          if (jdumpd .ne. jdumpi .or. jstepd .ne. jstepi) then
    	      write (itty, 22) jdumpd, jdumpi, jstepd, jstepi
              stop 'start15'
          end if
      end if
c
      return
      end if
c
c+++++++++++++++++++++++++++++++ linit .ne. 1 +++++++++++++++++++++++end
c
c=======================================================================
c
   20 format(/' model jdump ='i6' jstep ='i6' archived')
   21 format(/' unit idump. irec = 'i3' ier ='i20, o22, a8)
   22 format( ' jdumpo ='i7' jdumpi ='i7' jstepo ='i7' jstepi ='i7)
c
      end
