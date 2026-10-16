      subroutine start10
c
c***********************************************************************
c     setup for RADIATING BLAST WAVE problem
c.......................................................................
c     calling sequence: start > start10 > eos, writout, gridinit
c***********************************************************************
c     comparison results  exist in the  literature for times  t = 3.0e3,
c     4.0e3,  5.0e3, 5.5e3, 5.8e3,  6.0e3, 6.2e3, 6.8e3,  8.02e3, 1.0e4,
c     1.1e4, 1.2e4, 1.7e4, 2.0e4, and 2.2e4 seconds.
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
      dimension xx(mgr)
c
c=======================================================================
c     general setup
c=======================================================================
c-----------------------------------------------------------------------
c     read header and initial values
c-----------------------------------------------------------------------
c     rmin = inner radius of shell (cm)
c     rmax = outer radius of shell (cm)
c     xlum = imposed luminosity (erg/s)
c     rho0 = density (gm/cm**3)
c     xmin = mass of outermost zone
c     ecnt = energy deposit at inner bnd
c-----------------------------------------------------------------------
c
      read (iin,'(a80)') header
      read (iin,'(3e10.3,3e10.1)') rmin, rmax, xlum, rho0, xmin, ecnt
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
      cgrav = 0.0D0              ! geometry is spherical; gravity is off
       g    = 0.0D0              ! that's ok too
      xmint = 0.0D0
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
c     declare equations to be solved
c     ir > 0 => grid         eqn
c     im > 0 => mass defin.  eqn;  id > 0 => continuity eqn
c     iu > 0 => gas momentum eqn;  it > 0 => gas energy eqn
c     if > 0 => rad momentum eqn;  ie > 0 => rad energy eqn
c
c     neqn  = number of equations at each grid-point
c     ngrs  = index of first (left /inner-most) cell
c     ngre  = index of last  (right/outer-most) cell
c     ncell = number of cells in the domain
c-----------------------------------------------------------------------
c
      ir = 1
      im = 2
      id = 3
      iu = 4
      it = 5
      ie = 6
      if = 7
      jr = ir
      jm = im
      jd = id
      ju = iu
      jt = it
      je = ie
      jf = if
      jb = iu
c
      neqn  = 7
      ngrs  = 5
      ncell = 100                            ! Ensman's computation: 350
      ngre  = ngrs + ncell - 1
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
c     set hydro parameters: lhydr = 1  =>    hydro
c     thet   = time centering parameter
c     ql0    = fixed length for pseudoviscosity (planar)
c     ql1    = fixed relative length for pseudoviscosity (spherical)
c     cq1    = coefficient of linear     pseudoviscosity
c     cq2    = coefficient of quadratic  pseudoviscosity
c     cqvis  = coefficient of artificial stress tensor
c     q0     = floor on pseudoviscous pressure ratio
c     cadv   = order of advection. 0 => donor cell
c                                  2 => van leer
c     epsadv = overflow protection in advection switch
c     sigd   = diffusion coefficient in continuity eqn
c     sige   = diffusion coefficient in gas energy eqn
c-----------------------------------------------------------------------
c
      lhydr = 1
      thet  = 1.0D0
c
      ql0   = 0.0D0
      ql1   = 1.0D-2
      cq2   = 4.0D0
      cq1   = 0.0D0
c*****cq1   = 0.1D0 * cq2
      cqvis = 4.0D0 / 3.0D0
      q0    = 1.0D-30
c
      cadv   = 2.0D0
      epsadv = 1.0D-50
c
      sigd = 0.0D0
      sige = 0.0D0
c
c-----------------------------------------------------------------------
c     hydrodynamic boundary conditions 
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
      llobc = 1
      leibc = 1
      leobc = 0
c
      phil0 = 0.0D0
      phil1 = 0.0D0
      phil2 = 0.0D0
      phir0 = 0.0D0
      phir1 = 0.0D0
      phir2 = 0.0D0
c
      ul    = 0.0D0
      uextl = 0.0D0
      uextr = 0.0D0
      pextr = 0.0D0
      omega = 1.0D30
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
	    alph  = 1.5D0
	    ibet  = 1
	    tau   = 1.0D-20
            ql1   = 1.0D-4
c
	    ladx     = 1
	    lady( 1) = 2
	    lady( 2) = 2
	    lady( 3) = 0
	    lady( 4) = 2
	    lady( 5) = 0
	    lady( 6) = 0
	    lady( 7) = 0
	    lady( 8) = 0
	    lady( 9) = 0
	    lady(10) = 1
	    lady(11) = 0
c
            if (ladx     .eq. 1) xscale   = rmax 
	    if (lady( 1) .gt. 0) yscl( 1) = 0.0D0
	    if (lady( 2) .eq. 1) yscl( 2) = xindef
	    if (lady( 3) .eq. 1) yscl( 3) = 1.0D0
	    if (lady( 4) .eq. 1) yscl( 4) = 1.0D0
	    if (lady( 5) .eq. 1) yscl( 5) = xindef
	    if (lady( 6) .eq. 1) yscl( 6) = 1.0D0
	    if (lady( 7) .eq. 1) yscl( 7) = xindef 
	    if (lady( 8) .eq. 1) yscl( 8) = 1.0D0
	    if (lady( 9) .eq. 1) yscl( 9) = xindef
	    if (lady(10) .eq. 1) yscl(10) = 1.0D8
	    if (lady(11) .eq. 1) yscl(11) = xindef
c
	    do 66  l = 1, mad
	    lsum = lsum + lady(l)
   66       continue
c
	    wt( 1) = 1.0D0
	    wt( 2) = 1.0D0
	    wt( 3) = 0.0D0
	    wt( 4) = 1.0D0
	    wt( 5) = 0.0D0
	    wt( 6) = 0.0D0
	    wt( 7) = 0.0D0
	    wt( 8) = 0.0D0
	    wt( 9) = 0.0D0
	    wt(10) = 1.0D0
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
      niter = 15
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
      ndump = 100
      nout  = 100
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
      gam = 5.0D0 / 3.0D0
      gmu = 0.5D0
c
c-----------------------------------------------------------------------
c     opacity: lopac = 4  => constant opacity
c-----------------------------------------------------------------------
c
      lopac = 4
      ratio = 1.0D0
c
c=======================================================================
c     generate initial solution
c=======================================================================
c
c++++++++ linit .eq. 1 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (linit .eq. 1) then
c-----------------------------------------------------------------------
c     relax to initial temperature profile
c-----------------------------------------------------------------------
c.......................................................................
c     set up radiative transport only
c.......................................................................
c     declare equations to be solved (radiation only in static medium)
c     ir > 0 => grid         eqn;  it > 0 => gas energy eqn 
c     if > 0 => rad momentum eqn;  ie > 0 => rad energy eqn 
c.......................................................................
c
      iro = ir
      imo = im
      ido = id
      iuo = iu
      ito = it
      ieo = ie
      ifo = if
c
      ir = 1
      it = 2
      ie = 3
      if = 4
      id = 5
      im = 6
      iu = 7
      jr = ir
      jm = im
      jd = id
      ju = iu
      jt = it
      je = ie
      jf = if
c
      neqo = neqn
      neqn = 4
c
c.......................................................................
c     radiation transport: lrad  = 1  => radiation
c.......................................................................
c
      lrado = lrad
      lrad  = 1
c
c.......................................................................
c     hydrodynamics: lhydr = 0  => no hydro
c     thet = time centering parameter: thet = 1 => fully implicit 
c.......................................................................
c
      lhydro = lhydr
      lhydr  = 0
      theto  = thet
      thet   = 1.0D0
c
c.......................................................................
c     grid specification: lgrid = 3  => lagrangean grid
c.......................................................................
c
      lgrido = lgrid
      lgrid  = 3
c
c.......................................................................
c     radial grid and mass distribution
c.......................................................................
c     
      xmtot = rho0 * rmin**mup1                   ! total mass is 1.0D30
c
      rn (ngre+2) = rmax  * 1.01D0
      rn (ngre+1) = rmax
      xmn(ngre+1) = xmtot * 0.25D0* ( 1.D0 - (rmin/rmax)**4 )
      xx (ngre  ) = xmin  / (4.D0 * cpi)
      xx (ngre+1) = 0.0D0
c
                        f = 1.125D0
      do 1 k = ngre, ngrs - 1, -1
      xmn(k  ) = xmn(k+1) - xx(k)
      xx (k-1) =        f * xx(k)
      rn (k  ) = 1.0D0 / rmin**4 * ( 1.D0 - 4.D0 * xmn(k) / xmtot )
      rn (k  ) = 1.0D0 / rn(k)**0.25D0
    1 continue
c
      xmn(ngrs-1) = 1.0D0
c
      do 2 k = ngrs - 1, ngre + 2
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
    2 continue
c
      do 3 k = ngrs - 1, ngre + 1
      dvoln(k) = (rmup1n(k+1) - rmup1n(k)) / xmup1
    3 continue
c
c.......................................................................
c     density
c.......................................................................
c
      dn(ngre+1) = rho0 * (rmin / rmax)**7
c
      do 4 k  = ngre, ngrs, -1
      rmid    =   0.5D0 * (rn(k) + rn(k+1))
      dn  (k) =    rho0 * (rmin / rmid)**7
    4 continue
c
c.......................................................................
c     velocity
c.......................................................................
c
      do 5 k = ngrs - 1, ngre + 2
      un (k) = 0.0D0
    5 continue
c
c.......................................................................
c     radiation field and temperature
c.......................................................................
c   
c     initial luminosity, surface flux, and surface radiation energy
c     density
c
      frn0  = xlum / (4.0D0 * cpi * rmun(ngre+2))
      ern0  = sqrt(3.0D0) * frn0 / cc
      teff  = (ern0 / car)**0.25D0
c
c     compute initial radiation energy density assuming constant 
c     luminosity and using diffusion approximation
c
      con = 3.0D0 * xlum * csige / (4.0D0 * cpi * cc * cm0)
c
      do 6 k = ngrs - 1, ngre + 1
c
      ern (k) =  ern0 + con * dn(k) * ( 2.D0/(rn(k) + rn(k+1))
     .                                - 1.D0/ rn(ngre+2)     )
      tn  (k) = (ern(k) / car)**0.25D0
      plfn(k) =  ern(k) * cc / ( 4.0D0 * cpi ) 
      dplfdltn(k) = 4.0D0 * plfn(k)
    6 continue
           k = ngre + 2
      ern (k) =  ern0 
      tn  (k) = (ern0   / car)**0.25D0
c
c     radiation flux
c
      do 7 k = ngrs - 1, ngre + 2
      frn (k) = xlum / (4.0D0 * cpi * rmun(k))
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
      ern(ngre  ) = frn(ngre) / (cc * geddr)
      ern(ngre+1) = frn(ngre) / (cc * geddr)
      ern(ngre+2) = frn(ngre) / (cc * geddr)
c
      do 8 k = ngre - 1, ngrs - 1, -1
      ern(k) = (chifn(k+1) * dn(k+1) * frn(k+1) / cc
     .          +  2.0D0 * fedd(k+1) * ern(k+1) / (rn(k+2) - rn(k))
     .  + 0.5D0 * (3.0D0 * fedd(k+1) - 1.0D0) * ern(k+1) / rn(k+1) )
     .          / (2.0D0 * fedd(k  )            / (rn(k+2) - rn(k)) 
     .  - 0.5D0 * (3.0D0 * fedd(k  ) - 1.0D0)            / rn(k+1) )
    8 continue
c
      do 9 k = ngrs - 1, ngre + 2
      tn  (k) = (ern(k) / car)**0.25D0
    9 continue
c
c.......................................................................
c     phantom zones
c.......................................................................
c
      k = ngrs - 1
c
      tn (k) =   tn (k+1)
      dn (k) =   dn (k+1)
      ern(k) =   ern(k+1)
      frn(k) = ( rn (k+1) / rn(k) )**mu * frn(k+1)
c
      k = ngre + 1
c
      tn (k  ) =   tn (k-1)
      tn (k+1) =   tn (k  )
      dn (k+1) =   dn (k  )
      ern(k  ) =   ern(k-1)
      ern(k+1) =   ern(k  )
      frn(k+1) = ( rn (k) / rn(k+1) )**mu * frn(k)
      xmn(k+1) =   xmn(k) + delmr
c
      delml =  dn(ngrs) * dvoln(ngrs)
      delmr =  dn(ngre) * dvoln(ngre)
c
      geddl = 0.0D0
      geddr = 0.0D0
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
c     eos
c.......................................................................
c
      call eos
c
      do 10 k = ngrs - 1, ngre + 2
      egrn (k) = ern(k) / dn(k) + egn(k)
   10 continue
c
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 39 k = ngrs, ngre + 1
      xm0   (k) = xmn  (k)
      r0    (k) = rn   (k)
      d0    (k) = dn   (k)
      u0    (k) = un   (k)
      egrt0 (k) = egrn (k) + etn(k)
      dvol0 (k) = dvoln(k)
      tescr (k) = 0.0D0
      tescr0(k) = 0.0D0
   39 continue
c
      tee   = 0.0D0
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
c
c.......................................................................
c     initialize old time solution
c.......................................................................
c
      do 12 k = ngrs - 1, ngre + 2
      ro   (k) = rn(k)
      to   (k) = tn(k)
      do   (k) = dn(k)
      xnto (k) = 0.0D0
   12 continue
c
c.......................................................................
c     zero unneeded variables
c.......................................................................
c
      do 13  k = ngrs - 1, ngre + 1
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
   13 continue
c
      do 14 k = ngrs - 1, ngre + 2
      usn  (k) = 0.0D0
      urel (k) = 0.0D0
      dmdt (k) = 0.0D0
      drdt (k) = 0.0D0
      frnom(k) = frn(k)
       unom(k) = 1.0D0
   14 continue
c
c.......................................................................
c     initialize time and timestep
c.......................................................................
c
      timen = 0.0D0
      timeo = 0.0D0
      tfac  = 1.0D0
c
      dtime = 1.0D0
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
c     relaxation to initial temperature profile
c.......................................................................
c
      write(itty,'(/" relaxation to initial temperature profile")')
c
      do 15 jstep = 1,50
      call   step
   15 continue
c
c-----------------------------------------------------------------------
c     set up for radiation hydrodynamics
c-----------------------------------------------------------------------
c.......................................................................
c     declare equations to be solved
c     ir > 0 => grid         eqn
c     im > 0 => mass defin.  eqn;  id > 0 => continuity eqn
c     iu > 0 => gas momentum eqn;  it > 0 => gas energy eqn
c     if > 0 => rad momentum eqn;  ie > 0 => rad energy eqn
c
c     neqn  = number of equations at each grid-point
c.......................................................................
c
      ir = iro
      im = imo
      id = ido
      iu = iuo
      it = ito
      ie = ieo
      if = ifo
      jr = ir
      jm = im
      jd = id
      ju = iu
      jt = it
      je = ie
      jf = if
c
      neqn  = neqo
c
c.......................................................................
c     radiation transport
c.......................................................................
c
      lrad = lrado
c
c.......................................................................
c     hydrodynamics
c.......................................................................
c
      lhydr = lhydro
      thet  = theto
c
c.......................................................................
c     grid specification
c.......................................................................
c
      lgrid  = lgrido
c
c-----------------------------------------------------------------------
c     zero unneeded variables
c-----------------------------------------------------------------------
c
      do 17  k = ngrs - 1, ngre + 1
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
   17 continue
c
      do 18 k = ngrs - 1, ngre + 2
      usn  (k) = 0.0D0
      urel (k) = 0.0D0
      dmdt (k) = 0.0D0
      drdt (k) = 0.0D0
      frnom(k) = frn(k)
       unom(k) = 1.0D0
   18 continue
c
c-----------------------------------------------------------------------
c     radial grid
c-----------------------------------------------------------------------
c
      do 19 k = ngrs - 1, ngre + 1
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
   19 continue
c
      do 30 k = ngrs - 1, ngre
      dvoln(k) = ( rmup1n(k+1) - rmup1n(k) ) / xmup1
   30 continue
c
c-----------------------------------------------------------------------
c     density
c-----------------------------------------------------------------------
c
      dn(ngre+1) = rho0 * (rmin / rmax)**7
c
      do 31 k  = ngre, ngrs, -1
      rmid     =   0.5D0 * (rn(k) + rn(k+1))
      dn   (k) =    rho0 * (rmin / rmid)**7
   31 continue
c
c-----------------------------------------------------------------------
c     mass and velocity
c-----------------------------------------------------------------------
c
      xmn(ngrs-1) = 0.0D0
      un (ngrs  ) = 0.0D0
c
      do 32 k = ngrs, ngre + 1
      xmn  (k) = xmn(k-1) + dn(k-1) * dvoln(k-1)
      un   (k) = 0.0D0
   32 continue
c
      delml = xmn(ngrs)
      delmr =  dn(ngre) * dvoln(ngre)
      xmn(ngre+2) = xmn(ngre+1) + delmr
c
c.......................................................................
c     exterior mass
c.......................................................................
c
      xmass = xmn(ngre+2)
      xmext = 1.0D-10
      xmtot = xmass + xmext
      xmen(ngre+2) = 1.D-30
c
      do 155 k = ngrs - 1, ngre + 1
      xmen(k) = xmtot - xmn(k)
  155 continue
c
c-----------------------------------------------------------------------
c     dump internal energy in first zone; recalculate temperature
c     material energy, radiation energy, and pressure
c
      etot = ecnt / ( dn(ngrs) * dvoln(ngrs) )
      write(itty,'("energy deposition   ="1pe12.5)')etot
      write(itty,'("mass of inner shell ="1pe12.5)')dn(ngrs)*dvoln(ngrs)
c-----------------------------------------------------------------------
c
c     NEWTON-RAPHSON iteration to get correct temperature
c
      tn (ngrs) = (etot * dn(ngrs) / car)**0.25D0
      do 16 kount = 1,5
      egn(ngrs) = (ck / (cm0 * gmu)) * tn(ngrs) / (gam - 1.D0)
      ern(ngrs) =              car   * tn(ngrs)**4
      delt = (etot - egn(ngrs) - ern(ngrs) / dn(ngrs))
     .                      / (egn(ngrs) + 4.0D0 * ern(ngrs) / dn(ngrs))
      tn(ngrs) = (1.0D0 + delt) * tn(ngrs)
   16 continue
c
      pgn (ngrs) = (ck / (cm0 * gmu)) * tn(ngrs) * dn(ngrs) 
      egn (ngrs) = (ck / (cm0 * gmu)) * tn(ngrs) / (gam - 1.D0)
      ern (ngrs) = car * tn(ngrs)**4
c
c-----------------------------------------------------------------------
c     radiation field
c-----------------------------------------------------------------------
c   
      do 33 k = ngrs, ngre + 2
      ern (k) =  car * tn(k)**4
      plfn(k) =  ern(k) * cc / ( 4.0D0 * cpi ) 
      dplfdltn(k) = 4.0D0 * plfn(k)
   33 continue
c
c     radiation flux
c
      do 34 k = ngrs - 1, ngre + 2
      frn (k) = xlum / (4.0D0 * cpi * rmun(k))
   34 continue
c
      call opac
c
      call eddfac
c
      ern(ngre  ) = frn(ngre) / (cc * geddr)
      ern(ngre+1) = frn(ngre) / (cc * geddr)
      ern(ngre+2) = frn(ngre) / (cc * geddr)
c
C     do 35 k = ngre - 1, ngrs - 1, -1
C     ern(k) = (chifn(k+1) * dn(k+1) * frn(k+1) / cc
C    .          +  2.0D0 * fedd(k+1) * ern(k+1) / (rn(k+2) - rn(k))
C    .  + 0.5D0 * (3.0D0 * fedd(k+1) - 1.0D0) * ern(k+1) / rn(k+1) )
C    .          / (2.0D0 * fedd(k  )            / (rn(k+2) - rn(k)) 
C    .  - 0.5D0 * (3.0D0 * fedd(k  ) - 1.0D0)            / rn(k+1) )
C  35 continue
c
c-----------------------------------------------------------------------
c     phantom zones
c-----------------------------------------------------------------------
c
      k = ngrs - 1
      tn (k  ) =  tn (k+1)
      dn (k  ) =  dn (k+1)
      un (k  ) =  un (k+1)
      egn(k  ) =  egn(k+1)
      xmn(k  ) =  xmn(k+1) - dn(k)/xmup1 * (rn(k+1)**mup1 -
     *                                      rn(k  )**mup1)
c
      k = ngre + 1
      tn (k  ) =   tn (k-1)
      tn (k+1) =   tn (k  )
      dn (k  ) =   dn (k-1)
      dn (k+1) =   dn (k  )
      rn (k+1) = ( rn (k  )**mup1 + xmup1*delmr/dn(k) )**rxm1
      xmn(k+1) =   xmn(k  )             + delmr
c
      geddl = 0.0D0
      geddr = 0.0D0
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
c     eos
c-----------------------------------------------------------------------
c
      call eos
c
      do 36 k = ngrs - 1, ngre + 2
      egrn (k) = ern(k) / dn(k) + egn(k)
   36 continue
c
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 37 k = ngrs, ngre
      tescr (k) = dn(k) * dvoln(k) * egrn(k)
      tescr0(k) = 0.0D0
   37 continue
c
      tee   = ssum(ngre - ngrs + 1, tescr(ngrs), 1)
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
      write(iout, '(/" total energy =" 1pe13.6)') tetot
      write(itty, '(/" total energy =" 1pe13.6)') tetot
c
c-----------------------------------------------------------------------
c     initialize old time solution
c-----------------------------------------------------------------------
c
      do 38 k = ngrs - 1, ngre + 2
      ro   (k) = rn(k)
      to   (k) = tn(k)
      do   (k) = dn(k)
      xnto (k) = 0.0D0
   38 continue
c
c-----------------------------------------------------------------------
c     initialize time and timestep
c-----------------------------------------------------------------------
c
      timen = 0.0D0
      timeo = 0.0D0
      tfac  = 1.0D0
c
      dtime = 1.0D-25
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
     .                      header, xmtot, rmin, rmax, xlum, rho0, xmin
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
      irec = 2
      if (idump .gt. 0) then
          write (idump, rec = irec, iostat = ier) jdump, jstep
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start07' 
          end if
      end if
c
      irec = jdump + 2
      if (idump .gt. 0) then
          write (idump, rec = irec, iostat = ier)
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
      irec = 1
      if (idump .gt. 0) then
          read (idump, rec = irec, iostat = ier)
     .                  yheader, xmtot, yrmin, yrmax, ylum, yrho0, ymin
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
          if (ylum .ne. xlum .or.
     .        ymin .ne. xmin .or. yrho0 .ne. rho0) then
              write(itty,'(1p6e10.3)') 
     .                               xmin, ymin, xlum, ylum, rho0, yrho0
              write(iout,'(1p6e10.3)') 
     .                               xmin, ymin, xlum, ylum, rho0, yrho0
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
      irec = 2
      if (idump .gt. 0) then
          read (idump, rec = irec, iostat = ier) jdumpo, jstepo
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
      end if
c
c-----------------------------------------------------------------------
c     read model
c-----------------------------------------------------------------------
c
      irec = jdumpi + 2
      if (idump .gt. 0) then
          read (idump, rec = irec, iostat = ier)
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
          frn0  = xlum / (4.0D0 * cpi * rmun(ngre+2))
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
