      subroutine start04
c
c***********************************************************************
c     setup for WOODWARD-COLLELA BLAST WAVE problem
c.......................................................................
c     calling sequence: start > start04 > pardoc, eos, writout, gridinit
c***********************************************************************
c     comparison results  exist in literature for times t = 0.01, 0.016, 
c     0.026, 0.028, 0.030, 0.032, 0.034, and 0.038 seconds.
c
c     AA, 266, 266-292 (1992): 400 grid points; alph = 2; tau = 0.1;
c                              rho0 = 1; pg0 = {1.e+3; 1.e-2; 1.e+2};
c                              solution at time t=0.03826 yields 
c                              contact at x=0.59, refl. shock at x=0.65,
c                              contact at x=0.76, contact at x=0.78, and
c                              refl. shock at x=0.85.
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
c     read in header
c-----------------------------------------------------------------------
c
      read (iin,'(a80)') header
c
c-----------------------------------------------------------------------
c     define parameters ("l" => "left"; "r" => "right")
c-----------------------------------------------------------------------
c     ydl1 = initial downstream density  at first  discontinuity
c     ypl1 = initial downstream pressure at first  discontinuity
c     yul1 = initial downstream velocity at first  discontinuity
c     ydr1 = initial upstream   density  at first  discontinuity
c     ypr1 = initial upstream   pressure at first  discontinuity
c     yur1 = initial upstream   velocity at first  discontinuity
c     ydl2 = initial downstream density  at second discontinuity
c     ypl2 = initial downstream pressure at second discontinuity
c     yul2 = initial downstream velocity at second discontinuity
c     ydr2 = initial upstream   density  at second discontinuity
c     ypr2 = initial upstream   pressure at second discontinuity
c     yur2 = initial upstream   velocity at second discontinuity
c-----------------------------------------------------------------------
c
c     density
c
      ydl1 = 1.0D0
      ydr1 = 1.0D0
      ydl2 = 1.0D0
      ydr2 = 1.0D0
c
c     pressure
c
      ypl1 = 1000.00D0
      ypr1 =    0.01D0
      ypl2 =    0.01D0
      ypr2 =  100.00D0
c
c     velocity
c
      yul1 = 0.0D0
      yur1 = 0.0D0
      yul2 = 0.0D0
      yur2 = 0.0D0
c
c-----------------------------------------------------------------------
c     geometry: lgeom = 0  => planar geometry
c-----------------------------------------------------------------------
c
      lgeom = 0
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
       g    = 0.0D0                 ! geometry is planar; gravity is off
      cgrav = 0.0D0                 ! that's ok too
c     xmext = 0.0D0                 ! set in exterior mass definition
c
c-----------------------------------------------------------------------
c     radiation transport: lrad  = 0  => no radiation
c-----------------------------------------------------------------------
c
      lrad = 0
      teff = 0.0D0
c
c-----------------------------------------------------------------------
c     turbulent transport: ltur  = 0  => no turbulence
c-----------------------------------------------------------------------
c
      ltur = 0
c
c-----------------------------------------------------------------------
c     declare equations to be solved (hydro only for this problem)
c     ir > 0 => grid         eqn 
c     im > 0 => mass defin.  eqn;  id > 0 => continuity eqn
c     iu > 0 => gas momentum eqn;  it > 0 => gas energy eqn 
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
      jr = ir
      jm = im
      jd = id
      ju = iu
      jt = it
c
      neqn  = 5
      ngrs  = 4
c     ncell should be a multiple of 10 for this problem. note: stone and
c     norman uses 1200 cells; woodward and collela use up to 2400.
      ncell = 200
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
      thet = 0.55D0 
c
      ql0   = 1.0D-3
      ql1   = 0.0D0
      cq2   = 4.0D0
      cq1   = 0.0D0
c*****cq1   = 0.1D0 * cq2
      cqvis = 4.0D0 / 3.0D0
      q0    = 1.0D-30
c
      cadv   = 2.0D0
      epsadv = 1.0D-50
c
      sigd = 1.0D-3                              ! 1.D-3 smears too much
      sige = 1.0D-3                              ! but evolves very fast
c
c-----------------------------------------------------------------------
c     hydrodynamic boundary conditions
c     leibc = 1  => eulerian inner bdy with zero hydrodynamic flux
c     leobc = 1  => eulerian outer bdy with zero hydrodynamic flux
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
      phil0 = 0.0D0
      phil1 = 0.0D0
      phil2 = 0.0D0
      phir0 = 0.0D0
      phir1 = 0.0D0
      phir2 = 0.0D0
c
c-----------------------------------------------------------------------
c     grid specification: lgrid = 1  => adaptive grid
c                         lgrid = 2  => eulerian grid
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
c     choose a scale factor  "xscale" .   Set  "lady(m)"  to the desired
c     option for variables actually used in defining the grid;   for all
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
	    tau   = 1.0D-6
	    ql0   = 1.0D-4
c
	    ladx     = 1
	    lady( 1) = 0
	    lady( 2) = 2
	    lady( 3) = 0
	    lady( 4) = 0
	    lady( 5) = 2 
	    lady( 6) = 2
	    lady( 7) = 0
	    lady( 8) = 0
	    lady( 9) = 0
	    lady(10) = 0
	    lady(11) = 0
c
	    if (ladx     .eq. 1) xscale   = 1.0D0
	    if (lady( 1) .gt. 0) yscl( 1) = xindef
	    if (lady( 2) .eq. 1) yscl( 2) = 1.0D0 
	    if (lady( 3) .eq. 1) yscl( 3) = xindef 
	    if (lady( 4) .eq. 1) yscl( 4) = xindef 
	    if (lady( 5) .eq. 1) yscl( 5) = xindef
	    if (lady( 6) .eq. 1) yscl( 6) = 1.0D0
	    if (lady( 7) .eq. 1) yscl( 7) = xindef 
	    if (lady( 8) .eq. 1) yscl( 8) = 1.0D3 
	    if (lady( 9) .eq. 1) yscl( 9) = xindef
	    if (lady(10) .eq. 1) yscl(10) = xindef
	    if (lady(11) .eq. 1) yscl(10) = xindef
c
	    do 66  l = 1, mad
	    lsum = lsum + lady(l)
   66       continue
c
	    wt( 1) = 0.0D0
	    wt( 2) = 1.0D0
	    wt( 3) = 0.0D0
	    wt( 4) = 0.0D0
	    wt( 5) = 1.0D0
	    wt( 6) = 1.0D0
	    wt( 7) = 0.0D0
	    wt( 8) = 0.0D0
	    wt( 9) = 0.0D0
	    wt(10) = 0.0D0
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
      if (ncell .lt. 2                        ) stop 'start03'
      if (mod(ncell, 2) .ne. 0                ) stop 'start04'
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
      nout  =  10
c
c-----------------------------------------------------------------------
c     document all the parameter settings
c-----------------------------------------------------------------------
c
      if (idoc .gt. 0) call pardoc
c
c-----------------------------------------------------------------------
c     eos: leos = 2  => analytical formula
c     set gas gamma and mean molecular weight
c-----------------------------------------------------------------------
c
      leos = 2
      gam = 1.4D0
      gmu = 1.0D0
c
c-----------------------------------------------------------------------
c     temperature
c-----------------------------------------------------------------------
c
      ytl1 = (cm0 / ck) * gmu * ypl1 / ydl1
      ytr1 = (cm0 / ck) * gmu * ypr1 / ydr1
      ytl2 = (cm0 / ck) * gmu * ypl2 / ydl2
      ytr2 = (cm0 / ck) * gmu * ypr2 / ydr2
c
c=======================================================================
c     generate initial solution
c=======================================================================
c
c++++++++ linit .eq. 1 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (linit .eq. 1) then
c
c-----------------------------------------------------------------------
c     set initial spatial run of variables
c-----------------------------------------------------------------------
c.......................................................................
c     radial grid
c.......................................................................
c
      rn(ngrs  ) = 1.0D0
      rn(ngre+1) = 2.0D0
c
      if (ladx .eq. 1) xscale  = rn(ngre+1) - rn(ngrs)
c
      dr = (rn(ngre+1) - rn(ngrs)) / float(ncell)
c
      do 1 k = ngrs + 1, ngre
      rn(k) = rn(ngrs) + dr * float(k - ngrs)
    1 continue
c
      do 2 k = ngrs, ngre + 1
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
    2 continue
c
      do 3 k = ngrs, ngre
      dvoln(k) = rn(k+1) - rn(k)
    3 continue
c
c.......................................................................
c     density and temperature
c.......................................................................
c
      x1 = rn(ngrs) + 0.1D0 * (rn(ngre+1) - rn(ngrs))
      x2 = rn(ngrs) + 0.9D0 * (rn(ngre+1) - rn(ngrs))
      xw = (rn(ngre+1) - rn(ngrs)) / (1.0D3 * float(ncell))
c
      ad1 =   0.50D0 * (ydr1 + ydl1)
      bd1 =   0.75D0 * (ydr1 - ydl1)
      cd1 =   0.00D0
      dd1 = - 0.25D0 * (ydr1 - ydl1)
c
      at1 =   0.50D0 * (ytr1 + ytl1)
      bt1 =   0.75D0 * (ytr1 - ytl1)
      ct1 =   0.00D0
      dt1 = - 0.25D0 * (ytr1 - ytl1)
c
      ad2 =   0.50D0 * (ydr2 + ydl2)
      bd2 =   0.75D0 * (ydr2 - ydl2)
      cd2 =   0.00D0
      dd2 = - 0.25D0 * (ydr2 - ydl2)
c
      at2 =   0.50D0 * (ytr2 + ytl2)
      bt2 =   0.75D0 * (ytr2 - ytl2)
      ct2 =   0.00D0
      dt2 = - 0.25D0 * (ytr2 - ytl2)
c
      do 4 k = ngrs, ngre
c
      rmid = 0.5D0 * (rn(k) + rn(k+1))
      z1 = (rmid - x1) / xw
      z2 = (rmid - x2) / xw
c
      dn (k) = 1.0D0
      tn (k) = 0.0D0
c
      tn (k) = cvmgt(ytl1, tn(k), z1 .lt. - 1.0D0)
      tn (k) = cvmgt(at1 + z1 * (bt1 + z1 * (ct1 + z1 * dt1)), tn(k), 
     .                            z1 .ge. -1.0D0 .and. z1 .le. +1.0D0)
      tn (k) = cvmgt(ytr1, tn(k), z1 .gt. +1.0D0 .and. z2 .lt. -1.0D0)
      tn (k) = cvmgt(at2 + z2 * (bt2 + z2 * (ct2 + z2 * dt2)), tn(k), 
     .                            z2 .ge. -1.0D0 .and. z2 .le. +1.0D0)
      tn (k) = cvmgt(ytr2, tn(k), z2 .gt. +1.0D0)
c
    4 continue
c
c.......................................................................
c     mass and velocity
c.......................................................................
c
      xmn(ngrs) = 1.0D0 
      un (ngrs) = 0.0D0
c
      do 5 k = ngrs + 1, ngre + 1
      xmn(k) = xmn(k-1) + dn(k-1) * dvoln(k-1)
      un (k) = 0.0D0
    5 continue
c
c.......................................................................
c     phantom zones
c.......................................................................
c
      k = ngrs - 1
      rn (k  ) =   rn (k+1) * 2.0D0 - rn (k+2)
      dn (k  ) =   dn (k+1)
      xmn(k  ) =   xmn(k+1) - dn(k  )/xmup1*(rn(k+1)**mup1-rn(k)**mup1)
      un (k  ) = - un (k+2)
      tn (k  ) =   tn (k+1)
c
      k = ngre + 1
      rn (k+1) =   rn (k  ) * 2.0D0 - rn (k-1)
      dn (k  ) =   dn (k-1)
      dn (k+1) =   dn (k  )
      xmn(k+1) =   xmn(k  ) + dn(k+1)/xmup1*(rn(k+1)**mup1-rn(k)**mup1)
      un (k+1) = - un (k-1)
      tn (k  ) =   tn (k-1)
      tn (k+1) =   tn (k  )
c
c.......................................................................
c     eos
c.......................................................................
c
      call eos
c
c.......................................................................
c     initiate old time solution
c.......................................................................
c
      do 6 k = ngrs - 1, ngre + 2
      ro  (k) = rn (k)
      to  (k) = tn (k)
      do  (k) = dn (k)
      xnto(k) = 0.0D0
    6 continue
c
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 34 k = ngrs, ngre + 1
      r0    (k) = rn (k)
      xm0   (k) = xmn(k)
      u0    (k) = un (k)
      tescr (k) = dn (k) * dvoln(k) * egn(k)
      tescr0(k) = 0.0D0
   34 continue
c
      tee   = 0.0D0
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
c
c.......................................................................
c     zero unneeded variables
c.......................................................................
c
      do 8 k = ngrs, ngre
      dsn   (k) = 0.0D0
      egsn  (k) = 0.0D0
      etn   (k) = 0.0D0
      etsn  (k) = 0.0D0
      ern   (k) = 0.0D0
      ersn  (k) = 0.0D0
      egtn  (k) = 0.0D0
      egtsn (k) = 0.0D0
      egrn  (k) = 0.0D0
      egrsn (k) = 0.0D0
      egrtn (k) = 0.0D0
      egrtsn(k) = 0.0D0
      frn   (k) = 0.0D0
      frsn  (k) = 0.0D0
      chifn (k) = 0.0D0
      avchin(k) = 0.0D0
      xken  (k) = 0.0D0
      xkpn  (k) = 0.0D0
      xnen  (k) = 0.0D0
      plfn  (k) = 0.0D0
      fedd  (k) = 0.0D0
      frnom (k) = 0.0D0
       unom (k) = 1.0D0
    8 continue
c
      do 9 k = ngrs, ngre + 1
      r0  (k) = rn (k)
      xm0 (k) = xmn(k)
      usn (k) = 0.0D0
      urel(k) = 0.0D0
      frn (k) = 0.0D0
      dmdt(k) = 0.0D0
      drdt(k) = 0.0D0
      qx  (k) = 0.0D0
      xnu (k) = 0.0D0
    9 continue
c
      geddl = 0.0D0
      geddr = 0.0D0
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
c     initialize time and timestep
c.......................................................................
c
      timen = 0.0D0
      timeo = 0.0D0
      dtime = 1.0D-4
      tfac  = 1.0D0
c
c.......................................................................
c     initialize grid equation
c.......................................................................
c
      if (lgrid .eq. 1) call grid4nit (ydl1,ydl2,ytl1,ytl2,
     .                                 ydr1,ydr2,ytr1,ytr2,x1,x2,xw)
c
c-----------------------------------------------------------------------
c     zero unneeded variables
c-----------------------------------------------------------------------
c
      do 10 k = ngrs, ngre
      dsn   (k) = 0.0D0
      egsn  (k) = 0.0D0
      etn   (k) = 0.0D0
      etsn  (k) = 0.0D0
      ern   (k) = 0.0D0
      ersn  (k) = 0.0D0
      egtn  (k) = 0.0D0
      egtsn (k) = 0.0D0
      egrn  (k) = 0.0D0
      egrsn (k) = 0.0D0
      egrtn (k) = 0.0D0
      egrtsn(k) = 0.0D0
      frn   (k) = 0.0D0
      frsn  (k) = 0.0D0
      chifn (k) = 0.0D0
      avchin(k) = 0.0D0
      xken  (k) = 0.0D0
      xkpn  (k) = 0.0D0
      xnen  (k) = 0.0D0
      plfn  (k) = 0.0D0
      fedd  (k) = 0.0D0
      frnom (k) = 0.0D0
       unom (k) = 1.0D0
   10 continue
c
      do 11 k = ngrs, ngre + 1
      r0  (k) = rn (k)
      xm0 (k) = xmn(k)
      usn (k) = 0.0D0
      urel(k) = 0.0D0
      frn (k) = 0.0D0
      dmdt(k) = 0.0D0
      drdt(k) = 0.0D0
   11 continue
c
      geddl = 0.0D0
      geddr = 0.0D0
c
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 12 k = ngrs, ngre + 1
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
   12 continue
c
      do 13 k = ngrs, ngre
      dvoln (k) = ( rmup1n(k+1) - rmup1n(k) ) / xmup1
      tescr (k) = dn(k) * dvoln(k) * egn(k) 
      tescr0(k) = 0.0D0
   13 continue
c
      tee   = ssum(ngre - ngrs + 1, tescr(ngrs), 1)
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
c
c-----------------------------------------------------------------------
c     mass and velocity
c-----------------------------------------------------------------------
c
      xmn(ngrs) = 1.0D0 
      un (ngrs) = 0.0D0
c
      do 14 k = ngrs + 1, ngre + 1
      xmn(k) = xmn(k-1) + dn(k-1) * dvoln(k-1)
      un (k) = 0.0D0
   14 continue
c
c.......................................................................
c     exterior mass
c.......................................................................
c
      xmass = xmn(ngre+2)
      xmext = 1.0D-10
      xmtot = xmass + xmext
c
      do 155 k = ngrs - 1, ngre + 2
      xmen(k) = xmtot - xmn(k)
  155 continue
c
c-----------------------------------------------------------------------
c     initialize time and timestep
c-----------------------------------------------------------------------
c
      timen = 0.0D0
      timeo = 0.0D0
      dtime = 1.0D-6
      tfac  = 1.0D0
c
c=======================================================================
c     write initial model on dump file
c=======================================================================
c
      if (idump .gt. 0) then
                               irec = 1
          write ( idump, rec = irec, iostat = ier)
     .           header, ydl1, ypl1, yul1, ytl1, ydr1, ypr1, yur1, ytr1,
     .                   ydl2, ypl2, yul2, ytl2, ydr2, ypr2, yur2, ytr2
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
c++++++++ linit .eq. 1 +++++++++++++++++++++++++++++++++++++++++++++ end
c
      return
      end if
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
     .                              yheader, yydl1, yypl1, yyul1, yytl1,
     .                                       yydr1, yypr1, yyur1, yytr1,
     .                                       yydl2, yypl2, yyul2, yytl2,
     .                                       yydr2, yypr2, yyur2, yytr2 
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
c++++++++ linit .ne. 1 +++++++++++++++++++++++++++++++++++++++++++++ end
c
c=======================================================================
c
   20 format(/' model jdump ='i6' jstep ='i6' archived')
   21 format(/' unit idump. irec = 'i3' ier ='i20, o22, a8)
   22 format( ' jdumpo ='i7' jdumpi ='i7' jstepo ='i7' jstepi ='i7)
c
      end
