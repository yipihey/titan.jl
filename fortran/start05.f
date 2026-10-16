      subroutine start05
c
c***********************************************************************
c     setup for TAYLOR-SEDOV BLAST WAVE problem
c.......................................................................
c     calling sequence: start > start05 > pardoc, eos, writout, gridinit
c***********************************************************************
c     comparison results  exist in the  literature for times t = 1.08e3,
c     1.16e4, 4.4e4, 8.45e4, 1.46e5, and 2.88e5 seconds.
c
c     AA, 266, 266-282 (1992): 302 grid points; alph = 2; tau = 0
c                              solution at time t=0.02433 is given.
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
c     rho0 = initial density (gm/cm**3)
c     t0   = initial temperature (K)
c     rmax = initial outer radius (cm)
c-----------------------------------------------------------------------
c
      read (iin,'(a80)') header
      read (iin,'(3e10.0)') rho0, t0, rmax
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
      ngrs  = 5
      ncell = 100                            ! ENSMAN's computation: 260
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
      thet  = 0.55D0 
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
c     leibc = 1  => zero flux eulerian inner bdy
c     leobc = 1  => zero flux eulerian outer bdy
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
      uextl = 0.0D0
      pextr = 0.0D0
      omega = 1.0D30
c
c-----------------------------------------------------------------------
c     grid specification: lgrid = 1  => adaptive   grid
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
	    tau   = 1.0D-15
            ql1   = 1.0D-4
c
	    ladx     = 2
	    lady( 1) = 0
	    lady( 2) = 2
	    lady( 3) = 0
	    lady( 4) = 0
	    lady( 5) = 2
	    lady( 6) = 0
	    lady( 7) = 0
	    lady( 8) = 0
	    lady( 9) = 0
	    lady(10) = 0
	    lady(11) = 0
c
	    if (ladx     .eq. 1) xscale   = rmax
	    if (lady( 1) .gt. 0) yscl( 1) = xindef
	    if (lady( 2) .eq. 1) yscl( 2) = rho0
	    if (lady( 3) .eq. 1) yscl( 3) = xindef
	    if (lady( 4) .eq. 1) yscl( 4) = xindef 
	    if (lady( 5) .eq. 1) yscl( 5) = 1.0D0
	    if (lady( 6) .eq. 1) yscl( 6) = 1.0D0
	    if (lady( 7) .eq. 1) yscl( 7) = xindef 
	    if (lady( 8) .eq. 1) yscl( 8) = 1.0D8 
	    if (lady( 9) .eq. 1) yscl( 9) = xindef
	    if (lady(10) .eq. 1) yscl(10) = xindef
	    if (lady(11) .eq. 1) yscl(11) = xindef
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
	    wt( 6) = 0.0D0
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
      niter = 10
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
      nout  = 005
c
c-----------------------------------------------------------------------
c     eos: leos = 2  => analytical formula
c     set gas gamma and mean molecular weight
c-----------------------------------------------------------------------
c
      leos = 2
      gam = 5.0D0 / 3.0D0
      gmu = 0.5D0
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
c     radial grid
c.......................................................................
c
      rn(ngrs  ) = 1.0D0                       ! ENSMAN's inner boundary
      rn(ngrs  ) = rmax * 1.0D-3
      rn(ngrs-1) = rn(ngrs) * 0.999D0
      rn(ngre+1) = rmax
      rn(ngre+2) = rmax     * 1.001D0
c
      dr = (rn(ngre+1) - rn(ngrs)) / float(ncell)
c
      do 1 k = ngrs + 1, ngre
      rn(k) = rn(ngrs) + dr * float(k-ngrs)
    1 continue
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
c     density and temperature
c.......................................................................
c
      do 4 k = ngrs - 1, ngre + 1
      dn(k) = rho0
      tn(k) = t0
    4 continue
c
c.......................................................................
c     mass and velocity
c.......................................................................
c
      xmn(ngrs-1) = 0.0D0
      un (ngrs  ) = 0.0D0
c
      do 5 k = ngrs, ngre + 2
      xmn (k) = xmn(k-1) + dn(k-1) * dvoln(k-1)
      un  (k) = 0.0D0
    5 continue
c
c.......................................................................
c     dump internal energy in first zone
c     recalculate temperature
c.......................................................................
c
      egn(ngrs) = 1.0D50 / ( dn(ngrs) * dvoln(ngrs) ) 
      tn (ngrs) = (gam - 1.0D0) * (gmu * cm0 / ck) * egn(ngrs)
c
c.......................................................................
c     phantom zones
c.......................................................................
c
      k = ngrs - 1
      un (k  ) =  un (k+1)
      egn(k  ) =  egn(k+1)
c
      k = ngre + 1
      tn (k  ) =  tn (k-1)
      tn (k+1) =  tn (k  )
      dn (k+1) =  dn (k  )
c
      delml =  dn(ngrs) * dvoln(ngrs)
      delmr =  dn(ngre) * dvoln(ngre)
      delmr = xmn(ngre) - xmn(ngre-1)
c
c.......................................................................
c     eos
c.......................................................................
c
      call eos
c
c.......................................................................
c     initialize old time solution
c.......................................................................
c
      do 6 k = ngrs - 1, ngre + 2
      ro  (k) = rn(k)
      to  (k) = tn(k)
      do  (k) = dn(k)
      xnto(k) = 0.0D0
    6 continue
c
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 34  k = ngrs, ngre + 1
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
      xnu (k) = 0.0D0
      qx  (k) = 0.0D0
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
      dtime = 1.0D-2
      tfac  = 1.0D0
c
c.......................................................................
c     initialize grid equation
c.......................................................................
c
      if (lgrid .eq. 1) then
c
                        jsteps = 1
                        jstepe = 1
                        call writout
c
                                rw   = rn(ngrs) * .1D0
                                radc = rn(ngrs) + 5.D0 * rw  
                        rn(ngrs-1)   = rn(ngrs) - .5D0 * rw
                        rn(ngrs+1)   = rn(ngrs) + 1.D0 * rw
                        rn(ngrs+2)   = rn(ngrs) + 2.D0 * rw
c
                        ecnt = 1.0D50 / (rho0 * rn(ngrs)**mup1 ) * xmup1
                        ecnt = 0.5D00 / 1.209288D0 * ecnt
                        ytl  = (gam - 1.0D0) * (gmu * cm0 / ck) * ecnt
                        ytr  = tn(ngre)
c
                        call grid5nit (ytl, ytr, radc, rw)
      end if
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
      r0   (k) = rn (k)
      xm0  (k) = xmn(k)
      usn  (k) = 0.0D0
      urel (k) = 0.0D0
      frn  (k) = 0.0D0
      dmdt (k) = 0.0D0
      drdt (k) = 0.0D0
   11 continue
c
      geddl = 0.0D0
      geddr = 0.0D0
c
c-----------------------------------------------------------------------
c     radial grid
c-----------------------------------------------------------------------
c
      do 12 k = ngrs - 1, ngre + 2
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
   12 continue
c
      do 13 k = ngrs - 1, ngre + 1
      dvoln(k) = ( rmup1n(k+1) - rmup1n(k) ) / xmup1
   13 continue
c
c-----------------------------------------------------------------------
c     mass and velocity
c-----------------------------------------------------------------------
c
      xmn(ngrs-1) = 0.0D0
      un (ngrs-1) = 0.0D0
c
      do 14 k = ngrs, ngre + 2
      xmn  (k) = xmn(k-1) + dn(k-1) * dvoln(k-1)
      un   (k) = 0.0D0
   14 continue
c
c.......................................................................
c     exterior mass
c.......................................................................
c
      xmass = xmn(ngre+2)
      xmext = cmsol / (4.D0 * cpi)
      xmtot = xmass + xmext
c
      do 155 k = ngrs - 1, ngre + 2
      xmen(k) = xmtot - xmn(k)
  155 continue
      yscl(1) = xmtot
c
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 15 k = ngrs, ngre
      tescr (k) = dn(k) * dvoln(k) * egn(k)
      tescr0(k) = 0.0D0
   15 continue
c
      tee   = ssum(ngre - ngrs + 1, tescr(ngrs), 1) 
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
      write(iout, '(" total energy =" 1pe13.6)') tetot
c
c-----------------------------------------------------------------------
c     initialize time and timestep
c-----------------------------------------------------------------------
c
      timeo = 0.0D0
      timen = 0.0D0
      dtime = 1.0D-25
      tfac  = 1.0D0
c
c=======================================================================
c     write initial model on dump file
c=======================================================================
c
      if (idump .gt. 0) then
                               irec = 1
          write ( idump, rec = irec, iostat = ier) 
     .                                            header, rho0, t0, rmax
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
     .                                        yheader, yrho0, yt0, yrmax
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
          yscl(1) = xmtot
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
