      subroutine start01
c
c***********************************************************************
c     startup for YOUR first FAVOURITE PROBLEM
c.......................................................................
c     calling sequence: titan > start > start01
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
c
      read (iin ,'(a80)') header
c
      write(iout,'(a80)') header
      write(itty,'(a80)') header
c
c-----------------------------------------------------------------------
c     geometry: lgeom = 0  => planar    geometry
c                     = 2  => spherical geometry
c-----------------------------------------------------------------------
c
      lgeom = 1
       mu   = lgeom
       mum1 = mu - 1
       mup1 = mu + 1
      xmu   = float(mu  )
      xmum1 = float(mum1)
      xmup1 = float(mup1)
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
       g    = 0.0D0
      cgrav = 0.0D0
c     xmext = 0.0D0                               ! set in exterior mass
c
c-----------------------------------------------------------------------
c     radiation transport: lrad  = 0  => without radiation
c                                = 1  => with    radiation
c                          ltran = 1  => full transport
c                                = 2  => non-equilibrium diffusion
c                                = 3  =>     equilibrium diffusion
c                          lam   = 0  => without flux limiter
c                                = 1  => with    flux limiter
c                          ncor  = number of rays intersecting core
c-----------------------------------------------------------------------
c
      lrad  = 0
      ltran = 0
      lam   = 0
      ncor = 20
c
c-----------------------------------------------------------------------
c     radiation boundary conditions 
c     lribc = 1  => optically transmitting (inner/left) boundary
c           = 2  => optically reflecting          left  boundary 
c                   or Milne                inner       boundary
c           = 3  => imposed net flux at    (inner/left) boundary
c     lrobc = 1  => optically transmitting (outer/rite) boundary
c           = 2  => optically reflecting          rite  boundary 
c           = 3  => imposed net flux at    (outer/rite) boundary
c-----------------------------------------------------------------------
c
      lribc = 0
      lrobc = 0
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
c-----------------------------------------------------------------------
c     ir > 0 => grid         eqn 
c     im > 0 => mass defin.  eqn;  id > 0 => continuity eqn
c     iu > 0 => gas momentum eqn;  it > 0 => gas energy eqn 
c     if > 0 => rad momentum eqn;  ie > 0 => rad energy eqn
c
c     neqn  = number of equations at each grid-point
c     ngrs  = index of first (inner/left) cell
c     ngre  = index of last  (outer/rite) cell
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
c
      neqn = 0
c
      ngrs  = 5
      ncell = 100
      ngre  = ngrs + ncell - 1
c
c-----------------------------------------------------------------------
c     numerics 
c
c     lcray = 0  => nonCRAY computer; use linpack routines
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
c     set hydro parameters: lhydr = 0  => without hydrodynamics
c                                 = 1  => with    hydrodynamics
c     thet   = time centering parameter
c     ql0    = fixed absolute length for pseudoviscosity (planar)
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
c
      thet  = 0.55D0 
c
      ql0   = 0.0D0
      ql1   = 1.0D-3
      cq2   = 4.0D0
      cq1   = 0.0D0
      cqvis = 4.0D0 / 3.0D0
      q0    = 1.0D-30
c
        cadv = 2.0D0
      epsadv = 1.0D-50
c
      sigd = 0.0D0
      sige = 0.0D0
c
c-----------------------------------------------------------------------
c     hydrodynamic boundary conditions
c-----------------------------------------------------------------------
c     leibc = 1  => Eulerian     zero flux (inner/left) boundary
c           = 2  => Eulerian non-zero flux (inner/left) boundary
c     leobc = 1  => Eulerian     zero flux (outer/rite) boundary
c           = 2  => Eulerian non-zero flux (outer/rite) boundary
c           = 3  => Eulerian  transmitting (outer/rite) boundary
c
c     llibc = 1  => Lagrangean             (inner/left) boundary 
c			       with specified velocity
c     llobc = 1  => Lagrangean             (outer/rite) boundary 
c			       with specified velocity
c           = 2  => Lagrangean             (outer/rite) boundary
c			       with specified external pressure
c           = 3  => Lagrangean             (outer/rite) boundary
c			       transmitting
c.......................................................................
c     phil0 = mass     flux (cgs) through inner (left ) eulerian bdy
c     phir0 = mass     flux (cgs) through outer (right) eulerian bdy
c     phil1 = momentum flux (cgs) through inner (left ) eulerian bdy
c     phir1 = momentum flux (cgs) through outer (right) eulerian bdy
c     phil2 = energy   flux (cgs) through inner (left ) eulerian bdy
c     phir2 = energy   flux (cgs) through outer (right) eulerian bdy
c-----------------------------------------------------------------------
c
      llibc = 0
      llobc = 0
      leibc = 0
      leobc = 0
c
      phil0 = 0.0D0
      phil1 = 0.0D0
      phil2 = 0.0D0
      phir0 = 0.0D0
      phir1 = 0.0D0
      phir2 = 0.0D0
c
      tl    = 0.0D0
      dl    = 0.0D0
      ul    = 0.0D0
      tr    = 0.0D0
      dr    = 0.0D0
      ur    = 0.0D0
      delrl = 0.0D0
      delrr = 0.0D0
      delml = 0.0D0
      delmr = 0.0D0
      uextl = 0.0D0
      uextr = 0.0D0
      pextr = 0.0D0
      omega = 0.0D0
c
c-----------------------------------------------------------------------
c     grid specification: lgrid = 1  => adaptive   grid
c                               = 2  => Eulerian   grid
c                               = 3  => Lagrangean grid
c.......................................................................
c     for an adaptive grid:
c
c     ladx    = 2  => logarithmic resolution
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
c     choose a scale factor  "xscale" .  Set  "lady(m)"  to the desired
c     option for variables actually  used in defining the grid;  for all
c     other variables set lady(m) to zero. If linear resolution has been
c     chosen for a variable (lady(m) .eq. 1),   the scale factor yscl(m)
c     MUST  be set.  Choose nonzero weights  "wt(m)"  only for variables
c     used to define grid; take unity as default.Set all others to zero.
c     "yscl(1)" MUST be specified when "lady(1)" is turned on.
c.......................................................................
c     Values recommended by dorfi: alph = 2, bet = 1, tau = 1.e-7 * tff
c     for cepheids, and tau = 1.e-1 * tff for AGB stars, where tff is
c     the free-fall time of the core.  stein suggests bet = 2
c-----------------------------------------------------------------------
c
            lsum = 1
           lgrid = 0
c
      if ( lgrid .eq. 1 ) then
            alph  = 0.0D0
            ibet  = 0
            tau   = 0.0D0
c
	    ladx     = 0
	    lady( 1) = 0     ! need to specify mass scale xmscl for grid
	    lady( 2) = 0
	    lady( 3) = 0
	    lady( 4) = 0
	    lady( 5) = 0
	    lady( 6) = 0
	    lady( 7) = 0
	    lady( 8) = 0
	    lady( 9) = 0
	    lady(10) = 0
	    lady(11) = 0

            if (ladx     .eq. 1) xscale   = xindef 
	    if (lady( 1) .gt. 0) yscl( 1) = xindef  ! e.g. xmass,xmn(ns)
	    if (lady( 2) .eq. 1) yscl( 2) = xindef
	    if (lady( 3) .eq. 1) yscl( 3) = xindef
	    if (lady( 4) .eq. 1) yscl( 4) = xindef
	    if (lady( 5) .eq. 1) yscl( 5) = xindef
	    if (lady( 6) .eq. 1) yscl( 6) = xindef
	    if (lady( 7) .eq. 1) yscl( 7) = xindef 
	    if (lady( 8) .eq. 1) yscl( 8) = xindef
	    if (lady( 9) .eq. 1) yscl( 9) = xindef
	    if (lady(10) .eq. 1) yscl(10) = xindef
	    if (lady(11) .eq. 1) yscl(11) = xindef
c
	    do 66  l = 1, mad
	    lsum = lsum + lady(l)
   66       continue
c
	    wt( 1) = 0.0D0
	    wt( 2) = 0.0D0
	    wt( 3) = 0.0D0
	    wt( 4) = 0.0D0
	    wt( 5) = 0.0D0
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
      if ( lgeom + lcray + lband + lgrid    .le. 0 ) stop'start01'
      if ( lgrid .eq. 1 .and. lsum          .le. 0 ) stop'start02'
      if ( lgrid .eq. 1 .and. leibc + leobc .le. 0 ) stop'start03'
      if ( lgrid .eq. 3 .and. llibc * llobc .le. 0 ) stop'start04'
      if ( linit .ne. 1 .and. idump         .le. 0 ) stop'start05'
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
c.......................................................................
c     jdump = number of last dump file created
c     jstep = number of timestep  at dump number eq jdump
c     nstep = number of timesteps in current run
c=======================================================================
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
c-----------------------------------------------------------------------
c     ndump = number of timesteps between dumps     of models
c     nout  = number of timesteps between printouts of models
c=======================================================================
c
      ndump = 0
      nout  = 0
c
c=======================================================================
c     document all the parameter settings
c=======================================================================
c
      if ( idoc .gt. 0 ) call pardoc
c
c-----------------------------------------------------------------------
c     eos: leos = 1  => tables
c               = 2  => perfect gas
c               = 3  => STELLINGWERF formula
c-----------------------------------------------------------------------
c
      leos = 0
      gam   = 5.0D0 / 3.0D0
      gmu   = 0.50D0
      xabun = 0.700D0
      yabun = 0.299D0
      zabun = 0.001D0
c
c-----------------------------------------------------------------------
c     opacity: lopac = 1  => tables
c                    = 2  => STELLINGWERF opacity fit
c                    = 3  => constant opacity 
c                    = 4  => THOMSON free electron opacity
c-----------------------------------------------------------------------
c
      lopac = 0
      chif0 = 1.0D0
      ratio = 1.0D0
c
c=======================================================================
c     generate initial model
c=======================================================================
c
c+++++++++++++++++++++++++++++++ linit .eq. 1 ++++++++++++++++++++ begin
      if ( linit .eq. 1 ) then
c
c=======================================================================
c     set useful quantitites
c=======================================================================
c-----------------------------------------------------------------------
c     zero unneeded variables
c-----------------------------------------------------------------------
c
      do 16 k = ngrs - 1, ngre + 1
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
      tescr0(k) = 0.0D0
   16 continue
c
      do 17 k = ngrs - 1, ngre + 1
      r0   (k) = rn (k)
      xm0  (k) = xmn(k)
      usn  (k) = 0.0D0
      urel (k) = 0.0D0
      dmdt (k) = 0.0D0
      drdt (k) = 0.0D0
      frnom(k) = 1.0D0
       unom(k) = 1.0D0
   17 continue
c
c-----------------------------------------------------------------------
c     radial grid
c-----------------------------------------------------------------------
c
      do 18 k = ngrs - 1, ngre + 1
      rmun  (k) = rn(k)**mu
      rmup1n(k) = rn(k)**mup1
      rmum1n(k) = rn(k)**mum1
   18 continue
c
      do 19 k = ngrs - 1, ngre
      dvoln(k) = ( rmup1n(k+1) - rmup1n(k) ) / xmup1
   19 continue
c
c-----------------------------------------------------------------------
c     mass and velocity
c-----------------------------------------------------------------------
c
      do 30 k = ngre, ngrs, -1
      xmn  (k) = xmn(k+1) - dn(k) * dvoln(k) 
      un   (k) = 0.0D0
   30 continue
c
      delrl = 0.0D0
      delrr = 0.0D0
      delml = 0.0D0
      delmr = 0.0D0
c
c-----------------------------------------------------------------------
c     radiation field
c-----------------------------------------------------------------------
c
      do 31 k = ngrs, ngre + 1
      ern  (k) =  car * tn(k)**4
      plfn (k) =  cc  / (4.0D0 * cpi) * ern(k)
      dplfdltn(k) =      4.0D0 * plfn(k)
      frn  (k) = xlum / (4.0D0 * cpi * rmun(k) ) 
   31 continue
c
      call opac
c
      call eddfac
c
      ern(ngre  ) = frn(ngre) / (cc * geddr)
      ern(ngre+1) = frn(ngre) / (cc * geddr)
      ern(ngre+2) = frn(ngre) / (cc * geddr)
c
      do 32 k = ngre - 1, ngrs, -1
      ern(k) = (chifn(k+1) * dn(k+1) * frn(k+1) / cc
     .          +  2.0D0 * fedd(k+1) * ern(k+1) / (rn(k+2) - rn(k))
     .  + 0.5D0 * (3.0D0 * fedd(k+1) - 1.0D0) * ern(k+1) / rn(k+1) )
     .          / (2.0D0 * fedd(k  )            / (rn(k+2) - rn(k))
     .  - 0.5D0 * (3.0D0 * fedd(k  ) - 1.0D0)            / rn(k+1) )
   32 continue
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
      xmn(k  ) =  xmn(k+1) * 2.0D0 - xmn(k+2)
c
      k = ngre + 1
      tn (k  ) =  tn (k-1)
      tn (k+1) =  tn (k  )
      dn (k  ) =  dn (k-1)
      dn (k+1) =  dn (k  )
      rn (k+1) =  rn (k  ) * 2.0D0 - rn (k-1)
      xmn(k+1) =  xmn(k  ) * 2.0D0 - xmn(k-1)
c
      geddl = 0.0D0
      geddr = 0.0D0
c
c-----------------------------------------------------------------------
c     exterior mass
c-----------------------------------------------------------------------
c
      xmass = xmn(ngre+2)
      xmext = 0.0D0
      xmtot = xmass + xmext
c
      do 155 k = ngrs - 1, ngre + 2
      xmen(k) = xmtot - xmn(k)
  155 continue
c
c-----------------------------------------------------------------------
c     eos
c-----------------------------------------------------------------------
c
      call eos
c
      do 33 k = ngrs - 1, ngre + 2
      egrn (k) = ern(k) / dn(k) + egn(k)
   33 continue
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
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 34 k = ngrs, ngre
      tescr(k) = dn(k) * dvoln(k) * egrn(k) 
   34 continue
c
      tee   = ssum(ngre - ngrs + 1, tescr(ngrs), 1)
      tes   = 0.0D0
      tew   = 0.0D0
      tel   = 0.0D0
      teq   = 0.0D0
      tetot = tee
c
c-----------------------------------------------------------------------
c     initialize old time solution
c-----------------------------------------------------------------------
c
      do 35 k = ngrs - 1, ngre + 2
      ro  (k) = rn(k)
      to  (k) = tn(k)
      do  (k) = dn(k)
      xnto(k) = 0.0D0
   35 continue
c
c-----------------------------------------------------------------------
c     initialize time and timestep
c-----------------------------------------------------------------------
c
      timen = 0.0D0
      timeo = 0.0D0
      tfac  = 1.0D0
c
      dtime = 1.0D-15
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
c     write initial model on dump file; ignore if not file specified
c=======================================================================
c
      if (idump .gt. 0 ) then
                               irec = 1
          write ( idump, rec = irec, iostat = ier)
     .                                   header, xmass, xlum, xrad, teff
c
          if (ier .gt. 0) then
              write (itty, 21) irec, ier, ier, ier
              write (iout, 21) irec, ier, ier, ier
              stop'start15' 
          end if
      end if
c
      read (iin,'(3i5,e10.3)') jdump, jstep, nstep, tmax
c
      jsteps = jstep + 1
      jstepe = jsteps + nstep - 1
      if ( jdump .ne. 1 .or. jstep .ne. 1 ) stop 'start16'
c
      write (itty, 20) jdump, jstep
      write (iout, 20) jdump, jstep
c
      if (idump .gt. 0 ) then
                               irec = 2
          write ( idump, rec = irec, iostat = ier) jdump, jsteps, jstepe
c
          if (ier .gt. 0) then
              write (itty, 21) irec, ier, ier, ier
              write (iout, 21) irec, ier, ier, ier
              stop 'start17' 
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
              stop 'start18' 
          end if
      end if
c
      call writout
c
c++++++++++++++++++++++++++++++ linit .eq. 1 ++++++++++++++++++++++++end
c
      return
      end if
c
c=======================================================================
c     recover model from dump file; error if no file specified
c=======================================================================
c
c+++++++++++++++++++++++++++++++ linit .ne. 1 +++++++++++++++++++++begin
      if ( linit .ne. 1 ) then
c
c-----------------------------------------------------------------------
c     compare model header
c-----------------------------------------------------------------------
c
      if (idump .gt. 0 ) then
                               irec = 1
          read  ( idump, rec = irec, iostat = ier)
     .                                 yheader, ymass, ylum, yrad, yteff
c
          if (ier .gt. 0) then
              write (itty, 21) irec, ier, ier, ier
              write (iout, 21) irec, ier, ier, ier
              stop 'start20' 
          end if
c
          if (yheader .ne. header .or. ymass .ne. xmass .or.
     .        ylum    .ne. xlum   .or. yrad  .ne. xrad  .or.
     .        yteff   .ne. teff      ) then
c
              write(itty, '(a80/a80/1p4e15.7/1p4e15.7)') 
     .                                 header, yheader,
     .                                 xmass , xlum, xrad,  teff, 
     .                                 ymass , ylum, yrad, yteff
              stop 'start21' 
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
      if (idump .gt. 0 ) then
                               irec = 2
          read  ( idump, rec = irec, iostat = ier) jdumpo, jstepo
c
          if (ier .gt. 0) then
              write (itty, 21) irec, ier, ier, ier
              write (iout, 21) irec, ier, ier, ier
              stop 'start22' 
          end if
c
          if (jdumpo .lt. jdumpi .or. jstepo .lt. jstepi) then
              write (itty, 22) jdumpo, jdumpi, jstepo, jstepi
              stop 'start23'
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
              stop 'start24' 
          end if
c
c-----------------------------------------------------------------------
c     check that we got the right one
c-----------------------------------------------------------------------
c
          if (jdumpd .ne. jdumpi .or. jstepd .ne. jstepi) then
              write (itty, 22) jdumpd, jdumpi, jstepd, jstepi
              stop 'start25'
          end if
      end if
c
      return
      end if
c
c++++++++++++++++++++++++++++++ linit .ne. 1 ++++++++++++++++++++++++end
c
c=======================================================================
c
   20 format(/' model jdump ='i6' jstep ='i6' archived')
   21 format(/' unit idump. irec = 'i3' ier ='i20, o22, a8)
   22 format( ' jdumpo ='i7' jdumpi ='i7' jstepo ='i7' jstepi ='i7)
c
      end
