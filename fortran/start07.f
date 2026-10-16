      subroutine start07
c
c***********************************************************************
c     setup for COOLING of grey sphere from RADIATIVE EQUILIBRIUM 
c.......................................................................
c     calling sequence: start > start07 > pardoc, eddfac, eos, writout
c***********************************************************************
c     for nonequilibrium diffusion comparison results exist in the lite-
c     rature for times t = 9.8e5, 8.0e6, 3.2e7, 1.4e8, 3.1e8 seconds.
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
c           planar    geometry
c                                 off: set g = 0.0D0
c                                 on : set g = desired value
c           spherical geometry
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
      jb = if
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
      phil0 = 0.0D0
      phil1 = 0.0D0
      phil2 = 0.0D0
      phir0 = 0.0D0
      phir1 = 0.0D0
      phir2 = 0.0D0
c
      uextl = 0.0D0
      uextr = 0.0D0
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
           lsum = 1
          lgrid = 2
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
      if (lgeom .ne. 2                        ) stop 'start01'
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
      nout  = 10
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
      gam  = 5.0D0 / 3.0D0
      gmu  = 0.5D0
c
c-----------------------------------------------------------------------
c     opacity: lopac = 3  => constant opacity law
c-----------------------------------------------------------------------
c
      lopac = 3
c     ratio is read in from titan.in
c
c=======================================================================
c     generate initial solution
c=======================================================================
c
c++++++++ linit .eq. 1 ++++++++++++++++++++++++++++++++++++++++++++begin
      if (linit .eq. 1) then
c
c.......................................................................
c     read in final model from start06
c.......................................................................
c
      if (idump .eq. 0) stop 'start: no access to COOL.dmp'
c
      if (idump .gt. 0) then
                               irec = 1
          read  ( idump, rec = irec, iostat = ier) 
     .                         yheader, xmass, ylum, yrho0, ychif0, ximr
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start03' 
          end if
c
          write(itty, '(a80/a80)') yheader, header
          write(iout, '(a80/a80)') yheader, header
c
          if (ylum   .ne. xlum  .or.
     .        ychif0 .ne. chif0 .or. yrho0 .ne. rho0 ) then
	      write(itty,'(1p6e10.3)') 
     .                            xlum, ylum, rho0, yrho0, chif0, ychif0
    	      write(iout,'(1p6e10.3)') 
     .                            xlum, ylum, rho0, yrho0, chif0, ychif0
	      stop 'start04' 
          end if
      end if
c
c.......................................................................
c     read dump number and timestep number of initial model; check that
c     specified starting model exists
c.......................................................................
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
	      stop 'start05' 
          end if
c
c.......................................................................
c     read model
c.......................................................................
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
	      stop 'start06' 
          end if
c
c.......................................................................
c     check that we got the right one
c.......................................................................
c
          if (jdumpd .ne. jdumpi .or. jstepd .ne. jstepi) then
	      write (itty, 22) jdumpd, jdumpi, jstepd, jstepi
	      stop 'start07'
          end if
      end if
c
c-----------------------------------------------------------------------
c     radiation flux. apply zero luminosity at inner boundary
c     and use full flux everywhere else
c-----------------------------------------------------------------------
c
      frn0 = xlum / (4.0D0 * cpi * rmun(ngre+2))
      ern0 = sqrt(3.0D0) * frn0 / cc
      teff = (ern0 / car)**0.25D0
c
      xlum0 = xlum
      xlum = 1.0D0
c
           k = ngrs - 1
      frn (k) = xlum  / (4.0D0 * cpi * rmun(k))
           k = ngrs
      frn (k) = xlum  / (4.0D0 * cpi * rmun(k))
      do 1 k = ngrs + 1, ngre + 2
      frn (k) = xlum0 / (4.0D0 * cpi * rmun(k))
    1 continue
c
           k = ngrs - 1
      ern (k  ) =   ern(k+1)
      frn (k  ) = ( rn (k+1) / rn(k) )**mu * frn(k+1)
           k = ngre + 1
      ern (k  ) =   ern(k-1)
      ern (k+1) =   ern(k  )
      frn (k+1) = ( rn (k) / rn(k+1) )**mu * frn(k  )
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
c-----------------------------------------------------------------------
c     initialize terms in total energy equation
c-----------------------------------------------------------------------
c
      do 2  k = ngrs, ngre
      egrn  (k) = ern (k) / dn(k) + egn(k)
      tescr (k) = egrn(k) * dn(k) * dvoln(k)
      tescr0(k) = 0.0D0
    2 continue
c
      tee   = ssum(ngre - ngrs + 1, tescr(ngrs), 1)
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
      do 3   k = ngrs - 1, ngre + 2
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
      frnom (k) = frn(k)
       unom (k) = 0.0D0
c
      un  (k) = 0.0D0
      u0  (k) = 0.0D0
      r0  (k) = rn (k)
      xm0 (k) = xmn(k)
      usn (k) = 0.0D0
      urel(k) = 0.0D0
      dmdt(k) = 0.0D0
      drdt(k) = 0.0D0
    3 continue
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
              stop 'start08'
          end if
c
                               irec = 2
          write ( idump, rec = irec, iostat = ier) jdump, jstep
c
          if (ier .gt. 0) then
    	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start09' 
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
	      stop 'start10' 
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
	      stop 'start11' 
          end if
c
          if (yheader .ne. header) then
              write(itty, '(a80/a80)') header, yheader
              write(iout, '(a80/a80)') header, yheader
              stop 'start12' 
          end if
c
          ylum  = xlum
c
          if (ymass  .ne. xmass .or. ylum  .ne. xlum .or.
     .        ychif0 .ne. chif0 .or. yrho0 .ne. rho0 ) then
    	      write(itty,'(1p8e10.3)') 
     .              xmass, ymass, xlum, ylum, rho0, yrho0, chif0, ychif0
	      write(iout,'(1p8e10.3)') 
     .              xmass, ymass, xlum, ylum, rho0, yrho0, chif0, ychif0
	      stop 'start13' 
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
	      stop 'start14' 
          end if
c
          if (jdumpo .lt. jdumpi .or. jstepo .lt. jstepi) then
    	      write (itty, 22) jdumpo, jdumpi, jstepo, jstepi
	      stop 'start15'
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
          frn0 = xlum / (4.0D0 * cpi * rmun(ngre+2))
          ern0 = sqrt(3.0D0) * frn0 / cc
          teff = (ern0 / car)**0.25D0
          xlum = 1.0D0
c
          if (ier .gt. 0) then
	      write (itty, 21) irec, ier, ier, ier
	      write (iout, 21) irec, ier, ier, ier
	      stop 'start16' 
          end if
c
c-----------------------------------------------------------------------
c     check that we got the right one
c-----------------------------------------------------------------------
c
          if (jdumpd .ne. jdumpi .or. jstepd .ne. jstepi) then
	      write (itty, 22) jdumpd, jdumpi, jstepd, jstepi
c             stop 'start15'
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
