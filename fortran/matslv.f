      subroutine matslv
c
c***********************************************************************
c     compute solution of grand system as either a  block pentadiagonal 
c     or a banded system. both standard linpack solvers and T. JORDAN's 
c     cal-coded solvers can be used. note: except for linpack band sol-
c     ver, factorizers and solvers have no pivoting.              
c.......................................................................
c     called by: step
c***********************************************************************
c
      include 'titan.imp'
      include 'titan.par'
      include 'titan.com'
c
      dimension emat(meqn, meqn, mgr, 5)
c
      equivalence (em2, emat)
c
c=======================================================================
c     reverse sign of rhs and choose method of solution
c=======================================================================
c
      do 2 i = 1, meqn
      do 1 k = 1, mgr
      rhs(i, k) = - rhs(i, k)
    1 continue
    2 continue
c
c     choose method of solution
c
c+++++++++++++++++++++++++++ lband .eq. 0 ++++++++++++++++++++++++ begin
      if (lband .eq. 0) then 
c
c=======================================================================
c     pentadiagonal system
c=======================================================================
c-----------------------------------------------------------------------
c     forward elimination
c-----------------------------------------------------------------------
c
      do 7 k = ngrs - 2, ngre + 2
c
c.......................................................................
c     construct d
c.......................................................................
c
      call scopy( meqn * meqn, ep2(1, 1, k-2), 1, tt, 1 )
c
      call sgemm('n','n', neqn, neqn, neqn, 1.0D0, ep1(1, 1, k-2), meqn,
     .                                             ep1(1, 1, k-1), meqn,
     .                                      1.0D0, tt            , meqn)
c
      call sgemm('n','n', neqn, neqn, neqn, 1.0D0, em2(1, 1, k  ), meqn,
     .                                             tt            , meqn,
     .                                      1.0D0, e00(1, 1, k  ), meqn)
c
      call sgemm('n','n', neqn, neqn, neqn,-1.0D0, em1(1, 1, k  ), meqn,
     .                                             ep1(1, 1, k-1), meqn,
     .                                     -1.0D0, e00(1, 1, k  ), meqn)
c
c.......................................................................
c     choose factorizer appropriate to computer architecture
c.......................................................................
c
c     non-cray computer
c
      if (lcray .eq. 0) call sgefad(e00(1, 1, k), meqn, neqn, info)
c
c     cray computer
c
      if (lcray .eq. 1) call sgefaj(e00(1, 1, k), meqn, neqn, info)
c
      if (info .ne. 0) then
            write( itty, 100 ) info, k
            write( iout, 100 ) info, k
            stop 'matslv1'
      end if
c
c.......................................................................
c     construct a(k)
c.......................................................................
c
      call scopy( meqn * meqn, em1(1, 1, k), 1, tt, 1 )
c
      call sgemm('n','n', neqn, neqn, neqn, 1.0D0, em2(1, 1, k  ), meqn,
     .                                             ep1(1, 1, k-2), meqn,
     .                                      1.0D0, tt            , meqn)
c
      call sgemm('n','n', neqn, neqn, neqn, 1.0D0, tt            , meqn,
     .                                             ep2(1, 1, k-1), meqn,
     .                                      1.0D0, ep1(1, 1, k  ), meqn)
c
c.......................................................................
c     construct c(k)
c.......................................................................
c
      call scopy( meqn, rhs(1, k-2), 1, v, 1 )
c
      call sgemv( 'n', neqn, neqn,  1.0D0, ep1(1, 1, k-2), meqn, 
     .                                     rhs(1,    k-1), 1   , 
     .                              1.0D0, v             , 1   )
c
      call sgemv( 'n', neqn, neqn, -1.0D0, em2(1, 1, k  ), meqn, 
     .                                     v             , 1   ,
     .                              1.0D0, rhs(1,    k  ), 1   )
c
      call sgemv( 'n', neqn, neqn,  1.0D0, em1(1, 1, k  ), meqn, 
     .                                     rhs(1,    k-1), 1   ,
     .                             -1.0D0, rhs(1,    k  ), 1   )
c
c.......................................................................
c     solve all systems
c.......................................................................
c
c     non-cray computer
c
      if (lcray .eq. 0) then
c
	    do 5 j = 1, neqn
	    call sgesld(e00(1, 1, k), meqn, neqn, ep1(1, j, k), 0)
	    call sgesld(e00(1, 1, k), meqn, neqn, ep2(1, j, k), 0)
    5       continue
	    call sgesld(e00(1, 1, k), meqn, neqn, rhs(1,    k), 0)
      end if
c
c     cray computer
c
      if (lcray .eq. 1 .and. neqn .ne. meqn) then
c
        call sgeslmj(e00(1, 1, k), meqn, neqn, ep1(1, 1, k), meqn, neqn)
        call sgeslmj(e00(1, 1, k), meqn, neqn, ep2(1, 1, k), meqn, neqn)
	call sgesl j(e00(1, 1, k), meqn, neqn, rhs(1,    k))
      end if
c
      if (lcray .eq. 1 .and. neqn .eq. meqn) then
c
         call scopy( meqn*meqn, ep1(1, 1, k), 1, tt(1, 1, 1), 1 )
         call scopy( meqn*meqn, ep2(1, 1, k), 1, tt(1, 1, 2), 1 )
         call scopy( meqn     , rhs(1,    k), 1, tt(1, 1, 3), 1 )
c
         call sgeslmj(e00(1, 1, k), meqn, meqn, tt, meqn, 2 * meqn + 1)
c
         call scopy( meqn*meqn, tt(1, 1, 1), 1, ep1(1, 1, k), 1 )
         call scopy( meqn*meqn, tt(1, 1, 2), 1, ep2(1, 1, k), 1 )
         call scopy( meqn     , tt(1, 1, 3), 1, rhs(1,    k), 1 )
      end if
c
    7 continue
c
c-----------------------------------------------------------------------
c     back substitution
c-----------------------------------------------------------------------
c
      do 8 k = ngre + 1, ngrs - 2, -1
      call sgemv( 'n', neqn, neqn, 1.0D0, ep2(1, 1, k  ), meqn, 
     .                                    rhs(1,    k+2), 1   ,
     .                             1.0D0, rhs(1,    k  ), 1   )
c
      call sgemv( 'n', neqn, neqn, 1.0D0, ep1(1, 1, k  ), meqn, 
     .                                    rhs(1,    k+1), 1   ,
     .                             1.0D0, rhs(1,    k  ), 1   )
    8 continue
c
      return
c
c+++++++++++++++++++++++++++ lband .eq. 0 ++++++++++++++++++++++++++ end
c
      end if
c
c=======================================================================
c     band system
c=======================================================================
c
c+++++++++++++++++++++++++++ lband .eq. 1 ++++++++++++++++++++++++ begin
      if (lband .eq. 1) then
c
c-----------------------------------------------------------------------
c     map pentadiagonal system onto band system
c-----------------------------------------------------------------------
c
      npd = 3 * neqn
      nld = npd - 1
      nud = npd - 1
      ncl = neqn * (ngre - ngrs + 5)
c
      do 12 l = 1, 5
      ll = l - 3
c
         do 11 i = 1, neqn
         do 10 j = 1, neqn
c
c           non-cray computer
	    if(lcray .eq. 0) m = i - j - ll * neqn + 2 * npd - 1
c
c           cray computer
	    if(lcray .eq .1) m = i - j - ll * neqn +     npd
c
	                     kmin = ngrs - 2
	    if( ll .eq. -1 ) kmin = ngrs - 1
	    if( ll .eq. -2 ) kmin = ngrs
	                     kmax = ngre + 2
	    if( ll .eq.  1 ) kmax = ngre + 1
	    if( ll .eq.  2 ) kmax = ngre
c
            do 9 k = kmin, kmax
               n = j + neqn * (k - ngrs + 2 + ll)
               bm(m, n) = emat(i, j, k, l)
    9       continue
c
   10    continue
   11    continue
c
   12 continue
c
      do 14 i = 1, neqn
      do 13 k = ngrs - 2, ngre + 2
      j = i + neqn * (k - ngrs + 2)
      br(j) = rhs(i, k)
   13 continue
   14 continue
c
c-----------------------------------------------------------------------
c     factor band matrix
c-----------------------------------------------------------------------
c
c     non-cray computer
c
      if(lcray .eq. 0) call sgbfa(bm, 3*mpd-2, ncl, nld, nud, ipv, ifl)
c
c     cray computer
c
      if(lcray .eq. 1) call bglsdc(ncl, nld, nud, bm, 3*mpd-2, ifl)
c
      if( ifl .gt. 0 ) then
            write( itty, 200 ) ifl
            write( iout, 200 ) ifl
            stop 'matslv2'
      end if
c
c-----------------------------------------------------------------------
c     solve band system
c-----------------------------------------------------------------------
c
c     non-cray computer
c
      if(lcray .eq. 0) call sgbsl(bm, 3*mpd-2, ncl, nld, nud, ipv, br,0)
c
c     cray computer
c
      if(lcray .eq. 1) call bglssl(ncl, nld, nud, bm, 3*mpd-2, br, ifl)
c
      if( ifl .gt. 0 ) then
            write( itty, 300 ) ifl
            write( iout, 300 ) ifl
            stop 'matslv3'
      end if
c
c-----------------------------------------------------------------------
c     map band solution back onto pentadiagonal system
c-----------------------------------------------------------------------
c
      do 16 i = 1, neqn
      do 15 k = ngrs - 2, ngre + 2
      j = i + neqn * (k - ngrs + 2)
      rhs(i, k) = br(j)
   15 continue
   16 continue
c
      return
c
c+++++++++++++++++++++++++++ lband .eq. 1 ++++++++++++++++++++++++++ end
c
      end if
c
c-----------------------------------------------------------------------
c
  100 format(/' matslv100,  info ='i5' k ='i5 )
  200 format(/' matslv200,  iflag ='i5 )
  300 format(/' matslv300,  iflag ='i5 )
c
      end
