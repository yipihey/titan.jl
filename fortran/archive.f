      subroutine archive
c
c***********************************************************************
c     save the solution on the archive tape
c.......................................................................
c     called by: titan
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      common /hdf/ hdfname
      character*8  hdfname
c
          jmod = mod( jstep - jsteps + 1, ndump )
      if (jmod .eq. 0 .or. jstep .eq. jstepe 
     .                .or. timen .ge. tmax  ) then
c
	    jdump = jdump + 1
c
            if (idump .gt. 0) then
	                             irec = 2
	        write ( idump, rec = irec, iostat = ier) jdump, jstep
c
                if (ier .gt. 0) then
                    write (itty, 21) irec, ier, ier, ier
                    write (iout, 21) irec, ier, ier, ier
                    stop 'archive01' 
                end if
	                             irec  = jdump + 2
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
                    stop 'archive02' 
                end if
            end if
c
	    write (itty, 20) jdump, jstep
            write (iout, 20) jdump, jstep
c
            call hstore(hdfname)
c
       end if
c
c=======================================================================
c
   20 format(/' model jdump ='i6' jstep ='i6' archived')
   21 format(/'  unit idump. irec = 'i3' ier ='i20, o22, a8)
c
      return
      end
