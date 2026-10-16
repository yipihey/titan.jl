      subroutine readeos(fname, f, fx, fy, fxy)
c
c***********************************************************************
c     read eos tables
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      logical lxst
c
      character*8 fname, dname, pad
c
      dimension f(mxe, mye), fx(mxe, mye), fy(mxe, mye), fxy(mxe, mye)
      dimension pad(12)
c
c=======================================================================
c     open the file
c=======================================================================
c
      inquire( file = fname, exist = lxst )
      if (.not. lxst) go to 10
c
      open( unit = ieos, iostat = ier, file = fname )
      if (ier .gt. 0) go to 12
c
c=======================================================================
c     read the data                                                    
c=======================================================================
c
      read(ieos,'(13a8,3e20.12)') dname, pad, xabun, yabun, zabun
      if (dname .ne. fname) go to 14
c
      read(ieos,'(2(2e20.12,i5))') xxmin, dxx, nx, yymin, dyy, ny
      if (nx    .ne. mxe   .or. ny    .ne. mye  ) go to 16
      if (xxmin .ne. xmine .or. yymin .ne. ymine) go to 18
      if (dxx   .ne. dxe   .or. dyy   .ne. dye  ) go to 20
c
      read(ieos,'(4e20.12)') f, fx, fy, fxy
c
c=======================================================================
c     close the file
c=======================================================================
c
      close(unit = ieos, iostat = ier)
      if (ier .gt. 0) go to 22
c
      return
c
c=======================================================================
c     error exits                                                      
c=======================================================================
c
   10 write (itty, 11) fname
      write (iout, 11) fname
   11 format(' eos dataset 'a10' fails to exist')
      stop 'readeos1'
c
   12 write (itty, 13) fname, ier, ier, ier
      write (iout, 13) fname, ier, ier, ier
   13 format(' error opening eos file 'a8' ier =' i20, o22, a8)
      stop 'readeos2'
c
   14 write (itty, 15) fname, dname
      write (iout, 15) fname, dname
   15 format(' file-name discrepancy. fname = 'a8' dname = 'a8)
      stop 'readeos3'
c
   16 write (itty, 17) mxe, nx, mye, ny
      write (iout, 17) mxe, nx, mye, ny
   17 format(' discrepancy in dimensions. mx ='i3' nx ='i3
     .                                  ' my ='i4' ny ='i4)
      stop 'readeos4'
c
   18 write (itty, 19) xmine, xxmin, ymine, yymin
      write (iout, 19) xmine, xxmin, ymine, yymin
   19 format(' discrepant origin. xmine ='1pe20.12' xxmin ='e20.12
     .                          ' ymine ='1pe20.12' yymin ='e20.12)
      stop 'readeos5'
c
   20 write (itty, 21) dxe, dxx, dye, dyy
      write (iout, 21) dxe, dxx, dye, dyy
   21 format(' discrepant spacing. dxe ='1pe20.12' dxx ='e20.12
     .                           ' dye ='1pe20.12' dyy ='e20.12)
      stop 'readeos6'
c
   22 write (itty, 23) fname, ier, ier, ier
      write (iout, 23) fname, ier, ier, ier
   23 format(' error closing eos file 'a8' ier =' i20, o22, a8)
      stop 'readeos7'
c
      end
