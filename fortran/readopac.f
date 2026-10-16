      subroutine readopac(fname, f, fx, fy, fxy)
c
c***********************************************************************
c     read opacity tables
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      logical lxst
c
      character*8 fname, dname
c
      dimension f(mxo, myo), fx(mxo, myo), fy(mxo, myo), fxy(mxo, myo)
c
c=======================================================================
c     open the file
c=======================================================================
c
      inquire( file = fname, exist = lxst )
      if (.not. lxst) go to 10
c
      open( unit = iopac, iostat = ier, file = fname )
      if (ier .gt. 0) go to 12
c
c=======================================================================
c     read the data                                                    
c=======================================================================
c
      read(iopac,'(13a8)') dname
      if (dname .ne. fname) go to 14
c
      read(iopac,'(2(2e20.12,i5))') xxmin, dxx, nx, yymin, dyy, ny
      if (nx    .ne. mxo   .or. ny    .ne. myo  ) go to 16
      if (xxmin .ne. xmino .or. yymin .ne. ymino) go to 18
      if (dxx   .ne. dxo   .or. dyy   .ne. dyo  ) go to 20
c
      read(iopac,'(4e20.12)') f, fx, fy, fxy
c
c=======================================================================
c     close the file
c=======================================================================
c
      close(unit = iopac, iostat = ier)
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
      stop 'readopa1'
c
   12 write (itty, 13) fname, ier, ier, ier
      write (iout, 13) fname, ier, ier, ier
   13 format(' error opening eos file 'a8' ier =' i20, o22, a8)
      stop 'readopa2'
c
   14 write (itty, 15) fname, dname
      write (iout, 15) fname, dname
   15 format(' file-name discrepancy. fname = 'a8' dname = 'a8)
      stop 'readopa3'
c
   16 write (itty, 17) mxo, nx, myo, ny
      write (iout, 17) mxo, nx, myo, ny
   17 format(' discrepancy in dimensions. mx ='i3' nx ='i3
     .                                  ' my ='i4' ny ='i4)
      stop 'readopa4'
c
   18 write (itty, 19) xmino, xxmin, ymino, yymin
      write (iout, 19) xmino, xxmin, ymino, yymin
   19 format(' discrepant origin. xmino ='1pe20.12' xxmin ='e20.12
     .                          ' ymino ='1pe20.12' yymin ='e20.12)
      stop 'readopa5'
c
   20 write (itty, 21) dxo, dxx, dyo, dyy
      write (iout, 21) dxo, dxx, dyo, dyy
   21 format(' discrepant spacing. dxo ='1pe20.12' dxx ='e20.12
     .                           ' dyo ='1pe20.12' dyy ='e20.12)
      stop 'readopa6'
c
   22 write (itty, 23) fname, ier, ier, ier
      write (iout, 23) fname, ier, ier, ier
   23 format(' error closing opac file 'a8' ier =' i20, o22, a8)
      stop 'readopa7'
c
      end
