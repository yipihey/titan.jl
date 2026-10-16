      subroutine opnfil
c
c***********************************************************************
c     open i/o files
c.......................................................................
c     called by: titan
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      parameter( nlabels = 12 )
c
      common /hdf/ hdfname
c
      character*8  hdfname, dname, fname, pname, status
      character*70 annote(nlabels)
c
      logical lxst
c
c=======================================================================
c     initialize i/o unit numbers, set direct-access record lengths
c=======================================================================
c
      itty  = 4  ! 6 on SGIs
      iin   = 5  ! 1 on SGIs
      iout  = 0  ! 4 on SGIs
      idoc  = 0
      ieos  = 0
      iopac = 0
      idump = 0
      ihist = 0
      inrg  = 0
c
      ldump = 8 * (45 * mgr + 19)
c-----lhist = ##########
c
      dname = '********'
      fname = '********'
      pname = '********'
c
c=======================================================================
c     open terminal and input files
c=======================================================================
c-----------------------------------------------------------------------
c     open terminal
c-----------------------------------------------------------------------
c
      istp = 1
      dname = '/dev/tty'
      open(unit = itty, iostat = ier, file = dname)
      if (ier .gt. 0) go to 4
c
c-----------------------------------------------------------------------
c     open input file
c-----------------------------------------------------------------------
c
      istp = 2
      dname = 'titan.in'
      inquire(file = dname, exist = lxst)
      if (.not. lxst) go to 6
c
      istp = 3
      open(unit = iin , iostat = ier, file = dname)
      if (ier .gt. 0) go to 4
c
c=======================================================================
c     read file specifications, activate i/o unit numbers, open files
c=======================================================================
c-----------------------------------------------------------------------
c     read number of datasets to be handled
c-----------------------------------------------------------------------
c
      read (iin,'(i5)') nds
c
c-----------------------------------------------------------------------
c     loop over all datasets
c-----------------------------------------------------------------------
c
      do 3 i = 1, nds
c
      read (iin,'(i5,3(2x,a8))') iunit, dname, pname, status
c
      if (iunit .le. 0) go to 3
c
c-----------------------------------------------------------------------
c     open output file (sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. '  output') then
            iout = iunit
c
            istp = 4
            open(unit = iout , iostat = ier, file = pname)
            if (ier .gt. 0) go to 4
      end if
c
c-----------------------------------------------------------------------
c     open and read eos files (old, sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. '     eos') then
            ieos = iunit
c
            istp = 5
            if (status .ne. '     old') go to 10
	    fname = pname//'.pg'
	    call readeos(fname, fpg, fpgx, fpgy, fpgxy)
	    fname = pname//'.eg'
	    call readeos(fname, feg, fegx, fegy, fegxy)
	    fname = pname//'.pe'
	    call readeos(fname, fpe, fpex, fpey, fpexy)
      end if
c
c-----------------------------------------------------------------------
c     open and read opacity files (old, sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. ' opacity') then
            iopac = iunit
c
            istp = 6
            if (status .ne. '     old') go to 10
	    fname = pname//'.cf'
	    call readopac(fname, fcf, fcfx, fcfy, fcfxy)
c#######################################################################
c     skip PLANCK mean tables until available
c           fname = pname//'.ke'
c           call readopac(fname, fke, fkex, fkey, fkexy)
c#######################################################################
      end if
c
c-----------------------------------------------------------------------
c     open dump file (direct access)
c-----------------------------------------------------------------------
c
      if (dname .eq. '    dump') then
            idump = iunit
c
            istp = 7
            if (status .eq. '     old') then
		  inquire(file = pname, exist = lxst)
		  if (.not. lxst) go to 6
	    end if
c
            istp = 8
            open(unit = idump, iostat = ier, file = pname, 
     .            access = 'direct', recl = ldump, form = 'unformatted')
            if (ier .gt. 0) go to 4
      end if
c
c-----------------------------------------------------------------------
c     open history file (sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. ' history') then
            ihist = iunit
c
            istp = 9
            open(unit = ihist, iostat = ier, file = pname)
            if (ier .gt. 0) go to 4
      end if
c
c-----------------------------------------------------------------------
c     open documentation file (sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. '     doc') then
            idoc = iunit
c
            istp = 10
            open(unit = idoc , iostat = ier, file = pname)
            if (ier .gt. 0) go to 4
      end if
c
c-----------------------------------------------------------------------
c     open energy file (sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. '  energy') then
            inrg = iunit
c
            istp = 11
            open(unit = inrg , iostat = ier, file = pname)
            if (ier .gt. 0) go to 4
      end if
c
c-----------------------------------------------------------------------
c     open tty rerouting file (sequential)
c-----------------------------------------------------------------------
c
      if (dname .eq. '    rtty') then
            itty = iunit
c
            istp = 12
            open(unit = itty , iostat = ier, file = pname)
            if (ier .gt. 0) go to 4
      end if
c
c-----------------------------------------------------------------------
c     open HDF file
c-----------------------------------------------------------------------
c
      if (dname .eq. '     hdf') then
c
      annote( 1) = '%code TITAN radiation hydrodynamic simulations'
      annote( 2) = '%author '
      annote( 3) = '%date '
      annote( 4) = '%0desc '
      annote( 5) = '%0desc '
      annote( 6) = '%0desc '
      annote( 7) = '%1desc '
      annote( 8) = '%1desc '
      annote( 9) = '%1desc '
      annote(10) = '%1desc '
      annote(11) = '%2desc '
      annote(12) = '%2desc '
c
                   hdfname = pname
      call hdfinit(hdfname, annote, nlabels, ier)
c
      if (ier .ne. 1) go to 12
      end if
c
    3 continue
c
c-----------------------------------------------------------------------
c     check that output file was assigned
c-----------------------------------------------------------------------
c
      if (iout .le. 0) go to 8
c
      return
c
c=======================================================================
c     error exits                                                      
c=======================================================================
c
    4 write (itty, 5) istp, ier, ier, ier, dname, pname
    5 format(' error return in step' i3 ' of file activation'/
     .             ' ier =' i20, o22, a8/ ' dname = ' a8 ' pname = ' a8)
      if (iout .gt. 0) write (iout, 5) istp, ier, ier, ier, dname, pname
      stop 'opnfil1'
c
    6 write (itty, 7) dname, pname, istp
    7 format(' dataset ' 2a10 ' fails to exist. istp =' i3)
      if (iout .gt. 0) write (iout, 7) dname, pname, istp
      stop 'opnfil2'
c
    8 write (itty, 9)
    9 format(' no output file specified')
      stop 'opnfil3'
c
   10 write (itty, 11) status, dname, pname, istp
   11 format(' status of dataset ' 2a10 ' is not "old". istp =' i3)
      stop 'opnfil4'
   12 write (itty, 13) pname
   13 format(' initialization of' a10 'failed')
c
      end
