c-----------------------------------------------------------------------
c     PORT phase 2: the console/.out messages and STOPs of step.f, so a
c     Julia re-implementation of step's control flow writes identical
c     output.  icode: 1 divergence, 2 failure (writout, clzfil, STOP),
c     3 changes too large, 4 too many repeats (clzfil, STOP), 5 repeat.
c-----------------------------------------------------------------------
      subroutine port_msg(icode, xs)
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      if (icode .eq. 1) then
      write( itty,'(/34x " divergence. jtry = " i3)') jtry
      write( iout,'(/34x " divergence. jtry = " i3)') jtry
      end if
      if (icode .eq. 2) then
          write(iout,'(" failure to converge after"i3" attempts")') ntry
          call writout
          call clzfil
          stop 'step01'
      end if
      if (icode .eq. 3) then
      write(itty, 101) xs
      write(iout, 101) xs
      end if
      if (icode .eq. 4) then
          write(itty, 102) nback
          write(iout, 102) nback
          call clzfil
          stop 'step02'
      end if
      if (icode .eq. 5) then
      write(itty,103) dtime
      write(iout,103) dtime
      end if
      return
  101 format(/' changes over timestep'1pe9.2' times larger than limit')
  102 format( ' total changes still too large after'i2' attempts')
  103 format( ' repeat integration step with dtime ='1pe9.2)
      end
c
c-----------------------------------------------------------------------
c     PORT phase 3: optional override of ncell from the environment
c     variable TITAN_NCELL (no effect when unset).
c-----------------------------------------------------------------------
      subroutine port_env_ncell(n)
      integer n
      character*16 s
      s = ' '
      call getenv('TITAN_NCELL', s)
      if (s .ne. ' ') read(s, *) n
      return
      end
