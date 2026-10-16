      subroutine pardoc
c
c***********************************************************************
c     print the values of all parameters set for the current run
c.......................................................................
c     calling sequence: titan > start > start{01, 02, ....} > pardoc
c***********************************************************************
c
      include'titan.par'
      include'titan.imp'
c
      character*80 header
c
      common /agrid / x2(((3*mad + 6)*meqn + mad + 21)*mgr + 2*mad + 3),
     .                i2
c
      common /bc    / x4(26)
c
      common /const / x5(12)
c
      common /edding/ x8(4*mang + 7*mgr + 8*(2 * mgr + mcor) + 1), 
     .                i8(2)
c
      common /energy/ x9(8*mgr + 10)
c
      common /eostab/ x10(12*mgr + 12*mxe*mye + 5)
c
      common /geom  / x11(4), i11(3)
c
      common /hydro / x12(11)
c
      common /index / i13(17)
c
      common /io    / idoc, i15(14)
c
      common /logic / i17(mad + 24)
c
      common /matrix/ x18(meqn*(mgr*(5*meqn + 3*mpd) + (3*meqn + 1))),
     .                i18(meqn* mgr + 3)
c
      common /rad   / x24(mgr + 3)
c
      common /star  / x26(18), header 
c
      write (idoc,'(//" agrid"/)')
      write (idoc, 100) (x2(j), j = 1, 2 * mad + 3)
      write (idoc, 101) i2
c
      write (idoc,'(//" bc"/)')
      write (idoc, 100) x4
c
      write (idoc,'(//" const"/)')
      write (idoc, 100) x5
c
      jj = 3 * mang + 7 * mv + 8 * (2 * mv + mcor) + 1
      write (idoc,'(//" edding"/)')
      write (idoc, 100) (x8(j), j = jj, jj + mang)
      write (idoc, 101) i8
c
      write (idoc,'(//" energy"/)')
      write (idoc, 100) (x9(j), j = 1, 5)
c
      write (idoc,'(//" eostab"/)')
      write (idoc, 100) (x10(j), j = 1, 5)
c
      write (idoc,'(//" geom"/)')
      write (idoc, 100) x11
      write (idoc, 101) i11
c
      write (idoc,'(//" hydro"/)')
      write (idoc, 100) x12
c
      write (idoc,'(//" index"/)')
      write (idoc, 101) i13
c
      write (idoc,'(//" io"/)')
      write (idoc, 101) idoc, i15
c
      write (idoc,'(//" logic"/)')
      write (idoc, 101) i17
c
      jj = meqn * mgr
      write (idoc,'(//" matrix"/)')
      write (idoc, 101) (i18(j), j = jj + 1, jj + 3)
c
      write (idoc,'(//" rad"/)')
      write (idoc, 100) (x24(j), j = mgr + 1, mgr + 3)
c
      write (idoc,'(//" star"/)')
      write (idoc, 100) x26
      write (idoc, 102) header
c
  100 format(1p5e16.8)
  101 format(5i16)
  102 format(a80)
c
      return
      end
