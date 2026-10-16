      subroutine reset
c
c***********************************************************************
c     return to solution before last timestep in case of divergence or
c     decision to integrate over with smaller timestep
c.......................................................................
c     called by: step
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
c=======================================================================
c     replace solution at new time with solution at old time
c=======================================================================
c
      timen = timeo
c
      tes = teso
      tew = tewo
      tel = telo
      teq = teqo
c
      do 1 k = ngrs - 1, ngre + 2
      tn    (k) = to    (k)
      dn    (k) = do    (k)
      rn    (k) = ro    (k)
      rmun  (k) = rmuo  (k)
      rmup1n(k) = rmup1o(k)
      rmum1n(k) = rmum1o(k)
      dvoln (k) = dvolo (k)
      xmen  (k) = xmeo  (k)
      xmn   (k) = xmo   (k)
      un    (k) = uo    (k)
      egn   (k) = ego   (k)
      etn   (k) = eto   (k)
      ern   (k) = ero   (k)
      frn   (k) = fro   (k)
      drdt  (k) = drdto (k)
      dmdt  (k) = dmdto (k)
      urel  (k) = urelo (k)
      xnt   (k) = xnto  (k)
    1 continue
c
      do 2 k = ngrs - 1, ngre + 1
      dsn   (k) = dso   (k)
      usn   (k) = uso   (k)
      pgn   (k) = pgo   (k)
      asn   (k) = aso   (k)
      egsn  (k) = egso  (k)
      etsn  (k) = etso  (k)
      ersn  (k) = erso  (k)
      egrn  (k) = egro  (k)
      egrsn (k) = egrso (k)
      egtn  (k) = egto  (k)
      egtsn (k) = egtso (k)
      egrtn (k) = egrto (k)
      egrtsn(k) = egrtso(k)
      frsn  (k) = frso  (k)
      chifn (k) = chifo (k)
      xken  (k) = xkeo  (k)
      xkpn  (k) = xkpo  (k)
      xnen  (k) = xneo  (k)
      plfn  (k) = plfo  (k)
      cpn   (k) = cpo   (k)
      betan (k) = betao (k)
      delsn (k) = delso (k)
    2 continue
c
      do 3 k = ngrs, ngre + 1
      avchin(k) = avchio(k)
    3 continue
c
      return
      end
