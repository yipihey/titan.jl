      program titantr
c-----------------------------------------------------------------------
c     PORT phase 1: identical to titan.f, plus a raw per-step dump of
c     common /new/ (all 41 arrays, declaration order) for the first
c     ntrace steps to trace.bin (stream access).  Used only to verify
c     the Julia-driven library run step by step.
c-----------------------------------------------------------------------
      include'titan.imp'
      include'titan.par'
      include'titan.com'
      include'titan.dat'
      parameter (ntrace = 50, itr = 77)
      xindef = 0.0D0
      call opnfil
      call start
      open(unit=itr, file='trace.bin', access='stream',
     .     form='unformatted', status='replace')
      do 1 jstep = jsteps, jstepe
      call step
      if (jstep - jsteps + 1 .le. ntrace) write(itr) jstep, iter,
     .  timen, dtime,
     .  rn,rmun,rmup1n,rmum1n,rcn,xmn,xmen,xmcn,dvoln,dn,dsn,
     .  un,usn,tn,asn,egn,egsn,pgn,xnen,ern,ersn,egrn,egrsn,
     .  etn,etsn,egtn,egtsn,frn,frsn,egrtn,egrtsn,
     .  avchin,betan,delsn,cpn,chifn,xken,xkpn,plfn,bri,gamman
      call archive
      call writout
      call writsty
      if (timen .ge. tmax) then
                           close(itr)
                           call clzfil
                           stop
      end if
    1 continue
      close(itr)
      call clzfil
      stop
      end
