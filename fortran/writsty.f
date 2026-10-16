      subroutine writsty
c
c***********************************************************************
c     evolution of selected physical variables
c.......................................................................
c     called by: titan, start{01,02,03,...,10}
c***********************************************************************
c
      include'titan.imp'
      include'titan.par'
      include'titan.com'
c
      real*8  odepth(mgr), flux(mgr), trad(mgr)
c
c=======================================================================
c     auxiliary variables
c=======================================================================
c
            odepth(ngre+2) = 0.0D0
            odepth(ngre+1) = 0.0D0
            do 20 k = ngre, ngrs - 1, -1
            odepth(k) = odepth(k+1) + dn(k) * chifn(k) * (rn(k+1)-rn(k))
            flux  (k) = frn(k) + un(k) * ern(k) * (1.D0 + fedd(k))
            trad  (k) = sqrt(sqrt(abs(ern(k)/car)))*ern(k)/abs(ern(k))
   20       continue
c
      if (lprob .ge. 8 .and. lprob .le. 10) then
 	    kq   = isamax(ngre - ngrs + 1, qf(ngrs), 1)
 	    rshk =   rn  (ngrs + kq - 1)
 	    dshk =   dn  (ngrs + kq - 1)
 	    tshk =   tn  (ngrs + kq - 1)
 	    pshk =   pgn (ngrs + kq - 1)
 	    vshk =   un  (ngrs + kq - 1)
 	    fshk =   frn (ngrs + kq - 1) * rshk**2 * 4.0D0 * cpi
 	    flsh =   flux(ngrs + kq - 1) * rshk**2 * 4.0D0 * cpi
      end if
c
c=======================================================================
c     problem 1: first FAVOURITE PROBLEM
c=======================================================================
c
      if (lprob .eq. 1) return
c
c=======================================================================
c     problem 2 - 5:
c=======================================================================
c
      if (lprob .ge. 2 .and. lprob .le. 5) return
c
c=======================================================================
c     problem 6 & 7: RADIATIVE EQUILIBRIUM
c=======================================================================
c
      if (lprob .eq. 6 .or. lprob .eq. 7) then
c
          if (jstep .eq. jsteps) write(ihist, 1) header 
          if (jstep .eq. jsteps) write(ihist, 3) 
	  write(ihist, 4)   jstep, timen, dtime,
     .                      4.0D0*cpi*rmun(ngre)*frn(ngre), 
     .                      ern(ngre-1), trad(ngre-1), tn(ngre-1), 
     .                      odepth(ngre-1)
      return
      end if
c
c=======================================================================
c     problem 8 & 9: SUB- AND SUPERCRITIAL SHOCK WAVES
c=======================================================================
c
      if (lprob .eq. 8 .or. lprob .eq. 9) then
c
          if (jstep .eq. jsteps) write(ihist, 1) header 
          if (jstep .eq. jsteps) write(ihist, 5) 
	  write(ihist, 4)   jstep, timen, 
     .                      rn(ngre), (xmass-xmen(ngre))*4.D0*cpi,
     .                      un(ngre), tn(ngre), pgn(ngre), 
     .                      4.0D0*cpi*rmun(ngre)* frn(ngre)/clsol
      return
      end if
c
c=======================================================================
c     problem 10: RADIATIVE BLAST WAVE
c=======================================================================
c
      if (lprob .eq. 10) then
c
          if (jstep .eq. jsteps) write(ihist, 1) header 
          if (jstep .eq. jsteps) write(ihist, 6) 
	  write(ihist, 8)   jstep, timen, 
     .                      4.0D0*cpi*rmun(ngre)* frn(ngre)/clsol,
     .                      un(ngre),   tn(ngre),  dn(ngre), 
     .                      odepth(ngre), dmdt(ngre)/cmsol,
     .                      rshk, dshk, vshk, tshk, fshk/clsol, 
     .                      flsh/clsol
      return
      end if
c
c=======================================================================
c     problem 99: STELLAR PULSATION
c=======================================================================
c
      if (lprob .ge. 11 .and. lprob .lt. 99) return
c
      if (lprob .eq. 99) then
c
           if (jstep.eq.2) then
                             write(ihist, 1) header
                             write(ihist, 9)
           end if
           write(ihist, 10)   jstep, timen,
     .                       4.0D0*cpi*rmun(ngre)* frn(ngre)/clsol,
     .                       un(ngre),   tn(ngre),  dn(ngre),
     .                       odepth(ngre), dmdt(ngre)/cmsol
      return
      end if
c
c=======================================================================
c
    1 format(/a80)
    2 format(/'initial free fall time is'1pe11.4' seconds')
    3 format(/'jstep       timen       dtime  luminosity         E_r',
     .           '         T_r         T_g         tau'              )
    4 format(i5,1p7e12.4)
    5 format(/'jstep       timen       R           M           u'    ,
     .                     '           T           P           L'    )
    6 format(/'jstep       timen       L           u           T'    ,
     .                     '           d         tau       dM/dt'    ,
     . '       R_s       d_s       u_s       T_s       L_s       L_c')
    7 format(/'jstep       timen/tff0  dtime       L_com       L_lab',
     .           '         u           T_r         mass'             )
    8 format(i5,1p7e12.4,6e10.3)
    9 format(/'jstep       timen       L           u           T'    ,
     .                     '           d         tau       dM/dt'    )
   10 format(i5,1p7e12.4)
c
      end
