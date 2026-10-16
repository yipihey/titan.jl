c
c--------------------------------------------------------------titan.com
c
      character*80  header
c
      common /advec/ 
     .       adv(mgr, madv)
c
      common /agrid/
     .       xnc     (mgr), xnt     (mgr), xnu     (mgr), ss      (mgr), 
     .                      xnto    (mgr), xnuo    (mgr), rr      (mgr),
     .
     .                      dnudlr00(mgr), dnudlrp1(mgr),
     .       dncdlrm1(mgr), dncdlr00(mgr), dncdlrp1(mgr), dncdlrp2(mgr),
     .       dntdlrm1(mgr), dntdlr00(mgr), dntdlrp1(mgr), dntdlrp2(mgr),
     .       dabdlrm1(mgr), dabdlr00(mgr), dabdlrp1(mgr), dabdlrp2(mgr),
     .
     .       dcddlx00(mgr, meqn), dssdlx00(mgr, meqn),
     .       dcddlxp1(mgr, meqn), dssdlxp1(mgr, meqn),
     .       dcddlxp2(mgr, meqn), dssdlxp2(mgr, meqn),
     .
     .        cs     (mgr, mad      ),
     .       dcsdlx00(mgr, mad, meqn), dcsdlxp1(mgr, mad, meqn),
     .       dcsdlxp2(mgr, mad, meqn),
     .
     .       alph, tau, xscale, yscl(mad), wt(mad), ibet
c
      common /avuor/ 
     .       aur     (mgr), 
     .       dardlrp1(mgr), dardlr00(mgr), dardlup1(mgr), dardlu00(mgr),
     .       dudr    (mgr), 
     .       durdlrp1(mgr), durdlr00(mgr), durdlup1(mgr), durdlu00(mgr)
c
      common /bc/
     .       phil0, phil1, phil2, phil3, delrl, delml, 
     .       phir0, phir1, phir2, phir3, delrr, delmr, 
     .       dl, ul, egasl, tl, etl, uextl,
     .       dr, ur, egasr, tr, etr, uextr, pextr,
     .       omega
c
      common /const/
     .       car  , cc   , cgas , cgrav, ck   , clsol, cm0  , cmsol,
     .       cpi  , crsol, csige, csigr
c
      common /diffus/
     .       qd      (mgr), qdn     (mgr), df      (mgr),
     .       ddfdlrm1(mgr), ddfdlr00(mgr), ddfdlrp1(mgr),
     .       ddfdlqm1(mgr), ddfdlq00(mgr),
     .       ddfdldm1(mgr), ddfdld00(mgr),
     .       ddfdltm1(mgr), ddfdlt00(mgr),
     . 
     .       sig, zet
c
      common /dot/
     .       drdt (mgr), dmdt (mgr), urel (mgr),
     .       drdto(mgr), dmdto(mgr), urelo(mgr)
c
      common /edding/
     .       ang(mang), cjk(mang), ejk(mang), sjk(mang),
     . 
     .       xmom0(mgr), xmom2(mgr), dtau (mgr), sf   (mgr),
     .       asf  (mgr), bsf  (mgr), dsf  (mgr),
     . 
     .       xim(2*mgr + mcor), xi0 (2*mgr + mcor), xip(2*mgr + mcor),
     .       sl (2*mgr + mcor), sr  (2*mgr + mcor), rp (2*mgr + mcor),
     .       p  (2*mgr + mcor), xang(2*mgr + mcor),
     .
     .       epsedf, nang, ncor
c
      common /energy/
     .       r0   (mgr), xm0  (mgr), d0    (mgr), u0(mgr), egrt0(mgr),
     .       dvol0(mgr), tescr(mgr), tescr0(mgr),
     . 
     .       tetot, tee, tes, teso, tew, tewo, tel, telo, teq, teqo
c
      common /eostab/
     .       feg  (mxe, mye), fpg  (mxe, mye), fpe  (mxe, mye),
     .       fegx (mxe, mye), fpgx (mxe, mye), fpex (mxe, mye),
     .       fegy (mxe, mye), fpgy (mxe, mye), fpey (mxe, mye),
     .       fegxy(mxe, mye), fpgxy(mxe, mye), fpexy(mxe, mye),
     .
     .       dlpgdldn(mgr), dlpgdltn(mgr), dlegdldn(mgr), dlegdltn(mgr),
     .       dlpedldn(mgr), dlpedltn(mgr), dlcpdldn(mgr), dlcpdltn(mgr),
     .       dlbedldn(mgr), dlbedltn(mgr), dldsdldn(mgr), dldsdltn(mgr),
     .
     .       gam, gmu, xabun, yabun, zabun
c
      common /geom/
     .       xmu, xmup1, xmum1, rxm1, mu, mup1, mum1
c
      common /hydro/
     .       thet , ql0 , ql1   , cq1 , cq2 , q0, 
     .       cqvis, cadv, epsadv, sigd, sige
c
      common /index/
     .       ir, im, id, iu, it, ie, if, ic,
     .       jr, jm, jd, ju, jt, je, jf, jc, jb
c
      common /integrat/ 
     .       sx(mgr, meqn), smax  , stol, timeo, timen, dtime, tfac,
     .       jstep, jsteps, jstepe, jext, next , jback, nback
c
      common /io/
     .       idoc, idump, ihist, iin  , inrg , ieos , iopac, iout , 
     .       itty, jdump, jhist, ldump, lhist, ndump, nout
c
      common /iterate/ 
     .       dx(meqn), cx(mgr), dmax, dtol, cmax, ctol, conv,
     .       kx(meqn), iter, niter, jtry, ntry
c
      common /logic/
     .       lady(mad)   , ladx ,
     .       leos , lopac, lcray, lband, lboos, ltur ,
     .       lprob, linit, lgeom, lhydr, lrad , ltran, lam  , 
     .       leibc, leobc, llibc, llobc, lribc, lrobc, lgrid, 
     .       ldriv, lhsre, ltalk
c
      common /matrix/
     .       em2(meqn, meqn, mgr), em1(meqn, meqn, mgr),
     .       e00(meqn, meqn, mgr), ep1(meqn, meqn, mgr),
     .       ep2(meqn, meqn, mgr), rhs(meqn,       mgr),
     .
     .       bm(3*mpd - 2, meqn*mgr), br(meqn*mgr),
     .
     .       tt(meqn, meqn,3), v(meqn), ipv(meqn*mgr),
     .
     .       ngrs, ngre, neqn
c
      common /mid/ 
     .       r     (mgr), rmu   (mgr), rmup1 (mgr), rmum1 (mgr),
     .       rc    (mgr), xm    (mgr), xme   (mgr), xmc   (mgr),
     .       dvol  (mgr), d     (mgr), ds    (mgr), db    (mgr),
     .       u     (mgr), ub    (mgr), t     (mgr), as    (mgr),
     .       eg    (mgr), egb   (mgr), pg    (mgr), xne   (mgr),
     .       er    (mgr), erb   (mgr), egr   (mgr), egrb  (mgr),
     .       et    (mgr), etb   (mgr), egt   (mgr), egtb  (mgr),
     .       fr    (mgr), frb   (mgr), egrt  (mgr), egrtb (mgr),
     .       frnom (mgr), avchi (mgr), pt    (mgr), ambda (mgr),
     .       ft    (mgr), fc    (mgr), sink  (mgr), source(mgr),
     .       ftnom (mgr), fcnom (mgr), floor (mgr), acous (mgr),
     .       chif  (mgr), xke   (mgr), xkp   (mgr), plf   (mgr),
     .       beta  (mgr), dels  (mgr), cp    (mgr), fedd  (mgr),
     .       unom  (mgr), 
     .
     .       geddl, geddr
c
      common /new/ 
     .       rn    (mgr), rmun  (mgr), rmup1n(mgr), rmum1n(mgr),
     .       rcn   (mgr), xmn   (mgr), xmen  (mgr), xmcn  (mgr),
     .       dvoln (mgr), dn    (mgr), dsn   (mgr),
     .       un    (mgr), usn   (mgr), tn    (mgr), asn   (mgr),
     .       egn   (mgr), egsn  (mgr), pgn   (mgr), xnen  (mgr),
     .       ern   (mgr), ersn  (mgr), egrn  (mgr), egrsn (mgr),
     .       etn   (mgr), etsn  (mgr), egtn  (mgr), egtsn (mgr),
     .       frn   (mgr), frsn  (mgr), egrtn (mgr), egrtsn(mgr),
     .       avchin(mgr), betan (mgr), delsn (mgr), cpn   (mgr),
     .       chifn (mgr), xken  (mgr), xkpn  (mgr), plfn  (mgr),
     .       bri   (mgr), gamman(mgr)
c
      common /old/ 
     .       ro    (mgr), rmuo  (mgr), rmup1o(mgr), rmum1o(mgr),
     .       rco   (mgr), xmo   (mgr), xmeo  (mgr), xmco  (mgr),
     .       dvolo (mgr), do    (mgr), dso   (mgr),
     .       uo    (mgr), uso   (mgr), to    (mgr), aso   (mgr),
     .       ego   (mgr), egso  (mgr), pgo   (mgr), xneo  (mgr),
     .       ero   (mgr), erso  (mgr), egro  (mgr), egrso (mgr),
     .       eto   (mgr), etso  (mgr), egto  (mgr), egtso (mgr),
     .       fro   (mgr), frso  (mgr), egrto (mgr), egrtso(mgr),
     .       avchio(mgr), betao (mgr), delso (mgr), cpo   (mgr),
     .       chifo (mgr), xkeo  (mgr), xkpo  (mgr), plfo  (mgr)
c
      common /opactab/
     .       fcf  (mxo, myo), fke  (mxo, myo),
     .       fcfx (mxo, myo), fkex (mxo, myo),
     .       fcfy (mxo, myo), fkey (mxo, myo),
     .       fcfxy(mxo, myo), fkexy(mxo, myo),
     .
     .       dlkedldn(mgr), dlkpdldn(mgr), dlcfdldn(mgr),
     .       dlkedltn(mgr), dlkpdltn(mgr), dlcfdltn(mgr), dplfdltn(mgr)
c
      common /r3o2/ 
     .        r3(mgr), dr3dlrp1(mgr), dr3dlr00(mgr),
     .        r1(mgr), dr1dlr00(mgr)
c
      common /rad/
     .        drd(mgr), xlam, xipl, ximr
c
      common /rdif/ 
     .        dif   (mgr), 
     .       ddifdrp(mgr), ddifdr0(mgr), ddifdrm(mgr),
     .                     ddifde0(mgr), ddifdem(mgr),
     .                     ddifdd0(mgr), ddifddm(mgr),
     .                     ddifdt0(mgr), ddifdtm(mgr)
c
      common /star/
     .       xmass, xmext, xmtot, xmrat, g, xme0, tff0,
     .       xlum, xrad, teff, chif0, tau0, ratio,
     .       xindef, rho0, tmax, rin, rout,
     .       header
c
      common /turb1/
     .       dladld00(mgr), dladldm1(mgr), dladlt00(mgr), dladltm1(mgr),
     .       dladlrp1(mgr), dladlr00(mgr), dladlrm1(mgr), dladlx00(mgr),
     .       dladlup1(mgr), dladlu00(mgr),
     .
     .       dftdld00(mgr), dftdldm1(mgr), dftdlt00(mgr), dftdltm1(mgr),
     .       dftdlrp1(mgr), dftdlr00(mgr), dftdlrm1(mgr), dftdlx00(mgr),
     .       dftdlc00(mgr), dftdlcm1(mgr),
     .
     .       dskdld00(mgr), dskdldm1(mgr), dskdlt00(mgr), dskdltm1(mgr),
     .       dskdlrp1(mgr), dskdlr00(mgr), dskdlrm1(mgr), dskdlx00(mgr),
     .       dskdlc00(mgr),
     .
     .       dfldld00(mgr), dfldldm1(mgr), dfldlt00(mgr), dfldltm1(mgr),
     .       dfldlrp1(mgr), dfldlr00(mgr), dfldlrm1(mgr), dfldlx00(mgr),
     .       dfldlc00(mgr)
c
      common /turb2/
     .       dptdldn (mgr), dptdlcn (mgr),
     .
     .       dacdldn (mgr), dacdlcn (mgr), dacdltn (mgr),
     .
     .       cmlt, ctvis, cdrag, cdiss, cflux, etflr, ntur
c
      common /viscos1/ 
     .       qe      (mgr), ql      (mgr), div     (mgr),
     .       qp      (mgr), qx      (mgr), qm      (mgr),
     .
     .       dqpdlr00(mgr), dqpdlrp1(mgr), dqpdlu00(mgr), dqpdlup1(mgr),
     .       dqpdld00(mgr), dqpdldp1(mgr), dqpdlt00(mgr), 
     .
     .       dqxdlr00(mgr), dqxdlrp1(mgr), dqxdlu00(mgr), dqxdlup1(mgr),
     .       dqxdld00(mgr), dqxdldp1(mgr), dqxdlt00(mgr),
     .
     .       dqmdlrm1(mgr), dqmdlr00(mgr), dqmdlrp1(mgr), dqmdlx00(mgr),
     .       dqmdlc00(mgr), dqmdlu00(mgr), dqmdlup1(mgr),
     .       dqmdldm1(mgr), dqmdld00(mgr), dqmdltm1(mgr), dqmdlt00(mgr)
c
      common /viscos2/ 
     .       qu      (mgr), qv      (mgr), qh      (mgr), pk      (mgr),
     .       qf      (mgr),
     .
     .       dqudlr00(mgr), dqudlrp1(mgr), dqudlu00(mgr), dqudlup1(mgr),
     .
     .       dqvdlrm1(mgr), dqvdlr00(mgr), dqvdlrp1(mgr), dqvdlx00(mgr),
     .       dqvdlc00(mgr), dqvdlu00(mgr), dqvdlup1(mgr), 
     .       dqvdldm1(mgr), dqvdld00(mgr), dqvdltm1(mgr), dqvdlt00(mgr),
     .
     .       dqhdlr00(mgr), dqhdlrp1(mgr), 
     .
     .       dpkdld00(mgr), dpkdlu00(mgr), dpkdlup1(mgr), dpkdlt00(mgr),
     .
     .       dqfdld00(mgr), dqfdlu00(mgr), dqfdlup1(mgr), dqfdlt00(mgr),
     .                      dqfdlr00(mgr), dqfdlrp1(mgr)
c
c--------------------------------------------------------------titan.com
c
