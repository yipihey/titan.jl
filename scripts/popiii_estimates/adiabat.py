import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from saha import *
e=4.80320e-10
def T_at(rho,u_target,Tg):
    lo,hi=Tg*0.5,Tg*3
    for _ in range(100):
        m=np.sqrt(lo*hi); 
        if state(rho,m)[1]>u_target: hi=m
        else: lo=m
    return np.sqrt(lo*hi)
for T0 in (2000.,3000.):
    print(f"\n# equilibrium adiabat from n_H = 1e16 cm^-3, T = {T0:.0f} K (ideal gas, Saha; He neutral)")
    print("  log n_H   rho[g/cc]    T[K]   x_H2    x_e      n(H-)/n_H   Gamma1   t_ff[s]   Gamma_coul")
    lr=np.log(1e16*mH/X); T=T0; dl=0.002
    P,u,c=state(np.exp(lr),T); prev=(lr,np.log(P))
    nxt=16.0
    while lr < np.log(10.0):
        rho=np.exp(lr); P,u,c=state(rho,T)
        # midpoint step: du = P/rho^2 drho
        rho2=np.exp(lr+dl); Tg=T_at(rho2,u+P/rho**2*(rho2-rho),T)
        Pm,_,_=state(np.sqrt(rho*rho2),np.sqrt(T*Tg)); rhom=np.sqrt(rho*rho2)
        T=T_at(rho2,u+Pm/rhom**2*(rho2-rho),Tg); lr+=dl
        lgn=np.log10(np.exp(lr)*X/mH)
        if lgn>=nxt-1e-9:
            P2,_,c=state(np.exp(lr),T); P1,_,_=state(np.exp(lr-dl),T_at(np.exp(lr-dl),state(np.exp(lr),T)[1]-P2/np.exp(lr)**2*(np.exp(lr)-np.exp(lr-dl)),T))
            g1=(np.log(P2)-np.log(P1))/dl
            tff=np.sqrt(3*np.pi/(32*G*np.exp(lr)))
            ni=c['np']+c['nH']+c['nH2']+c['nHe']; a=(3/(4*np.pi*ni))**(1/3); Gc=e*e/(a*kB*T)
            print(f"  {lgn:6.2f}  {np.exp(lr):9.2e} {T:9.0f}  {2*c['nH2']/c['nHt']:6.3f}  {c['ne']/c['nHt']:8.2e}  {c['nHm']/c['nHt']:9.2e}  {g1:6.3f}  {tff:8.2e}  {Gc:7.3f}")
            nxt+=0.5
