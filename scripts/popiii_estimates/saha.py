# Equilibrium (Saha) chemistry of primordial gas and its adiabat: H2 <-> 2H, H <-> p + e, H- ; He neutral spectator.
# Ideal gas, rigid-rotor + harmonic-oscillator H2 (standard molecular constants, to be replaced by verified data).
import numpy as np
kB=1.380649e-16; h=6.62607015e-27; me=9.1093837e-28; mH=1.6735575e-24; eV=1.602176634e-12; G=6.674e-8
X=0.76; Y=0.24
D0=4.478*eV; chi=13.598*eV; chim=0.754*eV; thr=85.3; thv=6332.
def Zint(T):
    J=np.arange(0,40); g=np.where(J%2==0,1,3)
    zr=np.sum((2*J+1)*g*np.exp(-thr*J*(J+1)/T)); zv=1/(1-np.exp(-thv/T))
    # mean internal energy (rot+vib) per molecule
    er=np.sum((2*J+1)*g*thr*J*(J+1)*np.exp(-thr*J*(J+1)/T))/zr*kB
    ev=kB*thv/(np.exp(thv/T)-1)
    return zr*zv, er+ev
def lam(m,T): return (2*np.pi*m*kB*T/h**2)**1.5
def comp(rho,T):
    nHt=X*rho/mH; nHe=Y*rho/(4*mH)
    Z,_=Zint(T)
    Kd=(lam(mH,T)*4)**2/(lam(2*mH,T)*Z)*np.exp(-D0/(kB*T))   # nH^2/nH2
    Ki=lam(me,T)*np.exp(-chi/(kB*T))                          # np ne/nH
    f=lambda a: a+2*a*a/Kd+np.sqrt(a*Ki)-nHt
    lo,hi=0.0,nHt
    for _ in range(200):
        m=0.5*(lo+hi)
        if f(m)>0: hi=m
        else: lo=m
    a=0.5*(lo+hi); nH2=a*a/Kd; npr=np.sqrt(a*Ki); ne=npr
    nHm=a*ne/(4*lam(me,T)*np.exp(-chim/(kB*T)))
    return dict(nH=a,nH2=nH2,np=npr,ne=ne,nHm=nHm,nHe=nHe,nHt=nHt)
def state(rho,T):
    c=comp(rho,T); _,eint=Zint(T)
    n=c['nH']+c['nH2']+c['np']+c['ne']+c['nHe']
    P=n*kB*T
    u=(1.5*n*kB*T + c['nH2']*eint + c['nH']*D0/2 + c['np']*(D0/2+chi))/rho
    return P,u,c
if __name__=="__main__":
    print("# half-dissociation T (x_H2 = 2nH2/nHt = 0.5) and half-ionisation T vs n_H(tot)")
    for lg in [10,12,14,16,18,20,22,24]:
        rho=10**lg*mH/X
        Ts=np.logspace(2.8,5.5,4000); xm=[2*comp(rho,T)['nH2']/comp(rho,T)['nHt'] for T in Ts]
        xi=[comp(rho,T)['np']/comp(rho,T)['nHt'] for T in Ts]
        Td=Ts[np.argmin(np.abs(np.array(xm)-0.5))]; Ti=Ts[np.argmin(np.abs(np.array(xi)-0.5))]
        print(f"n_H=1e{lg:2d}  rho={rho:9.2e}  T_half_diss={Td:7.0f} K  T_half_ion={Ti:8.0f} K")
