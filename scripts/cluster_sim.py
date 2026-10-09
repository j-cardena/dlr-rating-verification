import numpy as np
from scipy.stats import chi2
rng=np.random.default_rng(1)
def ll(n0,n1,q):
    r=0.0
    if n0: r+=n0*np.log1p(-q)
    if n1: r+=n1*np.log(q)
    return r
def tests(I,a):
    T=len(I);N=I.sum();p=N/T
    luc=-2*(ll(T-N,N,a)-(ll(T-N,N,p) if 0<p<1 else 0))
    x0,x1=I[:-1],I[1:]
    n00=np.sum((x0==0)&(x1==0));n01=np.sum((x0==0)&(x1==1));n10=np.sum((x0==1)&(x1==0));n11=np.sum((x0==1)&(x1==1))
    p01=n01/max(n00+n01,1);p11=n11/max(n10+n11,1);p2=(n01+n11)/max(n00+n01+n10+n11,1)
    l1=ll(n00,n01,p01) if 0<p01<1 else 0
    l1+=ll(n10,n11,p11) if 0<p11<1 else 0
    l0=ll(n00+n10,n01+n11,p2) if 0<p2<1 else 0
    lind=-2*(l0-l1)
    return luc,lind,luc+lind
def sim(T,a,L):
    p11=1-1/L; p01=a*(1-p11)/(1-a)
    I=np.zeros(T,dtype=int); I[0]=rng.random()<a
    u=rng.random(T)
    for t in range(1,T):
        I[t]= u[t]<(p11 if I[t-1] else p01)
    return I
T=8760;a=0.01;R=2000
for L in [1/(1-0.01),2,3,6]:
    res=np.array([tests(sim(T,a,L),a) for _ in range(R)])
    r_uc=(res[:,0]>chi2.ppf(.95,1)).mean(); r_ind=(res[:,1]>chi2.ppf(.95,1)).mean(); r_cc=(res[:,2]>chi2.ppf(.95,2)).mean()
    print(f"L={L:.2f}h  calibrated a=1%: reject UC {r_uc:.2f}  IND {r_ind:.2f}  CC {r_cc:.2f}")
