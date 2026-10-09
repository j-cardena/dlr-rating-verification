import numpy as np
from scipy.stats import chi2, binom
def lr_uc(N,T,a):
    p=N/T
    def ll(q):
        return (T-N)*np.log1p(-q)+(N*np.log(q) if N>0 else 0)
    return -2*(ll(a)-(ll(p) if 0<p<1 else 0))
crit=chi2.ppf(0.95,1)
for T,a,label in [(8760,0.01,'1 yr hourly'),(8760,0.05,'1 yr hourly'),(8760,0.02,'1 yr hourly'),(365,0.01,'1 yr daily'),(365,0.05,'1 yr daily'),(2190,0.01,'3 months hourly'),(17520,0.01,'2 yr hourly')]:
    acc=[N for N in range(0,T+1) if lr_uc(N,T,a)<=crit]
    lo,hi=min(acc),max(acc)
    print(f"{label:16s} T={T:6d} a={a:.2f} E[N]={a*T:7.2f} accept N in [{lo},{hi}] -> rate [{lo/T*100:.2f}%,{hi/T*100:.2f}%]")
# power: true rate 2a vs a, T needed for 80% power
for a in [0.01,0.02,0.05]:
    for T in [365,730,2190,4380,8760,17520]:
        acc=[N for N in range(0,min(T,int(10*a*T)+60)+1) if lr_uc(N,T,a)<=crit]
        lo,hi=min(acc),max(acc)
        pw=1-(binom.cdf(hi,T,2*a)-binom.cdf(lo-1,T,2*a))
        print(f"a={a} T={T} power vs 2a = {pw:.2f}")
# effective sample size under persistence: if exceedances come in runs of mean length L, effective independent events ~ N/L
