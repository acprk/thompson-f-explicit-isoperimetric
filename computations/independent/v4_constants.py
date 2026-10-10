import sys; sys.path.insert(0,'.')
from fractions import Fraction as Fr
from v3_avg import closed
def best(kappa, ms):
    out=None
    for m in ms:
        D=2**m
        if kappa<=Fr(4,D): continue
        Nb,e0,M=closed(m); c=(kappa-Fr(4,D))/(4*(1-Fr(1,D))*M)
        if out is None or c>out[0]: out=(c,m,M)
    return out
def rep(name,k,ms):
    c,m,M=best(k,ms); print(f"{name}: kappa={float(k):.6e} best m={m} Mbar={M} ({float(M):.4f}) h_E>={float(c):.6e}")
    return c
# self-contained (Lean formula): rho=1/2048
rho=Fr(1,2048); LipG=max(Fr(1),30+Fr(7)/(4*rho)+3808*rho,5+Fr(7)/(8*rho)); L=8*LipG/rho
print("LipG",LipG,float(LipG)," L",L)
cP=rep("self-contained (delta=1/2, L=8LipG/rho)",(Fr(1,2)/L)**2,range(30,90))
# Hilbert normalization is (1/eps)-Lipschitz, not 2/eps -> L/2
rep("  ... with 1/eps normalization (L/2)",(Fr(1,2)/(L/2))**2,range(30,90))
for k in (29,33,Fr(3226,100)):
    d=Fr(2)/(2+Fr(k)); rep(f"retraction k={float(k)}: delta=2/(k+2), L=k",(d/Fr(k))**2,range(8,40))
# paper's own S (max over transports): kappa/4 sup
print("paper-S sup kappa/4 for k=29:", float(((Fr(2,31)/29)**2)/4))
