"""Independent cross-check of A007508/A080840 residuals against 2*C2*Li2(10^n) via Simpson quadrature."""
import numpy as np, json
C2=0.6601618158468696; K=2*C2
# OEIS values as of 2026-07-02 (raw A007508 responses preserved in research/raw/):
twin=[2,8,35,205,1224,8169,58980,440312,3424506,27412679,224376048,1870585220,
15834664872,135780321665,1177209242304,10304195697298,90948839353159,808675888577436,
7237516880334496]  # A007508 n=1..19 (a19 added Jun 02 2026; independent verification requested Jun 17 2026)
cousin=[1,8,41,203,1216,8144,58622,440258,3424680,27409999,224373161,1870585459,
15834656003,135779962760,1177207270204,10304191320777,90948823579815,808675898548206,
7237516937438708]  # A080840 n=1..19 (a19 added Jun 02 2026)

def Li2(X, steps=2_000_000):  # Simpson in t=ln x
    a,b=np.log(2.0),np.log(float(X))
    t=np.linspace(a,b,steps+1); f=np.exp(t)/t**2
    h=(b-a)/steps
    return float(h/3*(f[0]+f[-1]+4*f[1:-1:2].sum()+2*f[2:-2:2].sum()))
# convergence check
c1,c2=Li2(1e18,1_000_000),Li2(1e18,2_000_000)
rows=[]
for n in range(8,20):
    HL=K*Li2(10.0**n)
    row={"n":n,"HL":HL}
    if n<=len(twin):
        row["twin_dev"]=twin[n-1]-HL; row["twin_dev/sqrtHL"]=round((twin[n-1]-HL)/HL**0.5,2)
    if n<=len(cousin):
        row["cousin_dev"]=cousin[n-1]-HL; row["cousin_dev/sqrtHL"]=round((cousin[n-1]-HL)/HL**0.5,2)
    rows.append(row)
out={"Li2_convergence_relerr":abs(c1-c2)/c2,
     "note":"dev = OEIS value - 2*C2*Li2(10^n); dev/sqrtHL = Poisson-scale z proxy",
     "rows":[{k:(round(v,1) if isinstance(v,float) else v) for k,v in r.items()} for r in rows]}
HL19=K*Li2(1e19)
out["HL_1e19"]=round(HL19,1)
print(json.dumps(out,indent=1))
