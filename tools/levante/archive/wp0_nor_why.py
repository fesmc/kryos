import numpy as np, netCDF4 as nc
from scipy import ndimage as ndi
G="/work/ba1442/robinson/kryos-gate/ref/nor"
r0=f"{G}/yx-918d996/v24/cont"; r2=f"{G}/yx-5ac89c8/S4/cont"
src=nc.Dataset("/work/ba1442/fesmdata/North/NH-32KM/NH-32KM_TOPO-GEBCO2026.nc")
H=src["H_ice"][:].filled(np.nan); H=np.squeeze(H)
print("source H_ice shape",H.shape,"cells: 0<H<0.1",int(((H>0)&(H<.1)).sum()),"0.1<=H<1",int(((H>=.1)&(H<1)).sum()),"1<=H<10",int(((H>=1)&(H<10)).sum()),">=10",int((H>=10).sum()))
def f(d,v,k=0): return nc.Dataset(f"{d}/yelmo.nc")[v][k].filled(np.nan)
for v in ("H_ice","H_ice_ref"):
    a,b=f(r0,v),f(r2,v); d=np.abs(a-b); m=d>0
    print(v,"diff cells",int(m.sum()),"max",np.nanmax(d),">1m",int((d>1).sum()),">10m",int((d>10).sum()),"R0 sum",np.nansum(a),"R2 sum",np.nansum(b))
    thin=(H>=.1)&(H<1); thin2=(H>=1)&(H<10); anythin=(H>0)&(H<10)
    for lab,t in (("0.1<=H<1",thin),("1<=H<10",thin2),("0<H<10",anythin)):
        for r in (0,2,6):
            dil=ndi.binary_dilation(t,iterations=r) if r else t
            print(f"   {v}: diff cells within {r} cells of source {lab}: {int((m&dil).sum())} of {int(m.sum())}")
# direct comparison: R0-style cleaned initial ice vs R2'-style
R0i=np.where(H<0.1,0,np.where(H<10,10,H)); R2i=np.where(H<1,0,H)
dd=np.abs(R0i-R2i); print("initial ice, R0 init-path vs hub-prepared: diff cells",int((dd>0).sum()),"max",np.nanmax(dd))
hr=f(r2,"H_ice_ref"); hr0=f(r0,"H_ice_ref")
m=np.abs(hr-hr0)>0
print("diff H_ice_ref: both zero-ice cells(R0,R2)", int((m&(hr0==0)&(hr>0)).sum()), "R0 ice, R2 none", int((m&(hr0>0)&(hr==0)).sum()), "both ice", int((m&(hr0>0)&(hr>0)).sum()))
print("H_ice_ref sign: R0>R2", int(((hr0-hr)>0).sum()), "R0<R2", int(((hr0-hr)<0).sum()))
a=f(r0,"H_ice");b=f(r2,"H_ice"); m2=np.abs(a-b)>0
print("H_ice t0 diff cells also in H_ice_ref-diff dilated by 6:",int((m2&ndi.binary_dilation(m,iterations=6)).sum()),"; by 30:",int((m2&ndi.binary_dilation(m,iterations=30)).sum()),"of",int(m2.sum()))
print("H_ice t0 diff where R0 ice>0 & R2 ice>0",int((m2&(a>0)&(b>0)).sum()), "R0 ice only",int((m2&(a>0)&(b==0)).sum()),"R2 ice only",int((m2&(a==0)&(b>0)).sum()))
print("H_ice t0 diff >10 m cells:",int((np.abs(a-b)>10).sum()))
