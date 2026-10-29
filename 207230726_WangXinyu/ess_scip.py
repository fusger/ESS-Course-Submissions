"""Solve the exact YALMIP-exported MISOCP using SCIP without deleting constraints."""
from pathlib import Path
import sys, json, time
import numpy as np
from scipy.io import loadmat, savemat
from pyscipopt import Model, quicksum
src,dst=map(Path,sys.argv[1:3]);p=loadmat(src,simplify_cells=True)
F=p['F_struc'].tocsr();c=np.atleast_1d(p['c']).flatten()
lb=np.atleast_1d(p['lb']).flatten();ub=np.atleast_1d(p['ub']).flatten()
bins=set(np.atleast_1d(p['bins']).astype(int).flatten()-1)
m=Model('ESS_207230726_MISOCP')
m.setLogfile(str(dst.with_suffix('.log')))
m.setRealParam('limits/time',600);m.setRealParam('limits/gap',1e-6)
m.setRealParam('limits/absgap',1e-7);m.setRealParam('numerics/feastol',1e-9)
m.setIntParam('display/verblevel',4)
x=[m.addVar(f'x{i}',vtype='B' if i in bins else 'C',
 lb=None if not np.isfinite(lb[i]) else float(lb[i]),
 ub=None if not np.isfinite(ub[i]) else float(ub[i])) for i in range(len(c))]
def expr(row):
 r=F.getrow(row)
 return quicksum(float(v)*x[int(j)-1] for j,v in zip(r.indices,r.data) if j>0)+float(F[row,0])
nf=int(p['nf']);nl=int(p['nl'])
for r in range(nf):m.addCons(expr(r)==0)
for r in range(nf,nf+nl):m.addCons(expr(r)>=0)
top=nf+nl
for q in np.atleast_1d(p['q']).astype(int).flatten():
 if q==0:continue
 a=expr(top);m.addCons(a>=0)
 m.addCons(quicksum(expr(j)*expr(j) for j in range(top+1,top+q))<=a*a)
 top+=q
scale=1000.0
offset=float(p.get('objective_offset',0.0))
m.setObjective(scale*(quicksum(float(a)*v for a,v in zip(c,x) if a)+offset),'minimize')
if dst.exists():
 old=loadmat(dst,simplify_cells=True).get('x',[])
 if len(old)==len(x):
  warm=m.createSol()
  for var,val in zip(x,old):m.setSolVal(warm,var,float(val))
  m.addSol(warm)
t=time.time();m.optimize();elapsed=time.time()-t;status=str(m.getStatus())
if m.getNSols()==0:raise RuntimeError(f'No feasible solution: {status}')
sol=m.getBestSol();z=np.array([m.getSolVal(sol,v) for v in x])
info={'status':status,'seconds':elapsed,'lower_bound':m.getDualbound()/scale,
 'upper_bound':m.getPrimalbound()/scale,'gap':m.getGap(),'nodes':m.getNNodes(),
 'SCIP_version':f'{m.getMajorVersion()}.{m.getMinorVersion()}.{m.getTechVersion()}',
 'nvars':len(c),'binary_variables':len(bins),'cones':len(np.atleast_1d(p['q']))}
savemat(dst,dict(x=z,**info));dst.with_suffix('.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
print(json.dumps(info))
if status!='optimal' and m.getGap()>1e-6:raise RuntimeError('MIP gap is not closed')
