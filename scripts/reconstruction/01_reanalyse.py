"""Targeted post-review repair. Inputs read only; all outputs under v09."""
from pathlib import Path
import json, hashlib, platform, importlib.util, time, os
import numpy as np
import pandas as pd
import scipy
from scipy import stats, linalg

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('TROGO_OUTPUT',str(ROOT/'qa/recomputed/reconstruction_v09')));OUT.mkdir(parents=True,exist_ok=True)
REPO=ROOT
RAW=Path(os.environ['TROGO_RAW_INPUTS'])
TARGETS=['SCD','CH25H','ATF3','FADS2']
INPUTS={'selected':REPO/'baseline/v7/analysis/organoid/organoid_selected_libraries.csv',
 'lfc':RAW/'organoid/supplementary_table_6_revision.csv',
 'model':RAW/'depmap_26Q1/Model.csv','effect':RAW/'depmap_26Q1/CRISPRGeneEffect.csv'}
def hashfile(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def savejson(name,x): (OUT/name).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf8')
savejson('run_inputs.json',{'inputs':{k:{'path':str(p),'sha256':hashfile(p),'bytes':p.stat().st_size} for k,p in INPUTS.items()},'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'script_sha256':hashfile(Path(__file__))})
print('Input hashes recorded',flush=True)

# Organoid complete matrix and original QC design.
sel=pd.read_csv(INPUTS['selected']);assert len(sel)==85 and sel.sample_ID.is_unique
samples=sel.sample_ID.astype(str).tolist();sidx={s:i for i,s in enumerate(samples)}
pairs=set(zip(sel.sample_ID.astype(str),sel.library.astype(str)))
values={}; duplicates=0
for chunk in pd.read_csv(INPUTS['lfc'],sep=r'\s+',chunksize=250000):
 mask=[(str(s),str(l)) in pairs for s,l in zip(chunk.sample_ID,chunk.library)]
 for s,g,v in chunk.loc[mask,['sample_ID','gene','LFC']].itertuples(index=False):
  d=values.setdefault(str(g),{}); j=sidx[str(s)]
  if j in d: duplicates+=1; assert d[j]==float(v)
  d[j]=float(v)
genes=sorted(g for g,d in values.items() if len(d)==85 and all(np.isfinite(list(d.values()))))
M=np.array([[values[g][j] for j in range(85)] for g in genes]);assert M.shape==(16562,85)
C=np.column_stack([np.ones(85),stats.rankdata(np.median(M,axis=0)),stats.rankdata(sel.AUC_ROC),stats.rankdata(sel.AUC_PR)])
assert np.isfinite(C).all();rank=int(np.linalg.matrix_rank(C));df=85-rank-1;assert df==80
R=stats.rankdata(M,axis=1);Rc=R-R.mean(axis=1,keepdims=True)
res=R.T-C@np.linalg.lstsq(C,R.T,rcond=None)[0];res-=res.mean(axis=0)
rn=np.linalg.norm(Rc,axis=1);en=np.linalg.norm(res,axis=0)
all_tables=[];summaries=[]
for g in TARGETS:
 j=genes.index(g); raw=np.clip(Rc@Rc[j]/(rn*rn[j]),-1,1);adj=np.clip(res.T@res[:,j]/(en*en[j]),-1,1)
 keep=np.arange(len(genes))!=j
 def pval(r,d):return 2*stats.t.sf(np.abs(r)*np.sqrt(d/np.maximum(1-r*r,1e-15)),d)
 tab=pd.DataFrame({'anchor':g,'gene':np.array(genes)[keep],'rho_raw':raw[keep],'rho_adj':adj[keep], 'p_raw':pval(raw[keep],83),'p_adj':pval(adj[keep],df)})
 tab['q_raw']=stats.false_discovery_control(tab.p_raw);tab['q_adj']=stats.false_discovery_control(tab.p_adj)
 assert len(tab)==16561 and tab[['p_adj','q_adj']].notna().all().all()
 # Independent algebra check: target added to reduced rank regression.
 hit='VPS72' if g=='SCD' else 'SMYD2' if g=='ATF3' else 'SCD'
 y=R[genes.index(hit)];X=np.column_stack([C,R[j]]);b=np.linalg.lstsq(X,y,rcond=None)[0]
 e=y-X@b;se=np.sqrt((e@e/df)*np.linalg.inv(X.T@X)[-1,-1]);regp=2*stats.t.sf(abs(b[-1]/se),df)
 assert abs(regp-tab.set_index('gene').loc[hit,'p_adj'])<1e-10
 tab.sort_values('rho_adj',ascending=False).to_csv(OUT/f'organoid_{g}.csv',index=False)
 all_tables.append(tab)
 summaries.append({'anchor':g,'n':85,'df':df,'tests':len(tab),'median_LFC':float(np.median(M[j])),'proportion_LFC_below_minus1':float(np.mean(M[j]<-1)), 'significant':tab[tab.q_adj<.05].to_dict('records')})
 print(g,summaries[-1],flush=True)
pool=pd.concat(all_tables,ignore_index=True);pool['q_pooled_4anchors']=stats.false_discovery_control(pool.p_adj)
pool.to_csv(OUT/'organoid_all_anchors.csv.gz',index=False)
pd.DataFrame(M.T,index=samples,columns=genes).to_csv(OUT/'organoid_complete_matrix.csv.gz')
pd.DataFrame(C[:,1:],index=samples,columns=['rank_global_median','rank_AUC_ROC','rank_AUC_PR']).to_csv(OUT/'organoid_QC_design.csv')
savejson('organoid_audit.json',{'duplicate_pairs':duplicates,'design_rank':rank,'df':df,'target_summaries':summaries})

# Match original DepMap representation (float32 input, float64 calculations).
model=pd.read_csv(INPUTS['model'],low_memory=False);cell=model[model.ModelType.eq('Cell Line')]
crc_ids=set(cell.loc[cell.OncotreePrimaryDisease.eq('Colorectal Adenocarcinoma'),'ModelID'])
parts=[]
for c in pd.read_csv(INPUTS['effect'],index_col=0,chunksize=128):
 parts.append(c.loc[c.index.isin(set(cell.ModelID))].astype(np.float32))
ge=pd.concat(parts).astype(float);ge.columns=[x.split(' (')[0] for x in ge.columns]
ge=ge.loc[:,~ge.columns.duplicated(keep=False)];assert ge.shape==(1208,18531) and ge.index.is_unique
crc=ge.index.isin(crc_ids);assert sum(crc)==63
lineage=[];rng=np.random.default_rng(42)
for g in TARGETS:
 a=ge.loc[crc,g].dropna().to_numpy();b=ge.loc[~crc,g].dropna().to_numpy()
 boot=np.array([np.median(rng.choice(a,len(a)))-np.median(rng.choice(b,len(b))) for _ in range(10000)])
 lineage.append({'gene':g,'crc_n':len(a),'noncrc_n':len(b),'crc_median':np.median(a),'noncrc_median':np.median(b),'difference':np.median(a)-np.median(b),'ci_low':np.quantile(boot,.025),'ci_high':np.quantile(boot,.975),'p':stats.mannwhitneyu(a,b).pvalue,'crc_below_minus1':int(sum(a<-1)),'rank':int(ge.loc[crc].median().rank(method='min')[g])})
lt=pd.DataFrame(lineage);lt['q_4anchors']=stats.false_discovery_control(lt.p);lt.to_csv(OUT/'lineage_four_anchors.csv',index=False)
print('Lineage checks complete',flush=True)

# Import pure functions only from archived implementation.
sp=importlib.util.spec_from_file_location('original_depmap',REPO/'scripts/analysis/02_depmap_scd_vps72_replication.py');mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
mat=ge.to_numpy();completeness=np.mean(np.isfinite(mat),axis=0);med=np.nanmedian(mat,axis=0);dep=np.mean(mat<-.5,axis=0)
rowmed=np.nanmedian(mat,axis=1);mad=np.nanmedian(abs(mat-rowmed[:,None]),axis=1)
replication=[]
for a,b in [('SCD','VPS72'),('ATF3','SMYD2')]:
 rng=np.random.default_rng(20260901)
 mask=(completeness>=.9)&(med<=-.5)&(dep>=.8)&~ge.columns.isin([a,b])
 common=mat[:,mask].copy();mr,mc=np.where(~np.isfinite(common));common[mr,mc]=np.nanmedian(common,axis=0)[mc]
 sd=common.std(axis=0,ddof=1);variable=sd>1e-8; common=(common[:,variable]-common[:,variable].mean(axis=0))/sd[variable]
 u,s,_=linalg.svd(common,full_matrices=False);cov=np.column_stack([rowmed,mad,u[:,:5]*s[:5]])
 good=crc&np.isfinite(ge[a])&np.isfinite(ge[b])&np.isfinite(cov).all(axis=1)
 x=ge.loc[good,a].to_numpy();y=ge.loc[good,b].to_numpy();cc=cov[good]
 rho=mod.residual_correlation(x,y,cc);xr,yr=mod.adjusted_residuals(x,y,cc)
 perm=np.array([np.corrcoef(xr,rng.permutation(yr))[0,1] for _ in range(10000)])
 pp=(1+sum(abs(perm)>=abs(rho)))/10001
 boot=mod.bootstrap(x,y,cc,10000,rng)
 boot.to_csv(OUT/f'depmap_{a}_{b}_bootstrap.csv.gz',index=False)
 pd.DataFrame({'rho':perm}).to_csv(OUT/f'depmap_{a}_{b}_permutation.csv.gz',index=False)
 tab=pd.DataFrame({'ModelID':ge.index[good],a:x,b:y,**{f'cov{i}':cc[:,i] for i in range(7)}});tab.to_csv(OUT/f'depmap_{a}_{b}_source.csv',index=False)
 ci=np.quantile(boot.adjusted_partial_spearman,[.025,.975]);raw=stats.spearmanr(x,y)
 rr={'anchor':a,'partner':b,'n':len(x),'common_essential_n':int(sum(mask)),'rho_raw':raw.statistic,'p_raw':raw.pvalue,'rho_adj':rho,'ci_low':ci[0],'ci_high':ci[1],'p_perm':pp}
 replication.append(rr);print(rr,flush=True)
 if a=='SCD':assert abs(rho-(-.175))<.001
pd.DataFrame(replication).to_csv(OUT/'depmap_replication.csv',index=False)
savejson('run_completion.json',{'status':'completed','scope':'targeted organoid statistics, four-target lineage and two-pair DepMap checks; not full raw single-cell rerun','time_local':time.strftime('%Y-%m-%dT%H:%M:%S'),'script_sha256':hashfile(Path(__file__))})
