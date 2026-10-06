"""Bounded audit/repair: local candidate inputs read-only; no downloads or publication.
Reason: self gene retained at zero; GMT version and BH terminology require correction.
Validation: exact sparse/dense ES agreement, observed ES comparison, fixed 5,000
gene-set permutations, full tables and multiplicity, no significance-driven reruns.
"""
from pathlib import Path
import hashlib,json,zipfile,xml.etree.ElementTree as ET,time,platform
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import os
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('TROGO_GSEA_OUTPUT',str(ROOT/'qa/recomputed/reconstruction_gsea')));OUT.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
gmt=Path(os.environ['TROGO_HALLMARK_GMT'])
assert sha(gmt)=='4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81','Use the recorded GMT; do not silently substitute another release.'
rankpath=Path(os.environ.get('TROGO_RANKING',str(ROOT/'analysis/reconstruction_v09/organoid_all_anchors.csv.gz')))
manifest={'inputs':[{'filename':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in [gmt,rankpath]],'python':platform.python_version(),'numpy':np.__version__,'seed':42,'permutations':5000,'status':'running'}
sets={}
for line in gmt.read_text(encoding='utf8').splitlines():
    parts=line.split('\t');sets[parts[0]]=set(parts[2:])
assert len(sets)==50
d=pd.read_csv(rankpath);allrows=[];checks=[];start=time.time()
def sparse_es(w,pos):
    k=len(pos);n=len(w);hit=w[pos];cum=np.cumsum(hit)/hit.sum()
    after=cum-(pos-np.arange(k))/(n-k)
    before=np.r_[0,cum[:-1]]-(pos-np.arange(k))/(n-k)
    hi=after.max();lo=before.min()
    return float(hi if abs(hi)>=abs(lo) else lo)
def dense_es(w,pos):
    hit=np.zeros(len(w),dtype=bool);hit[pos]=True
    walk=np.cumsum(w*hit/w[hit].sum()-(~hit)/(~hit).sum())
    return float(walk[np.argmax(np.abs(walk))])
for anchor in ['SCD','CH25H','ATF3','FADS2']:
    a=d[d.anchor==anchor].copy();a=a[a.gene!=anchor].dropna(subset=['rho_adj']).sort_values(['rho_adj','gene'],ascending=[False,True])
    assert len(a)==16561 and a.gene.is_unique and anchor not in set(a.gene)
    genes=a.gene.to_numpy();w=np.abs(a.rho_adj.to_numpy());n=len(a)
    rng=np.random.default_rng(42);rows=[]
    for name,gs in sets.items():
        pos=np.flatnonzero(np.isin(genes,list(gs)));k=len(pos)
        if not 15<=k<=500: continue
        obs=sparse_es(w,pos); dense=dense_es(w,pos)
        assert np.isclose(obs,dense,atol=1e-11)
        null=np.empty(5000)
        for i in range(5000):null[i]=sparse_es(w,np.sort(rng.choice(n,k,replace=False)))
        same=null[null>=0] if obs>=0 else null[null<0]
        extreme=np.sum(same>=obs) if obs>=0 else np.sum(same<=obs)
        prob=(int(extreme)+1)/(len(same)+1)
        nes=obs/np.mean(np.abs(same))
        rows.append({'anchor':anchor,'pathway':name,'size':k,'ES':obs,'NES':nes,'p_nominal':prob,'same_sign_null_n':len(same),'extreme_n':int(extreme),'ranking_n':n})
    tab=pd.DataFrame(rows);tab['q_BH_within_anchor']=false_discovery_control(tab.p_nominal)
    allrows.append(tab)
    checks.append({'anchor':anchor,'n_sets':len(tab),'minimum_BH_q':float(tab.q_BH_within_anchor.min()),'n_BH_lt_005':int((tab.q_BH_within_anchor<.05).sum()),'top':tab.nsmallest(1,'q_BH_within_anchor').to_dict('records')[0]})
    print(anchor,json.dumps(checks[-1]),'elapsed',round(time.time()-start,1),flush=True)
full=pd.concat(allrows,ignore_index=True);full['q_BH_pooled_200']=false_discovery_control(full.p_nominal)
full.to_csv(OUT/'Table_S10_Hallmark_four_anchor_exploratory.csv',index=False)
manifest.update(status='completed',elapsed_seconds=time.time()-start,sets=50,anchors=4,rows=len(full),checks=checks,sparse_dense_ES='all 200 observed scores agree within 1e-11',GMT_claim='Local Hallmark 2020-labelled GMT; no verified equivalence to 2023.2; not redistributed',correction='BH of plus-one same-sign nominal P values, not pooled-NES GSEA FDR',script_sha256=sha(Path(__file__)))
(OUT/'gsea_run.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
