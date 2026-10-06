"""Source-driven reconstruction figures. No inference rerun or generated imagery."""
from pathlib import Path
import json, hashlib, shutil, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT
SRC=REPO/'figures/source_data'
OUT=Path(os.environ.get('TROGO_FIGURE_OUTPUT',str(ROOT/'qa/recomputed/reconstruction_figures')))
for folder in ['png','tiff','pdf','source_data']: (OUT/folder).mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.titlesize':11,'pdf.fonttype':42,'savefig.facecolor':'white'})
blue='#365F91'; red='#C44333'; pale='#BCD2E3'; grey='#75828B'
records=[]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def source(p):
    dest=OUT/'source_data'/p.name
    shutil.copy2(p,dest)
    records.append({'source':str(p),'copy':str(dest.relative_to(ROOT)),'sha256':sha(p)})
    return pd.read_csv(p)
def save(fig,name):
    for fmt,folder in [('png','png'),('tiff','tiff'),('pdf','pdf')]:
        kw={'pil_kwargs':{'compression':'tiff_lzw'}} if fmt=='tiff' else {}
        fig.savefig(OUT/folder/f'{name}.{fmt}',dpi=300,**kw)
    plt.close(fig)
def title(ax,label,text): ax.set_title(f'{label}  {text}',loc='left',pad=12)

# Five previously approved figures are preserved byte-for-byte.
for new,old in [('Figure_2',5),('Figure_3',6),('Figure_4',4),('Figure_S1',1),('Figure_S2',2)]:
    for fmt,folder in [('png','png'),('tiff','tiff'),('pdf','pdf')]:
        p=REPO/f'figures/frontiers_v8_6/Figure_{old}.{fmt}'
        d=OUT/folder/f'{new}.{fmt}'
        shutil.copy2(p,d)
        assert sha(p)==sha(d)
        records.append({'source':str(p),'copy':str(d.relative_to(ROOT)),'sha256':sha(d),'mode':'unchanged approved figure'})

null=source(SRC/'Fig7a_matched_null_k50_v8_1.csv')
closest=source(SRC/'Fig7a_closest10_sensitivity_v8_1.csv')
balance=source(SRC/'Fig7b_matching_balance_null_v8_1.csv')
frozen=source(SRC/'Fig7b_matching_balance_frozen_v8_1.csv')
patients=source(SRC/'Fig7d_gse132465_patient_scores_v8_1.csv')
fig,axs=plt.subplots(2,2,figsize=(12.8,9.3))
fig.subplots_adjust(left=.085,right=.965,bottom=.085,top=.94,hspace=.51,wspace=.32)
ax=axs[0,0]; title(ax,'A','Matched-null calibration in GSE178341')
xx=np.linspace(.30,.95,600)
yy=gaussian_kde(null.rho)(xx)
ax.fill_between(xx,yy,color=pale,alpha=.8,label='Primary matched null (10,000)')
ax.plot(xx,yy,color=blue,lw=1.7)
ax.plot(xx,gaussian_kde(closest.rho)(xx),color=red,lw=1.8,label='Closest 10% (post hoc)')
observed=json.loads((REPO/'analysis/matched_null_gse178341/gse178341_matched_null_audit.json').read_text())['observed_rho_complete_gene_cpm']
ax.axvline(observed,color=red,lw=1.3,ls='--')
ax.set(xlim=(.30,.95),ylim=(0,yy.max()*1.2),xlabel='Spearman ρ with T-cell fraction',ylabel='Density')
ax.text(.025,.93,'Frozen ρ = 0.767\nNull median = 0.766\nOne-sided P = 0.493\nEqual-tail two-sided P = 0.985',transform=ax.transAxes,va='top',fontsize=9.5,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
ax.legend(loc='lower left',fontsize=8.5,frameon=False)
ax=axs[0,1]; title(ax,'B','Matching balance')
features=list(frozen.feature)
vals=[balance.loc[balance.feature==f,'standardized_null_value'].values for f in features]
v=ax.violinplot(vals,positions=[1,2,3],showextrema=False,widths=.68)
for b in v['bodies']: b.set_facecolor(pale); b.set_edgecolor(blue); b.set_alpha(.9)
ax.axhline(0,color='#CCD2D8',lw=1)
for i,r in frozen.iterrows():
    ax.scatter(i+1,r.standardized_null_value,marker='D',s=40,color=red,zorder=4)
    ax.text(i+1,4.35,f'{r.percentile:.1f}%',ha='center',va='center',color=red,fontweight='bold',fontsize=10)
ax.set(xticks=[1,2,3],xticklabels=['Mean\nexpression','Cell\ndetection','TNK/ILC\nspecificity'],ylabel='Standardised matched-null feature',ylim=(-4.4,4.85))
ax.text(.5,1.015,'Red labels show frozen-module null percentiles',transform=ax.transAxes,ha='center',fontsize=8.7,color=grey)
stats=[('0.735','0.392 to 0.926','0.00060'),('0.376','−0.127 to 0.784','0.0892')]
for ax,comp,letter,st in zip(axs[1],['Tumour-wide pseudobulk','Epithelial pseudobulk'],['C','D'],stats):
    d=patients[patients.compartment==comp]
    assert len(d)==23,(comp,len(d),patients.compartment.unique())
    title(ax,letter,f'GSE132465 {comp.lower()}')
    ax.scatter(d.t_cell_fraction,d.score,s=35,color=blue,alpha=.85,edgecolor='white',linewidth=.4)
    x=np.array([d.t_cell_fraction.min(),d.t_cell_fraction.max()]); coef=np.polyfit(d.t_cell_fraction,d.score,1)
    ax.plot(x,np.polyval(coef,x),color=blue,lw=1.6)
    ax.set(xlabel='Author-annotated T-cell fraction',ylabel='Eight-gene RNA module score')
    yl=ax.get_ylim(); ax.set_ylim(yl[0],yl[1]+.42*(yl[1]-yl[0]))
    ax.text(.025,.975,f'n = 23 patients; ρ = {st[0]}\nBootstrap 95% CI {st[1]}\nPermutation P = {st[2]}',transform=ax.transAxes,va='top',fontsize=9.2,bbox={'facecolor':'white','edgecolor':'none','alpha':.94})
fig.text(.5,.013,'Absolute-correlation attenuation = 0.359 (bootstrap 95% CI −0.041 to 0.759)',ha='center',fontsize=9.5)
save(fig,'Figure_1')

lineage=source(ROOT/'analysis/reconstruction_v09/lineage_four_anchors.csv')
rep=source(ROOT/'analysis/reconstruction_v09/depmap_replication.csv')
fig,axs=plt.subplots(1,2,figsize=(13.5,5.4),gridspec_kw={'width_ratios':[1,1.13]})
fig.subplots_adjust(left=.075,right=.97,bottom=.23,top=.86,wspace=.40)
ax=axs[0]; title(ax,'A','CRC versus other cell-line lineages')
for i,r in lineage.iterrows():
    y=3-i
    ax.errorbar(r.difference,y,xerr=[[r.difference-r.ci_low],[r.ci_high-r.difference]],fmt='o',color=blue,capsize=4)
    ax.text(.095,y+.15,f'Δ = {r.difference:.3f}\n95% CI {r.ci_low:.3f}, {r.ci_high:.3f}\nP = {r.p:.3g}; q = {r.q_4anchors:.3g}',va='center',fontsize=8.5)
ax.axvline(0,color=grey,ls='--',lw=1)
ax.set(yticks=[3,2,1,0],yticklabels=lineage.gene,xlim=(-.90,.60),ylim=(-.60,3.7),xlabel='CRC minus non-CRC median Chronos effect')
ax=axs[1]; title(ax,'B','Adjusted DepMap CRC association checks')
for i,r in rep.iterrows():
    y=1-i
    ax.errorbar(r.rho_adj,y,xerr=[[r.rho_adj-r.ci_low],[r.ci_high-r.rho_adj]],fmt='o',color=blue,capsize=4)
    ax.text(-.59,y+.26,f'{r.anchor}–{r.partner}: ρ = {r.rho_adj:.3f}\n95% CI {r.ci_low:.3f}, {r.ci_high:.3f}; permutation P = {r.p_perm:.3f}',fontsize=9,va='bottom',bbox={'facecolor':'white','edgecolor':'none','pad':1})
ax.axvline(0,color=grey,ls='--',lw=1)
ax.set(xlim=(-.62,.62),ylim=(-.35,1.85),yticks=[],xlabel='QC/common-essentiality-adjusted rank correlation')
fig.text(.075,.065,'63 CRC and 1,145 non-CRC cell lines\nFour-target BH correction; 10,000 bootstrap resamples',fontsize=9,color=grey)
fig.text(.565,.065,'63 CRC cell lines; 1,179 common-essential genes\n10,000 refitted bootstraps and 10,000 residual permutations',fontsize=9,color=grey)
save(fig,'Figure_S3')

allcorr=source(ROOT/'analysis/reconstruction_v09/organoid_all_anchors.csv.gz')
gsea=source(REPO/'baseline/v7/analysis/scd_codependency_qc/scd_hallmark_raw_vs_qc.csv')
fig=plt.figure(figsize=(13.2,10.6)); grid=fig.add_gridspec(3,2,height_ratios=[1,1,1.0],left=.085,right=.97,bottom=.075,top=.95,hspace=.58,wspace=.28)
for k,anchor in enumerate(['SCD','CH25H','ATF3','FADS2']):
    ax=fig.add_subplot(grid[k//2,k%2]); title(ax,chr(65+k),f'{anchor} adjusted organoid associations')
    d=allcorr[allcorr.anchor==anchor]
    ax.scatter(d.rho_adj,-np.log10(d.q_adj),s=3,color=grey,alpha=.42,rasterized=True)
    ax.axhline(-np.log10(.05),color=red,ls='--',lw=1)
    sig=d[d.q_adj<.05]
    ax.set(xlim=(-.65,.65),ylim=(-.04,2.9),xlabel='Adjusted rank correlation',ylabel='−log10(single-anchor q)')
    for _,r in sig.iterrows():
        ax.scatter(r.rho_adj,-np.log10(r.q_adj),s=35,color=red,zorder=5)
        ax.annotate(f'{r.gene}: ρ = {r.rho_adj:.3f}\nP = {r.p_adj:.3g}; q = {r.q_adj:.4f}\nFour-anchor q = {r.q_pooled_4anchors:.4f}',xy=(r.rho_adj,-np.log10(r.q_adj)),xytext=(.42,.83),textcoords='axes fraction',fontsize=8.5,va='top',arrowprops={'arrowstyle':'-','color':red,'lw':.8})
    if sig.empty: ax.text(.97,.92,'No single-anchor q < 0.05',transform=ax.transAxes,ha='right',fontsize=9)
ax=fig.add_subplot(grid[2,:]); title(ax,'E','SCD Hallmark enrichment before and after QC adjustment')
sets=['HALLMARK_DNA_REPAIR','HALLMARK_CHOLESTEROL_HOMEOSTASIS']; labels=['DNA repair','Cholesterol\nhomeostasis']
for j,(key,color,offset) in enumerate([('raw_spearman',pale,-.17),('qc_adjusted',blue,.17)]):
    sub=gsea[gsea.ranking==key].set_index('pathway')
    assert len(sub)>0,gsea.ranking.unique()
    for i,path in enumerate(sets):
        r=sub.loc[path]; y=1-i+offset
        ax.barh(y,r.NES,height=.26,color=color,label='Raw ranks' if j==0 and i==0 else ('QC-adjusted ranks' if j==1 and i==0 else None))
        ax.text(.07,y,f'NES {r.NES:.3f}; P = {r.p_value:.4g}; q = {r.FDR_BH:.3f}',va='center',fontsize=9)
ax.axvline(0,color=grey,lw=1)
ax.set(yticks=[1,0],yticklabels=labels,xlim=(-1.85,1.5),ylim=(-.55,1.6),xlabel='Normalised enrichment score (NES)')
ax.legend(loc='upper right',frameon=False,fontsize=9)
fig.text(.5,.014,'A–D: 85 organoids; three QC covariates; df = 80; 16,561 tests per anchor. E: neither full Hallmark ranking has q < 0.05.',ha='center',fontsize=9)
save(fig,'Figure_S4')

(OUT/'provenance.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print('Completed 8 figures in PNG, TIFF and PDF; 5 figures reused unchanged.')
