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

null=source(SRC/'Fig7a_matched_null_k50_v8_1.csv')
closest=source(SRC/'Fig7a_closest10_sensitivity_v8_1.csv')
balance=source(SRC/'Fig7b_matching_balance_null_v8_1.csv')
frozen=source(SRC/'Fig7b_matching_balance_frozen_v8_1.csv')
patients=source(SRC/'Fig7d_gse132465_patient_scores_v8_1.csv')
fig,axs=plt.subplots(2,2,figsize=(12.8,9.3))
fig.subplots_adjust(left=.085,right=.965,bottom=.085,top=.925,hspace=.51,wspace=.32)
ax=axs[0,0]; ax.set_title('A  Matched-null calibration in GSE178341',loc='left',pad=30)
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
ax=axs[0,1]; ax.set_title('B  Matching balance',loc='left',pad=30)
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
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
title_boxes=[ax._left_title.get_window_extent(renderer) for ax in axs[0]]
assert abs(title_boxes[0].y0-title_boxes[1].y0)<1e-6
save(fig,'Figure_1')
(OUT/'figure1_title_alignment.json').write_text(json.dumps({'change':'Panel A title raised to align with panel B; panel B spacing unchanged','title_y0_pixels':[b.y0 for b in title_boxes],'status':'regenerated copy; archived approved exports remain unchanged','exports':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'png/Figure_1.png',OUT/'tiff/Figure_1.tiff',OUT/'pdf/Figure_1.pdf']}},indent=2),encoding='utf8')

(OUT/'source_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf8')
