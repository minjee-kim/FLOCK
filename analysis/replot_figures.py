"""Regenerate the report's figures from bundled historical CSV summaries only.
Does not run the original FLOCK simulations.
"""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'data'; F=ROOT/'figures'; F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
def graph(fn,which,title,output,xlabel):
    d=pd.read_csv(D/fn).sort_values(which)
    fig,ax=plt.subplots(figsize=(9.2,max(3.7,len(d)*.45)))
    ax.barh(np.arange(len(d)),d[which])
    ax.set_yticks(np.arange(len(d)),d['method'],fontsize=8)
    ax.set_xlabel(xlabel);ax.set_title(title,loc='left',fontweight='bold');ax.grid(axis='x',alpha=.18)
    for i,v in enumerate(d[which]):ax.text(v+max(d[which])*.015,i,f'{v:.3f}',va='center',fontsize=8)
    ax.set_xlim(0,max(d[which])*1.17)
    fig.tight_layout();fig.savefig(F/output,bbox_inches='tight');plt.close(fig)
graph('flock_experiment_3_summary.csv','mean_rmse','Experiment 3 | 5D, 100 replications','exp3_rmse.png','Mean RMSE (smaller is better)')
graph('flock_experiment_4_summary.csv','mean_rmse','Experiment 4 | 100 replications','exp4_rmse.png','Mean functional RMSE (smaller is better)')
graph('experiment5_complete_visible_reps_summary.csv','mean_rmse','Experiment 5 | first 25 complete visible replicates','exp5_partial_rmse.png','Mean RMSE; partial log, not full 100-replicate study')
graph('importance_pilot_reported_prior_analysis_UNVERIFIED.csv','mean_rmse','Importance-weighted pilot | 5D, 10 reported replicates','importance_pilot_reported.png','Mean RMSE (original archive unavailable for re-check)')
print('Regenerated figures from stored data in',F)