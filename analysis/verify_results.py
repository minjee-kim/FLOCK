"""Verify reported table values against the bundled historical CSV records.

This does not execute or validate the original R samplers.
"""
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'data'
for fn,col,n in [('flock_experiment_2_summary.csv','mean_rmse',100),('flock_experiment_3_summary.csv','mean_rmse',100),('flock_experiment_4_summary.csv','mean_rmse',100)]:
    df=pd.read_csv(D/fn)
    assert (df['n']==n).all(),f'Unexpected replication count in {fn}'
    assert np.isfinite(df[col]).all(),f'Nonfinite {col} in {fn}'
    print(f'{fn}: {len(df)} methods; n={n} per method')
r=pd.read_csv(D/'experiment5_visible_log_records.csv')
counts=r.groupby('rep')['method'].nunique()
complete=counts[counts==8].index
assert len(complete)==25, f'Expected 25 complete visible replicates, got {len(complete)}'
assert set(complete)==set(range(1,26)), 'Unexpected set of complete replicates'
sub=r[r['rep'].isin(complete)]
assert sub.groupby(['rep','method']).size().eq(1).all()
sum0=sub.groupby('method').agg(mean_rmse=('rmse','mean'),median_rmse=('rmse','median'),mean_weight_l1=('weight_l1','mean'))
sum1=pd.read_csv(D/'experiment5_complete_visible_reps_summary.csv').set_index('method')
for c in sum0.columns:
    assert np.allclose(sum0[c].sort_index(),sum1[c].sort_index()),c
pv=sub.pivot(index='rep',columns='method',values='rmse')
print('Exp 5 complete visible reps:',len(complete), 'Nov AM paired wins:',int((pv['FLOCK_Novelty']<pv['AdaptiveMixtureMH']).sum()))
print('Exp 5 FLOCK Novelty - AdaptiveMixture mean diff:',round(float((pv['FLOCK_Novelty']-pv['AdaptiveMixtureMH']).mean()),6))
assert (D/'importance_pilot_reported_prior_analysis_UNVERIFIED.csv').exists()
print('CAUTION: Importance-pilot values are SECONDARY historical reports, not verified from raw ZIPs.')
print('All checks of packaged summary data passed. Sampler correctness is NOT tested.')