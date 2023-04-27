# load python library
import pcDataLoader as pc
from scipy import stats
import numpy as np
import pandas as pd

# and off we go ...
ll_path = [
    ['24-12-1_a', 'output_24_12_1_off00', 'output_24_12_1_on00', 12],
    ['24-12-1_b', 'output_24_12_1_off01', 'output_24_12_1_on01', 12],
    ['24-12-1_c', 'output_24_12_1_off02', 'output_24_12_1_on02', 12],
    ['24-13-1_a', 'output_24_13_1_off00', 'output_24_13_1_on00', 13],
    ['24-13-1_b', 'output_24_13_1_off01', 'output_24_13_1_on01', 13],
    ['24-13-1_c', 'output_24_13_1_off02', 'output_24_13_1_on02', 13],
    ['24-14-1_a', 'output_24_14_1_off00', 'output_24_14_1_on00', 14],
    ['24-14-1_b', 'output_24_14_1_off01', 'output_24_14_1_on01', 14],
    ['24-14-1_c', 'output_24_14_1_off02', 'output_24_14_1_on02', 14],
    ['24-15-1_a', 'output_24_15_1_off00', 'output_24_15_1_on00', 15],
    ['24-15-1_b', 'output_24_15_1_off01', 'output_24_15_1_on01', 15],
    ['24-15-1_c', 'output_24_15_1_off02', 'output_24_15_1_on02', 15],
    ['24-16-1_a', 'output_24_16_1_off00', 'output_24_16_1_on00', 16],
    ['24-16-1_b', 'output_24_16_1_off01', 'output_24_16_1_on01', 16],
    ['24-16-1_c', 'output_24_16_1_off02', 'output_24_16_1_on02', 16],
]

ll_entropy = []
for s_process, s_random, s_ctrl, i_bin in ll_path:
    mcds_random = pc.pyMCDS(f'./{s_random}/output00000216.xml', microenv=False, graph=False)
    mcds_ctrl = pc.pyMCDS(f'./{s_ctrl}/output00000216.xml', microenv=False, graph=False)

    df_random = mcds_random.get_cell_df().loc[:,['cell_type','hamming_fract']].copy()
    df_random = pd.concat([df_random, pd.DataFrame([[0.0, n/i_bin] for n in range(1, i_bin+1)], columns=['cell_type','hamming_fract'])])

    df_ctrl = mcds_ctrl.get_cell_df().loc[:,['cell_type','hamming_fract']]
    df_ctrl = pd.concat([df_ctrl, pd.DataFrame([[0.0, n/i_bin] for n in range(1, i_bin+1)], columns=['cell_type','hamming_fract'])])

    ar_random = (df_random.groupby('hamming_fract').count() / df_random.shape[0]).loc[:,'cell_type'].values
    ar_ctrl = (df_ctrl.groupby('hamming_fract').count() / df_ctrl.shape[0]).loc[:,'cell_type'].values

    r_random = np.exp2(stats.entropy(ar_random, qk=None, base=2))
    r_ctrl = np.exp2(stats.entropy(ar_ctrl, qk=None, base=2))
    r_kldiv = np.exp2(stats.entropy(ar_ctrl, qk=ar_random, base=2)) - 1

    print(r_random, r_ctrl, r_kldiv)
    l_entropy = [s_process, i_bin, r_random, r_ctrl, r_kldiv]
    ll_entropy.append(l_entropy)

# plot box
df_entropy = pd.DataFrame(ll_entropy, columns=['run','hamming_complete','entropy_random_state','entropy_ctrl_state', 'kldiv_state'])
df_entropy.set_index('run',inplace=True)
df_entropy.plot(kind='box', by='hamming_complete', grid=True)


# plots  bar
r_random_bit = stats.entropy(ar_random, qk=None, base=2)
df_random.loc[:,'hamming_fract'].plot(
    kind='hist',
    bins=12,
    grid=True,
    title=f'random: {round(r_random_bit,3)}[bit] = {round(r_random,3)}[state] (out of 12[state])',
    color='green',
    ylim=(0,1000),
)

r_ctrl_bit = stats.entropy(ar_ctrl, qk=None, base=2)
df_ctrl.loc[:,'hamming_fract'].plot(
    kind='hist',
    bins=12,
    grid=True,
    title=f'controlled: {round(r_ctrl_bit,3)}[bit] = {round(r_ctrl,3)}[state] (out of 12[state])',
    color='maroon'
)


