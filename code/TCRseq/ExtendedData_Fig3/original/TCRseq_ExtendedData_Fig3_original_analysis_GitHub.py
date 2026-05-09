# %% [markdown]
# # TCR-seq analysis for Extended Data Fig. 3

# %% [markdown]
#
# This is the original full-analysis script for TCR gene fragment usage analysis.
#
# The original code was obtained from **Jones Salim** and adapted for the manuscript analysis.
#
# Run `TCR.R` to obtain the TCR extraction results for the two in-house samples, and integrate them with the reference TCR clone tables from Jones Salim.
#
# By default, the complete input workbook is expected at `data/TCR Clones Collated.xlsx`.
# This full workbook is not included in the demo repository. Modify `PATH_DATA_DIR` or `INPUT_WORKBOOK` below before running the full original analysis if your local paths differ.
# 

# %% [markdown]
# ## Import packages and set global properties

# %%
import pathlib
from pathlib import Path
import numpy as np
import pandas as pd
try:
    from IPython.display import display
except ImportError:
    def display(x):
        print(x)
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import PolyCollection
import seaborn as sns
from sklearn.metrics import r2_score
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.interpolate import splrep, BSpline

## settings
# User-configurable paths for the full original analysis.
# Place the complete TCR clone workbook at:
#   data/TCR Clones Collated.xlsx
# or modify PATH_DATA_DIR / INPUT_WORKBOOK according to your local environment.
PATH_DATA_DIR = pathlib.Path("data")
PATH_FIG_DIR = pathlib.Path("out")
INPUT_WORKBOOK = PATH_DATA_DIR / "TCR Clones Collated.xlsx"

# Create the output directory if it does not exist.
PATH_FIG_DIR.mkdir(parents=True, exist_ok=True)

if not INPUT_WORKBOOK.exists():
    raise FileNotFoundError(
        f"Full input workbook was not found: {INPUT_WORKBOOK}
"
        "This is the original full-analysis script. The complete input workbook "
        "is not included in the demo repository. Please place the workbook at "
        "data/TCR Clones Collated.xlsx or modify PATH_DATA_DIR / INPUT_WORKBOOK."
    )

print("Running original full TCR-seq analysis")
print(f"Input workbook: {INPUT_WORKBOOK}")
print(f"Output folder: {PATH_FIG_DIR}")

# Jupyter display setting; commented out so this .py file can run from the command line.
# %config InlineBackend.figure_formats = ['png']
plt.rcParams['figure.dpi'] = '300'     # Default is 100

# %% [markdown]
# ## Define constants/functions

# %% [markdown]
# ### Define TCR gene order
# 
# +/- pseudogenes

# %%
INCLUDE_PSEUDOGENES = False

TRBJ = ['TRBJ01-01', 'TRBJ01-02', 'TRBJ01-03', 'TRBJ01-04', 'TRBJ01-05', 'TRBJ01-06',
        'TRBJ02-01', 'TRBJ02-02', 'TRBJ02-03', 'TRBJ02-04', 'TRBJ02-05', 'TRBJ02-06', 'TRBJ02-07']

if INCLUDE_PSEUDOGENES:
    # All genes
    TRBV = ['TRBV01-01', 'TRBV02-01', 'TRBV03-01/03-02',
            'TRBV04-01', 'TRBV05-01', 'TRBV06-01', 'TRBV07-01',
            'TRBV04-02',
            'TRBV04-03',              'TRBV06-02/06-03', 'TRBV07-02', 'TRBV08-01',
                         'TRBV05-02', 'TRBV06-04', 'TRBV07-03', 'TRBV08-02',
                         'TRBV05-03',
            'TRBV09-01', 'TRBV10-01', 'TRBV11-01',
                         'TRBV10-02', 'TRBV11-02', 'TRBV12-02',
                                      'TRBV06-05', 'TRBV07-04',
                         'TRBV05-04', 'TRBV06-06', 'TRBV07-05',
                         'TRBV05-05', 'TRBV06-07', 'TRBV07-06',
                         'TRBV05-06', 'TRBV06-08', 'TRBV07-07',
                         'TRBV05-07', 'TRBV06-09', 'TRBV07-08',
                         'TRBV05-08',              'TRBV07-09', 'TRBV13-01',
                         'TRBV10-03', 'TRBV11-03', 'TRBV12-03/12-04',
                                                   'TRBV12-05',
            'TRBV14-01', 'TRBV15-01', 'TRBV16-01', 'TRBV17-01', 'TRBV18-01',
            'TRBV19-01', 'TRBV20-01', 'TRBV20-or9-02', 'TRBV21-01',
            'TRBV22-01', 'TRBV23-01', 'TRBV24-01', 'TRBV25-01',
            'TRBVA-01',  'TRBV26-01', 'TRBVB-01',  'TRBV27-01',
            'TRBV28-01', 'TRBV29-01', 'TRBV30-01']
else:
    # Sans pseudogenes
    TRBV = ['TRBV02-01', 'TRBV03-01/03-02',
                         'TRBV04-01', 'TRBV05-01', 'TRBV06-01', 'TRBV07-01',
            'TRBV04-02',
                         'TRBV04-03',              'TRBV06-02/06-03', 'TRBV07-02',
                                      'TRBV06-04', 'TRBV07-03',
                         'TRBV05-03',
            'TRBV09-01', 'TRBV10-01', 'TRBV11-01',
                         'TRBV10-02', 'TRBV11-02',
                                      'TRBV06-05', 'TRBV07-04',
                         'TRBV05-04', 'TRBV06-06',
                         'TRBV05-05', 'TRBV06-07', 'TRBV07-06',
                         'TRBV05-06', 'TRBV06-08', 'TRBV07-07',
                         'TRBV05-07', 'TRBV06-09', 'TRBV07-08',
                         'TRBV05-08',              'TRBV07-09', 'TRBV13-01',
                         'TRBV10-03', 'TRBV11-03', 'TRBV12-03/12-04',
                                                   'TRBV12-05',
                         'TRBV14-01', 'TRBV15-01', 'TRBV16-01', 'TRBV17-01', 'TRBV18-01',
            'TRBV19-01', 'TRBV20-01', 'TRBV21-01',
                         'TRBV23-01', 'TRBV24-01', 'TRBV25-01',
                         'TRBV27-01', 'TRBV28-01', 'TRBV29-01', 'TRBV30-01']

print(len(TRBV))


if INCLUDE_PSEUDOGENES:
    # All genes
    TRAV = ['TRAV01-01',     'TRAV01-02', 'TRAV02-01', 'TRAV03-01', 'TRAV04-01', 'TRAV05-01',
            'TRAV06-01',     'TRAV07-01', 'TRAV08-01', 'TRAV09-01', 'TRAV10-01', 'TRAV11-01',
            'TRAV12-01',     'TRAV08-02', 'TRAV08-03', 'TRAV13-01',
            'TRAV12-02',     'TRAV08-04', 'TRAV08-05', 'TRAV13-02',
            'TRAV14-01/DV4', 'TRAV09-02', 'TRAV15-01',
            'TRAV12-03',     'TRAV08-06', 'TRAV08-07',
            'TRAV16-01',     'TRAV17-01', 'TRAV18-01', 'TRAV19-01', 'TRAV20-01', 'TRAV21-01', 'TRAV22-01',
            'TRAV23-01/DV6', 'TRDV1-01',  'TRAV24-01', 'TRAV25-01',
            'TRAV26-01',     'TRAV08-07',
            'TRAV27-01',     'TRAV28-01', 'TRAV29-01', 'TRAV30-01', 'TRAV31-01', 'TRAV32-01', 'TRAV33-01',
            'TRAV26-02',     'TRAV34-01', 'TRAV35-01',
            'TRAV36-01/DV7', 'TRAV37-01', 'TRAV38-01',
            'TRAV38-02/DV8', 'TRAV39-01', 'TRAV40-01', 'TRAV41-01']

    TRAJ = ['TRAJ' + str(x) for x in range(61, 1, -1) + '-01']

else:
    # Sans pseudogenes
    TRAV = ['TRAV01-01',     'TRAV01-02', 'TRAV02-01', 'TRAV03-01', 'TRAV04-01', 'TRAV05-01',
            'TRAV06-01',     'TRAV07-01', 'TRAV08-01', 'TRAV09-01', 'TRAV10-01',
            'TRAV12-01',     'TRAV08-02', 'TRAV08-03', 'TRAV13-01',
            'TRAV12-02',     'TRAV08-04',              'TRAV13-02',
            'TRAV14-01/DV4', 'TRAV09-02',
            'TRAV12-03',     'TRAV08-06',
            'TRAV16-01',     'TRAV17-01', 'TRAV19-01', 'TRAV20-01', 'TRAV21-01', 'TRAV22-01',
            'TRAV23-01/DV6', 'TRAV24-01', 'TRAV25-01',
            'TRAV26-01',
            'TRAV27-01',     'TRAV29-01', 'TRAV30-01',
            'TRAV26-02',     'TRAV34-01', 'TRAV35-01',
            'TRAV36-01/DV7', 'TRAV38-01',
            'TRAV38-02/DV8', 'TRAV39-01', 'TRAV40-01', 'TRAV41-01']

    TRAJ = ['TRAJ58-01', 'TRAJ57-01', 'TRAJ56-01',
                         'TRAJ54-01', 'TRAJ53-01', 'TRAJ52-01',
            'TRAJ50-01', 'TRAJ49-01', 'TRAJ48-01', 'TRAJ47-01', 'TRAJ46-01',
            'TRAJ45-01', 'TRAJ44-01', 'TRAJ43-01', 'TRAJ42-01', 'TRAJ41-01',
            'TRAJ40-01', 'TRAJ39-01', 'TRAJ38-01', 'TRAJ37-01', 'TRAJ36-01',
                         'TRAJ34-01', 'TRAJ33-01', 'TRAJ32-01', 'TRAJ31-01',
            'TRAJ30-01', 'TRAJ29-01', 'TRAJ28-01', 'TRAJ27-01', 'TRAJ26-01',
            'TRAJ25-01', 'TRAJ24-01', 'TRAJ23-01', 'TRAJ22-01', 'TRAJ21-01',
            'TRAJ20-01',              'TRAJ18-01', 'TRAJ17-01',
            'TRAJ15-01',              'TRAJ13-01', 'TRAJ12-01', 'TRAJ11-01',
            'TRAJ10-01', 'TRAJ09-01', 'TRAJ08-01', 'TRAJ07-01', 'TRAJ06-01',
            'TRAJ05-01', 'TRAJ04-01', 'TRAJ03-01']

# %% [markdown]
# ## Process data

# %% [markdown]
# ### Load data

# %%
datasets = ['PBMC Collated 2 filtered', 'LD-SC40 Revision',
            'Fetal_w13_FCAImm', 'Fetal_w14_FCAImm', 'Fetal_thy_w17',
            'Postnatal_thy_10m', 'Postnatal_30m_T06', 'Adult_thy', 'DMSO-T', 'Inhibitor-T']

df_all = pd.read_excel(
    INPUT_WORKBOOK, sheet_name=datasets, header=0)

# %% [markdown]
# ### Label subsets

# %%
# Label different subsets
# - cell_type: standard cell types
# - cell_type_ext: extended mem/eff for iPSCs
for d in df_all:

    for index, row in df_all[d].iterrows():

        # NOTE: Don't do these as match-case since we're looking for substrings

        if ('gdT' in row['Sample']) or ('ILC' in row['Sample']) or ('CD8aa' in row['Sample']) or ('CD8αα' in row['Sample']):
            cell_type = 'Innate'
            cell_type_ext = 'Innate'
        elif ('Met-Act/Stress' in row['Sample']):
            cell_type = 'Stressed'
            cell_type_ext = 'Stressed'

        # elif 'DP (P)' in row['Sample']: # From old code
        #     cell_type = 'DP (P)'
        elif 'DP' in row['Sample']:
            cell_type = 'DP'
            cell_type_ext = 'DP'

        # From old code
        # elif 'CD4 immature' in row['Sample']:
        #     cell_type = 'CD4_immature'
        # elif 'CD4 mature' in row['Sample']:
        #     cell_type = 'CD4_mature'

        # This is if we want to break out different mem/eff types
        elif ('4SP_N/CM' in row['Sample']) or ('4SP_TM' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4_N/CM/TM'
        elif ('8SP_N/CM' in row['Sample']) or ('8SP_TM' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8_N/CM/TM'
        elif ('4SP_EM' in row['Sample']) or ('4SP_Mem(P)' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4_EM/Mem(P)'
        elif ('8SP_EM' in row['Sample']) or ('8SP_Mem(P)' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8_EM/Mem(P)'
        elif ('4SP_Eff/Ex' in row['Sample']) or ('4SP_Eff(P)' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4_Eff/Ex/Eff(P)'
        elif ('8SP_Eff/Ex' in row['Sample']) or ('8SP_Eff(P)' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8_Eff/Ex/Eff(P)'

        elif ('4SP' in row['Sample']) or ('CD4' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4'
        elif ('8SP' in row['Sample']) or ('CD8' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8'
        else:
            cell_type = None
            cell_type_ext = None

        df_all[d].loc[index, 'Cell Type'] = cell_type
        df_all[d].loc[index, 'Cell Type Ext'] = cell_type_ext

    if 'Fetal' in d:
        df_all[d]['Subset'] = 'Fetal Thym'
    elif 'Postnatal' in d:
        df_all[d]['Subset'] = 'Postnatal Thym'
    elif 'Adult' in d:
        df_all[d]['Subset'] = 'Adult Thym'
    elif 'PBMC' in d:
        df_all[d]['Subset'] = 'PBMC'
    elif 'DMSO' in d:
        df_all[d]['Subset'] = 'DMSO'
    elif 'Inhibitor' in d:
        df_all[d]['Subset'] = 'Inhibitor'
    else:
        df_all[d]['Subset'] = 'iPSC'

    df_all[d]['Subset Cell Type'] = df_all[d]['Subset'] + \
        ' ' + df_all[d]['Cell Type']
    df_all[d]['Subset Cell Type Ext'] = df_all[d]['Subset'] + \
        ' ' + df_all[d]['Cell Type Ext']

# %% [markdown]
# ### Label subsets
# 
# remove NA VB

# %%
for d in df_all:
    current_df = df_all[d]

    filtered_df = current_df.dropna(subset=['CDR3_B'], inplace=False)
    filtered_df = filtered_df.dropna(subset=['CDR3_A'], inplace=False)
    nan_count = filtered_df['CDR3_B'].isna().sum()
    print(f"dataset {d} filtered CDR3_B column number of NaN：{nan_count}\n")

    df_all[d] = filtered_df

    for index, row in df_all[d].iterrows():
        # NOTE: Don't do these as match-case since we're looking for substrings
        if ('gdT' in row['Sample']) or ('ILC' in row['Sample']) or ('CD8aa' in row['Sample']) or ('CD8αα' in row['Sample']):
            cell_type = 'Innate'
            cell_type_ext = 'Innate'
        elif ('Met-Act/Stress' in row['Sample']):
            cell_type = 'Stressed'
            cell_type_ext = 'Stressed'
        elif 'DP' in row['Sample']:
            cell_type = 'DP'
            cell_type_ext = 'DP'
        elif ('4SP_N/CM' in row['Sample']) or ('4SP_TM' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4_N/CM/TM'
        elif ('8SP_N/CM' in row['Sample']) or ('8SP_TM' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8_N/CM/TM'
        elif ('4SP_EM' in row['Sample']) or ('4SP_Mem(P)' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4_EM/Mem(P)'
        elif ('8SP_EM' in row['Sample']) or ('8SP_Mem(P)' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8_EM/Mem(P)'
        elif ('4SP_Eff/Ex' in row['Sample']) or ('4SP_Eff(P)' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4_Eff/Ex/Eff(P)'
        elif ('8SP_Eff/Ex' in row['Sample']) or ('8SP_Eff(P)' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8_Eff/Ex/Eff(P)'
        elif ('4SP' in row['Sample']) or ('CD4' in row['Sample']):
            cell_type = 'CD4'
            cell_type_ext = 'CD4'
        elif ('8SP' in row['Sample']) or ('CD8' in row['Sample']):
            cell_type = 'CD8'
            cell_type_ext = 'CD8'
        else:
            cell_type = None
            cell_type_ext = None

        df_all[d].loc[index, 'Cell Type'] = cell_type
        df_all[d].loc[index, 'Cell Type Ext'] = cell_type_ext

    if 'Fetal' in d:
        df_all[d]['Subset'] = 'Fetal Thym'
    elif 'Postnatal' in d:
        df_all[d]['Subset'] = 'Postnatal Thym'
    elif 'Adult' in d:
        df_all[d]['Subset'] = 'Adult Thym'
    elif 'PBMC' in d:
        df_all[d]['Subset'] = 'PBMC'
    elif 'DMSO' in d:
        df_all[d]['Subset'] = 'DMSO'
    elif 'Inhibitor' in d:
        df_all[d]['Subset'] = 'Inhibitor'
    else:
        df_all[d]['Subset'] = 'iPSC'

    df_all[d]['Subset Cell Type'] = df_all[d]['Subset'] + \
        ' ' + df_all[d]['Cell Type'].fillna('Unknown')
    df_all[d]['Subset Cell Type Ext'] = df_all[d]['Subset'] + \
        ' ' + df_all[d]['Cell Type Ext'].fillna('Unknown')

# %% [markdown]
# ### Clean up the names of TCR fragments
# 
# * Remove allele designations from TCR gene names
# * Add '-01' to the end of any fragment missing such an indicator
#     * Some sources do/don't include if only one e.g. TRAV05 vs TRAV05-01
# * Set all names to TRAJ, TRAV, TRBJ, TRBV, rather than TCRA/B/V/J 
# 
# Compute CDR3 lengths if not already in the table

# %% [markdown]
# #### Functions

# %%
# A few possible exceptions remaining include:
# - Tandem names: TRBV03-01/03-02, TRBV06-02/06-03, and TRBV12-03/12-04

def fixdash(x):
    if x in [None, np.nan]:
        return x

    if len(x) < 6 or x[5] == '-':
        # Fix missing '0' for single-digit fragment ID eg 1, 2, 3...
        x = x[:4] + '0' + x[4:]

    # Handle missing second part of fragment ID (often omitted)
    if len(x) < 7:
        x = x + '-01'
    else:
        # Handle presence of DV chain IDs
        x = x if x[6] == '-' else x[:6] + '-01' + x[6:]

    # Ensure after-dash part is two digits
    if len(x) < 9 or x[8] not in '0123456789':
        x = x[:7] + '0' + x[7:]

    return x


def fixdv(x):
    if x in [None, np.nan]:
        return x
    elif 'DV' in x:
        return x if '/' in x else x[:9] + '/' + x[9:]
    else:
        if 'TRAV14' in x:
            return x + '/DV4'
        elif 'TRAV23' in x:
            return x + '/DV6'
        elif 'TRAV29' in x:
            return x + '/DV5'
        elif 'TRAV36' in x:
            return x + '/DV7'
        elif 'TRAV38-02' in x:
            return x + '/DV8'
        else:
            return x


test = pd.DataFrame({'Chain': [None, np.nan, 'TRBJ2', 'TRAV01', 'TRBV01-01', 'TRBJ01-02',
                    'TRAJ06/DV6', 'TRAV08-01/DV4', 'TRAV29DV5', 'TRAV36/DV7', 'TRAV38-02', 'TRBV6-3']})
test['Chain_fixed_dash'] = test['Chain'].apply(fixdash)
test['Chain_fixed_dv'] = test['Chain_fixed_dash'].apply(fixdv)
display(test)

# %% [markdown]
# #### Apply functions

# %%
df_coll = pd.concat([sheet for _, sheet in df_all.items()]).reset_index()
# 'Cell Type' and 'Cell Type Ext' both always have paired data, so only need to check one
df = df_coll.loc[df_coll['Cell Type'].notna(), :].copy()

for chain in ['VB', 'JB', 'VA', 'JA']:
    # Strip allele information (most efficient to do these here)
    df[chain] = df[chain].str.replace(r'\*.*', '', regex=True)
    df[chain] = df[chain].str.replace(
        'TCR', 'TR', regex=True)  # Normalize naming convention
    df[chain] = df[chain].apply(fixdash)
    df[chain] = df[chain].apply(fixdv)

# Ensure CDR3 lengths are calculated
df['CDR3_B AA'] = df['CDR3_B'].apply(
    lambda x: len(x) if x not in [None, np.nan] else x)
df['CDR3_B bp'] = df['CDR3_B AA'].apply(lambda x: 3*x)
df['CDR3_A AA'] = df['CDR3_A'].apply(
    lambda x: len(x) if x not in [None, np.nan] else x)
df['CDR3_A bp'] = df['CDR3_A AA'].apply(lambda x: 3*x)

# display(df_coll)
display(df)

# %% [markdown]
# #### Drop unnamed columns

# %%
for col in df:
    if 'Unnamed' in col:
        df.drop([col], axis='columns', inplace=True)

display(df)

# %% [markdown]
# ### Compute counts for each V-J combination per subset

# %%
# Sum TRB V/J combos, then V and J gene totals separately and CDR3 lengths
# NOTE: adjusted for Template #s (reads)
# NOTE: using simple version for standard calculations here
df_counts_ab = df.pivot_table(values='Templates', index='Subset Cell Type', columns=[
                              'CDR3_A', 'VA', 'JA', 'CDR3_B', 'VB', 'JB'], aggfunc='sum')  # .reset_index()

df_counts_vjca = df.pivot_table(values='Templates', index='Subset Cell Type', columns=[
                                'CDR3_A', 'VA', 'JA'], aggfunc='sum')  # .reset_index()
df_counts_vja = df.pivot_table(values='Templates', index='Subset Cell Type', columns=[
                               'VA', 'JA'], aggfunc='sum', fill_value=0)  # .reset_index()
df_counts_va = df.pivot_table(values='Templates', index='Subset Cell Type',
                              columns='VA',        aggfunc='sum', fill_value=0)  # .reset_index()
df_counts_ja = df.pivot_table(values='Templates', index='Subset Cell Type',
                              columns='JA',        aggfunc='sum', fill_value=0)  # .reset_index()

df_counts_vjcb = df.pivot_table(values='Templates', index='Subset Cell Type', columns=[
                                'CDR3_B', 'VB', 'JB'], aggfunc='sum')  # .reset_index()
df_counts_vjb = df.pivot_table(values='Templates', index='Subset Cell Type', columns=[
                               'VB', 'JB'], aggfunc='sum', fill_value=0)  # .reset_index()
df_counts_vb = df.pivot_table(values='Templates', index='Subset Cell Type',
                              columns='VB',        aggfunc='sum', fill_value=0)  # .reset_index()
df_counts_jb = df.pivot_table(values='Templates', index='Subset Cell Type',
                              columns='JB',        aggfunc='sum', fill_value=0)  # .reset_index()

df_cdr3a = df[df['CDR3_A AA'].notna()].pivot_table(values='Templates', index='Subset Cell Type',
                                                   columns='CDR3_A AA', aggfunc='sum', sort=True)  # .reset_index()
df_cdr3b = df[df['CDR3_B AA'].notna()].pivot_table(values='Templates', index='Subset Cell Type',
                                                   columns='CDR3_B AA', aggfunc='sum', sort=True)  # .reset_index()

display(df_counts_ab)

display(df_counts_vjca)
display(df_counts_vja)
display(df_counts_va)
display(df_counts_ja)
display(df_cdr3a)

display(df_counts_vjcb)
display(df_counts_vjb)
display(df_counts_vb)
display(df_counts_jb)
display(df_cdr3b)

# %% [markdown]
# ### Set subset colors, linestyles, and markers

# %%
# Define colors
cmaps = {
    'blue':   plt.get_cmap('Blues'),
    'green':  plt.get_cmap('Greens'),
    'orange': plt.get_cmap('Oranges'),
    'red':    plt.get_cmap('Reds'),
    'purple': plt.get_cmap('Purples'),
    'gray':   plt.get_cmap('Greys')
}
colors = {
    'iPSC CD4_mature': cmaps['blue'](1.0), 'iPSC CD4_immature': cmaps['blue'](0.8), 'iPSC DP': cmaps['blue'](0.6), 'iPSC DP (P)': cmaps['blue'](0.4),
    'iPSC CD8': cmaps['blue'](0.2), 'iPSC Innate': cmaps['green'](0.5),
    'Fetal Thym CD4': cmaps['orange'](1.0), 'Fetal Thym DP': cmaps['orange'](0.66), 'Fetal Thym CD8': cmaps['orange'](0.33),
    'Postnatal Thym CD4': cmaps['orange'](1.0), 'Postnatal Thym DP': cmaps['orange'](0.66), 'Postnatal Thym CD8': cmaps['orange'](0.33),
    'Adult Thym CD4': cmaps['orange'](1.0), 'Adult Thym DP': cmaps['orange'](0.66), 'Adult Thym CD8': cmaps['orange'](0.33),
    'PBMC CD4': cmaps['red'](0.5), 'PBMC CD8': cmaps['red'](0.4),
    'DMSO CD8': cmaps['purple'](0.5), 'Inhibitor CD8': cmaps['gray'](0.4)}

plots = {
    'Simple': {'FigH': 2.5, 'Subset': 'Subset Cell Type'},
    # 'Ext':    {'FigH': 4.5, 'Subset': 'Subset Cell Type Ext'}
}

orders = {
    'Simple': [

        # 'iPSC DP', 'iPSC CD4', 'iPSC CD8',
        # 'Fetal Thym DP', 'Fetal Thym CD4',
        'Fetal Thym CD8',
        # 'Postnatal Thym DP', 'Postnatal Thym CD4',
        'Postnatal Thym CD8',
        # 'PBMC CD4',
        'PBMC CD8',
        'DMSO CD8',
        'Inhibitor CD8',
    ],
    'Ext': ['iPSC DP', 'iPSC CD4_N/CM/TM', 'iPSC CD4_EM/Mem(P)', 'iPSC CD4_Eff/Ex/Eff(P)',
            'iPSC CD8_N/CM/TM', 'iPSC CD8_EM/Mem(P)', 'iPSC CD8_Eff/Ex/Eff(P)',
            'Fetal Thym DP', 'Fetal Thym CD4', 'Fetal Thym CD8',
            'Postnatal Thym DP', 'Postnatal Thym CD4', 'Postnatal Thym CD8',
            'PBMC CD4', 'PBMC CD8']
}


def get_color(subname):

    # NOTE: Don't do these as match-case since we're looking for substrings
    if 'iPSC' in subname:
        cweight = 0.1
    elif 'Fetal' in subname:
        cweight = 0.2
    elif 'Postnatal' in subname:
        cweight = 0.4
    elif 'PBMC' in subname:
        cweight = 0.6
    elif 'DMSO' in subname:
        cweight = 0.8
    elif 'Inhibitor' in subname:
        cweight = 1.0

    if 'N/CM/TM' in subname:
        cweight -= 0.075  # slightly lighter
    elif 'EM/Mem(P)' in subname:
        cweight += 0     # stay same
    elif 'Eff/Ex/Eff(P)' in subname:
        cweight += 0.075  # slightly darker

    if 'CD4' in subname:
        color = plt.get_cmap('Blues')(cweight)
    elif 'DP' in subname:
        color = plt.get_cmap('Greens')(cweight)
    elif 'CD8' in subname:
        color = plt.get_cmap('Oranges')(cweight)

    return color


def get_linestyle(subname):

    if 'iPSC' in subname:
        linestyle = 'solid'
    elif 'Fetal' in subname:
        linestyle = 'dashed'
    elif 'Postnatal' in subname:
        linestyle = 'dashdot'
    elif 'PBMC' in subname:
        linestyle = 'dotted'
    elif 'DMSO' in subname:
        linestyle = (0, (5, 5))
    elif 'Inhibitor' in subname:
        linestyle = (0, (3, 1, 1, 1))

    return linestyle


def get_marker(subname):

    if 'iPSC' in subname:
        marker = 'o'
    elif 'Fetal' in subname:
        marker = '^'
    elif 'Postnatal' in subname:
        marker = 'v'
    elif 'PBMC' in subname:
        marker = 's'
    elif 'DMSO' in subname:
        marker = 'd'
    elif 'Inhibitor' in subname:
        marker = '*'

    return marker

# %% [markdown]
# ### Compare diversity metrics

# %% [markdown]
# #### Diversity functions
# 
# **Shannon Entropy** ($f$ = frequency of each clone): 
# 
# $H(f) = -\sum f_i \cdot log_2(f_i)$
# 
# Normalized for length ($T$ = total # TCRs clones):
# 
# $H(f) = \frac{-1}{log_2(T)} \sum f_i \cdot log_2(f_i)$
# 
# **Simpson's Diversity** ($n$ = number of templates for each clone -- https://www.statology.org/simpsons-diversity-index)
# * Interpreted as the unbiased probability of sampling the same thing twice
# * N = total number of templates (sequencing reads)
# 
# $D = \frac{1}{N \cdot (N - 1)} \cdot \sum_{i=0}^T n_i \cdot (n_i - 1)$
# 
# This can be inverted to **Simpson's Index** ($1-D$) to make it so larger numbers mean higher diversity, which is more intuitive and aligns with the format of Shannon entropy.
# 
# **Simpson clonality**, as specified by [Adaptive Bio](https://www.adaptivebiotech.com/wp-content/uploads/2020/06/immunoSEQ_Analyzer-Tech-Note_Clonality_WEB_MRK-00355.pdf), is the square root of Simpson's Diversity:
# 
# $C = \sqrt{D}$
# 
# The variance can be estimated in an unbiased manner by taking a higher moment (see Andreas Tiffeau-Meyer's [paper](doi.org/10.1103/PhysRevE.109.064411))
# 
# $var(D) = \frac{a}{1-b} \cdot D_T - \frac{b}{1-b} \cdot D^2 + \frac{c}{1-b} \cdot D$
# 
# $D_T = \frac{1}{N \cdot (N - 1) \cdot (N - 2)} \cdot \sum_{i=0}^T n_i \cdot (n_i - 1) \cdot (n_i - 2)$
# 
# $a = \frac{4(N-2)}{N(N-1)}, b = \frac{2(2N-3)}{N(N-1)}, c = \frac{2}{N(N-1)}$

# %%
def shannon(n, normalize=False):
    # n should be a vector of counts (optional: probabilities with sum = 1)
    # Optional normalization input normalizes for the number of elements sampled

    # convert n to  porbability vector
    if type(n) == list:
        n = np.array(n)
    f = n / sum(n)
    entropy = -sum(f * np.log2(f))

    # Normalize for length of n
    if normalize:
        entropy /= np.log2(len(f))

    return entropy


def simpson(n, metric='diversity'):
    # Unbiased estimator of sample diversity. Metrics:
    # - 'diversity' (Default): Simpson's Diversity Index -- smaller # = more diverse
    # - 'index':               Simpson's Index of Diversity (1 - Diversity) -- larger # = more diverse
    # - 'clonality':           Simpson Clonality (sqrt of Diversity) -- as computed by Adaptive Bio

    if type(n) == list:
        n = np.array(n)
    N = sum(n)
    D = sum(n * (n - 1)) / (N * (N - 1))

    metric = metric.lower()
    if metric in ['div', 'diversity', 'diversity index']:
        return D
    elif metric in ['ind', 'index', 'index of diversity']:
        return 1 - D
    elif metric in ['clon', 'clonality']:
        return np.sqrt(D)
    else:
        raise Exception(f'Invalid metric: {metric}')


def simpson_var(n, mode='variance'):
    # Higher moment of Simpson diversity index. Modes:
    # - 'variance' (Default): Variance estimate
    # - 'std':                Standard deviation estimate (sqrt(Variance))

    if type(n) == list:
        n = np.array(n)
    N = sum(n)
    D = sum(n * (n - 1)) / (N * (N - 1))
    D_T = sum(n * (n - 1) * (n - 2)) / (N * (N - 1) * (N - 2))
    a = 4 * (N - 2) / (N * (N - 1))
    b = 2 * (2*N - 3) / (N * (N - 1))
    c = 2 / (N * (N - 1))
    D_var = a / (1 - b) * D_T - b / (1 - b) * D**2 + c / (1 - b) * D

    mode = mode.lower()
    if mode in ['var', 'variance']:
        return D_var
    elif mode in ['std', 'stdev', 'standard deviation']:
        return np.sqrt(D_var)
    else:
        raise Exception(f'Invalid mode: {metric}')


# shannon(subsets['CVJ']['Thym']['Templates'], normalize = True)

# %% [markdown]
# #### Plots

# %%
from scipy.stats import mannwhitneyu

subplot_titles = ['TCRα chain', 'TCRβ chain']

for plotname, plotinfo in plots.items():

    order = np.flip(orders[plotname])

    figd, axsd = plt.subplots(1, 2, figsize=(4, plotinfo['FigH']), sharey=True)
    figd.tight_layout(w_pad=1)

    df_counts_vjca_curr = df.pivot_table(
        values='Templates',
        index=plotinfo['Subset'],
        columns=['CDR3_A', 'VA', 'JA'],
        aggfunc='sum'
    )

    df_counts_vjcb_curr = df.pivot_table(
        values='Templates',
        index=plotinfo['Subset'],
        columns=['CDR3_B', 'VB', 'JB'],
        aggfunc='sum'
    )

    for di, df_vjc in enumerate([df_counts_vjca_curr, df_counts_vjcb_curr]):

        axd = axsd[di]
        bar_counters = [0, 0]

        dmso_data = None
        inh_data = None

        for subname in order:

            color = get_color(subname)

            data = df_vjc.loc[subname, :]
            data = data[data.notna()]

            if subname == "DMSO CD8":
                dmso_data = data

            if subname == "Inhibitor CD8":
                inh_data = data

            entropy = shannon(data, normalize=True)
            sindex = simpson(data, metric='index')
            sdiversity = 1 - sindex

            print(f'{subname}: len: {len(data)}, sum: {sum(data)}, entropy: {entropy}, sindex: {sindex}, sdiversity: {sdiversity}')

            sindex_std = simpson_var(data, mode='std')

            axd.barh(bar_counters[di], sdiversity, color=color)

            axd.barh(
                bar_counters[di],
                sdiversity,
                xerr=sindex_std,
                capsize=3,
                error_kw={'elinewidth': 0.5, 'capthick': 0.5},
                color=color
            )

            bar_counters[di] += 1

        # ---------------------------
        # cal Wilcoxon rank-sum test
        # ---------------------------

        if dmso_data is not None and inh_data is not None:

            stat, p = mannwhitneyu(dmso_data, inh_data,
                                   alternative='two-sided')

            chain = subplot_titles[di]

            print(f"{chain} Wilcoxon rank-sum test (DMSO vs Inhibitor): p = {p}")

            y_dmso = list(order).index("DMSO CD8")
            y_inh = list(order).index("Inhibitor CD8")

            y_low = min(y_dmso, y_inh)
            y_high = max(y_dmso, y_inh)

            x_bracket = 1e-4

            axd.plot([x_bracket, x_bracket], [
                     y_low, y_high], color='black', lw=1)

            axd.plot([x_bracket*0.7, x_bracket],
                     [y_low, y_low], color='black', lw=1)
            axd.plot([x_bracket*0.7, x_bracket],
                     [y_high, y_high], color='black', lw=1)

            axd.text(
                x_bracket*1.3,
                (y_low + y_high) / 2,
                f"p={p:.2e}",
                va='center',
                fontsize=9
            )

        # ---------------------------
        # pos
        # ---------------------------

        axd.set(
            xscale='log',
            xticks=[1e-6, 1e-4, 1e-2, 1e0],
            xlim=[1e-6, 1],
            ylim=[-0.6, len(order) - 0.4]
        )

        axd.set_yticks(np.arange(len(order)))
        axd.set_yticklabels(order)

        axd.set_xlabel("Simpson's Diversity")

        axd.set_title(
            subplot_titles[di],
            fontsize=10,
            pad=8,
            ha='center'
        )

    figd.savefig(
        PATH_FIG_DIR /
        f'Fig.A - Simpson diversity.pdf',
        bbox_inches='tight',
        facecolor='#ffffff00'
    )

# %% [markdown]
# ### Compare TCR fragment usages
# 
# Compare fractional usages: PSC vs CB vs thymocytes vs PBMCs
# 
# Z-score across rows and cols to compare relative usages

# %% [markdown]
# #### PCA

# %%
cols = [TRBV, TRBJ, TRAV, TRAJ]
pca = PCA(n_components=2)

for plotname, plotinfo in plots.items():

    order = orders[plotname]

    df_counts_va_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VA', aggfunc='sum', fill_value=0)
    df_counts_ja_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JA', aggfunc='sum', fill_value=0)
    df_counts_vb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VB', aggfunc='sum', fill_value=0)
    df_counts_jb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JB', aggfunc='sum', fill_value=0)

    # Fit PCA on concatenated data for each cell type (rows = cell types, cols = all TCR segments)
    data = []
    for di, df_x in enumerate([df_counts_vb_curr, df_counts_jb_curr, df_counts_va_curr, df_counts_ja_curr]):

        # Normalize segment frequencies
        data.append(df_x.copy().reindex(order, columns=cols[di]).fillna(
            0).apply(lambda x: x / x.sum(), axis=1))

    seg_pca = pca.fit_transform(pd.concat(data, axis=1))

    print(f"Explained variance ratios for {plotname}:")
    for i, ratio in enumerate(pca.explained_variance_ratio_):
        print(f"  PC{i+1}: {ratio:.4f}")
    print(f"  Cumulative: {pca.explained_variance_ratio_.cumsum()[-1]:.4f}\n")

    df_seg_pca = pd.DataFrame(seg_pca, index=order)

    fig, axs = plt.subplots(1, 1, figsize=(3, 3))

    for subname in order:

        axs.scatter(df_seg_pca.loc[subname, 0], df_seg_pca.loc[subname, 1], marker=get_marker(
            subname), s=40, color=get_color(subname))

    # axs.set(xlabel = 'Principal Component 1', ylabel = 'Principal Component 2')

    plt.legend(order, loc='center left', bbox_to_anchor=(1.05, 0.5))

# %% [markdown]
# ### Modify PCA dim

# %%
cols = [TRBV, TRBJ, TRAV, TRAJ]
pca = PCA(n_components=3)

for plotname, plotinfo in plots.items():

    order = orders[plotname]

    df_counts_va_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VA', aggfunc='sum', fill_value=0)
    df_counts_ja_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JA', aggfunc='sum', fill_value=0)
    df_counts_vb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VB', aggfunc='sum', fill_value=0)
    df_counts_jb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JB', aggfunc='sum', fill_value=0)

    # Fit PCA on concatenated data for each cell type (rows = cell types, cols = all TCR segments)
    data = []
    for di, df_x in enumerate([df_counts_vb_curr, df_counts_jb_curr, df_counts_va_curr, df_counts_ja_curr]):

        # Normalize segment frequencies
        data.append(df_x.copy().reindex(order, columns=cols[di]).fillna(
            0).apply(lambda x: x / x.sum(), axis=1))

    seg_pca = pca.fit_transform(pd.concat(data, axis=1))
    df_seg_pca = pd.DataFrame(seg_pca, index=order)

    fig, axs = plt.subplots(1, 1, figsize=(3, 3))

    for subname in order:

        axs.scatter(df_seg_pca.loc[subname, 0], df_seg_pca.loc[subname, 2],
                    marker=get_marker(subname), s=40, color=get_color(subname))

    # axs.set(xlabel = 'Principal Component 1', ylabel = 'Principal Component 3')

    plt.legend(order, loc='center left', bbox_to_anchor=(1.05, 0.5))

    axs.set_xlabel('PC1')
    axs.set_ylabel('PC3')
    # Sets area outside ax to transparent
    # fig.savefig(
    #     PATH_FIG_DIR / f'TCR Segment PCA1-3 {plotname}.pdf', bbox_inches='tight', facecolor='#ffffff00')
    # plt.show()

cols = [TRBV, TRBJ, TRAV, TRAJ]
pca = PCA(n_components=3)

for plotname, plotinfo in plots.items():

    order = orders[plotname]

    df_counts_va_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VA', aggfunc='sum', fill_value=0)
    df_counts_ja_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JA', aggfunc='sum', fill_value=0)
    df_counts_vb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VB', aggfunc='sum', fill_value=0)
    df_counts_jb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JB', aggfunc='sum', fill_value=0)

    # Fit PCA on concatenated data for each cell type (rows = cell types, cols = all TCR segments)
    data = []
    for di, df_x in enumerate([df_counts_vb_curr, df_counts_jb_curr, df_counts_va_curr, df_counts_ja_curr]):

        # Normalize segment frequencies
        data.append(df_x.copy().reindex(order, columns=cols[di]).fillna(
            0).apply(lambda x: x / x.sum(), axis=1))

    seg_pca = pca.fit_transform(pd.concat(data, axis=1))

    # ratio
    print(f"Explained variance ratios for {plotname}:")
    for i, ratio in enumerate(pca.explained_variance_ratio_):
        print(f"  PC{i+1}: {ratio:.4f}")
    print(f"  Cumulative: {pca.explained_variance_ratio_.cumsum()[-1]:.4f}\n")

    df_seg_pca = pd.DataFrame(seg_pca, index=order)

    fig, axs = plt.subplots(1, 1, figsize=(3, 3))

    for subname in order:

        axs.scatter(df_seg_pca.loc[subname, 0], df_seg_pca.loc[subname, 1],
                    marker=get_marker(subname), s=40, color=get_color(subname))

    # axs.set(xlabel = 'Principal Component 1', ylabel = 'Principal Component 3')

    plt.legend(order, loc='center left', bbox_to_anchor=(1.05, 0.5))

    axs.set_xlabel('PC1')
    axs.set_ylabel('PC2')
    # Sets area outside ax to transparent
    # fig.savefig(
    #     PATH_FIG_DIR / f'TCR Segment PCA1-2 {plotname}.pdf', bbox_inches='tight', facecolor='#ffffff00')
    # plt.show()

cols = [TRBV, TRBJ, TRAV, TRAJ]
pca = PCA(n_components=3)

for plotname, plotinfo in plots.items():

    order = orders[plotname]

    df_counts_va_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VA', aggfunc='sum', fill_value=0)
    df_counts_ja_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JA', aggfunc='sum', fill_value=0)
    df_counts_vb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VB', aggfunc='sum', fill_value=0)
    df_counts_jb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JB', aggfunc='sum', fill_value=0)

    # Fit PCA on concatenated data for each cell type (rows = cell types, cols = all TCR segments)
    data = []
    for di, df_x in enumerate([df_counts_vb_curr, df_counts_jb_curr, df_counts_va_curr, df_counts_ja_curr]):

        # Normalize segment frequencies
        data.append(df_x.copy().reindex(order, columns=cols[di]).fillna(
            0).apply(lambda x: x / x.sum(), axis=1))

    seg_pca = pca.fit_transform(pd.concat(data, axis=1))

    # ratio
    print(f"Explained variance ratios for {plotname}:")
    for i, ratio in enumerate(pca.explained_variance_ratio_):
        print(f"  PC{i+1}: {ratio:.4f}")
    print(f"  Cumulative: {pca.explained_variance_ratio_.cumsum()[-1]:.4f}\n")

    df_seg_pca = pd.DataFrame(seg_pca, index=order)

    fig, axs = plt.subplots(1, 1, figsize=(3, 3))

    for subname in order:

        axs.scatter(df_seg_pca.loc[subname, 1], df_seg_pca.loc[subname, 2],
                    marker=get_marker(subname), s=40, color=get_color(subname))

    # axs.set(xlabel = 'Principal Component 1', ylabel = 'Principal Component 3')

    plt.legend(order, loc='center left', bbox_to_anchor=(1.05, 0.5))

    axs.set_xlabel('PC2')
    axs.set_ylabel('PC3')
    # Sets area outside ax to transparent
    fig.savefig(
        PATH_FIG_DIR / f'Fig.F - PCA.pdf', bbox_inches='tight', facecolor='#ffffff00')
    plt.show()

# %% [markdown]
# #### Densities
# 
# S: A smoothing condition. The amount of smoothness is determined by satisfying the conditions: sum((w * (y - g))**2,axis=0) <= s, where g(x) is the smoothed interpolation of (x,y). The user can use s to control the trade-off between closeness and smoothness of fit. Larger s means more smoothing while smaller values of s indicate less smoothing. Recommended values of s depend on the weights, w. If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range (m-sqrt(2*m),m+sqrt(2*m)), where m is the number of data points in x, y, and w.

# %%
ylims = [(0, 0.06), (0, 0.25), (0, 0.085), (0, 0.055)]
cols = [TRBV, TRBJ, TRAV, TRAJ]

for plotname, plotinfo in plots.items():

    order = orders[plotname]

    df_counts_va_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VA', aggfunc='sum', fill_value=0)
    df_counts_ja_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JA', aggfunc='sum', fill_value=0)
    df_counts_vb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='VB', aggfunc='sum', fill_value=0)
    df_counts_jb_curr = df.pivot_table(
        values='Templates', index=plotinfo['Subset'], columns='JB', aggfunc='sum', fill_value=0)

    f, axs = plt.subplots(1, 4, figsize=(13, 3), width_ratios=[1, 0.5, 1, 1])
    # f.tight_layout(w_pad = 0)

    for di, df_y in enumerate([df_counts_vb_curr, df_counts_jb_curr, df_counts_va_curr, df_counts_ja_curr]):

        ax = axs[di]
        data = df_y.copy().reindex(order, columns=cols[di]).fillna(
            0).apply(lambda x: x / x.sum(), axis=1)

        for subname in order:

            # if ('immature' in subname) or ('(P)' in subname):
            #     linestyle = 'dotted' # overwrite for these to separate 4_imm vs 4_mat and DP vs DP (P)

            y = data.loc[subname, :]
            w = [1 / np.std(y)] * len(y)
            x = np.arange(len(y))
            # , s = len(y) - np.sqrt(2*len(y)))
            tck = splrep(x, y, w, s=len(y))
            y_spline = BSpline(*tck)(x)

            ax.plot(x, y_spline, linewidth=2, linestyle=get_linestyle(
                subname), color=get_color(subname))

        ax.set(xlabel='', ylabel='', xticks=[], ylim=ylims[di])
        # ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)

    plt.legend(order, loc='center left', bbox_to_anchor=(1.05, 0.5))

    # Sets area outside ax to transparent
    f.savefig(PATH_FIG_DIR /
              f'Fig.E - Position density.pdf', bbox_inches='tight', facecolor='#ffffff00')
    plt.show()

# %% [markdown]
# #### Va + Vb

# %%
plotname = 'Simple'
order = orders[plotname]

figW = [10, 12]
figH = 1.7
cbar_pos = [.92, 0.1, .01, 0.8]

cols = [TRAV, TRBV]
labels = ['Va', 'Vb']

for di, df_v in enumerate([df_counts_va, df_counts_vb]):

    data = df_v.copy().reindex(order, columns=cols[di]).fillna(
        0).apply(lambda x: x / x.sum(), axis=1)

    # Fractional counts per sample
    f, ax = plt.subplots(figsize=(figW[di], figH), ncols=1, nrows=1)
    # cbar_ax = f.add_axes(cbar_pos)

    sns.heatmap(ax=ax, data=data, annot=False, linewidths=0.0, fmt='.3g', vmin=0, vmax=0.15,
                cmap='viridis', cbar=False)  # cbar_ax = cbar_ax,
    ax.set(xlabel='', ylabel='')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # Sets area outside ax to transparent
    f.savefig(PATH_FIG_DIR / f'Fig.C - {labels[di]} {plotname}.pdf',
              bbox_inches='tight', facecolor='#ffffff00')

    # Z-score by column
    f, ax = plt.subplots(figsize=(figW[di], figH), ncols=1, nrows=1)
    cbar_ax = f.add_axes(cbar_pos)

    sns.heatmap(ax=ax, data=data.apply(stats.zscore, axis=0), annot=False, linewidths=0.0, fmt='.3g', vmin=-2, vmax=2,
                cbar_ax=cbar_ax, cmap='RdBu_r')
    ax.set(xlabel='', ylabel='')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # Z-score by row
    # f, ax = plt.subplots(figsize = (figW[di], figH), ncols = 1, nrows = 1)
    # cbar_ax = f.add_axes(cbar_pos)

    # sns.heatmap(ax = ax, data = data.apply(stats.zscore, axis = 1), annot = False, linewidths = 0.0, fmt = '.3g', vmin = -2, vmax = 2, \
    #             cbar_ax = cbar_ax, cmap = 'RdBu_r')
    # ax.set(xlabel = '', ylabel = '')
    # ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)


plt.show()

# %% [markdown]
# #### Ja + Jb

# %%
plotname = 'Simple'
order = orders[plotname]

figW = [11, 3]
figH = 1.7
cbar_pos = [[.92, 0.1, .01, 0.8],
            [.96, 0.1, .04, 0.8]]

cols = [TRAJ, TRBJ]
labels = ['Ja', 'Jb']

for di, df_j in enumerate([df_counts_ja, df_counts_jb]):

    data = df_j.copy().reindex(order, columns=cols[di]).fillna(
        0).apply(lambda x: x / x.sum(), axis=1)

    # Fractional counts per sample
    f, ax = plt.subplots(figsize=(figW[di], figH), ncols=1, nrows=1)
    cbar_ax = f.add_axes(cbar_pos[di])

    sns.heatmap(ax=ax, data=data, annot=False, linewidths=0.0, fmt='.3g', vmin=0, vmax=0.15,
                cbar_ax=cbar_ax, cmap='viridis')
    ax.set(xlabel='', ylabel='', yticklabels=[])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    # Sets area outside ax to transparent
    f.savefig(PATH_FIG_DIR / f'Fig.D - {labels[di]} {plotname}.pdf',
              bbox_inches='tight', facecolor='#ffffff00')

    # Z-score by column
    f, ax = plt.subplots(figsize=(figW[di], figH), ncols=1, nrows=1)
    cbar_ax = f.add_axes(cbar_pos[di])

    sns.heatmap(ax=ax, data=data.apply(stats.zscore, axis=0), annot=False, linewidths=0.0, fmt='.3g', vmin=-2, vmax=2,
                cbar_ax=cbar_ax, cmap='RdBu_r')
    ax.set(xlabel='', ylabel='')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

plt.show()




print("Original full TCR-seq analysis completed.")
print(f"Output folder: {PATH_FIG_DIR}")
