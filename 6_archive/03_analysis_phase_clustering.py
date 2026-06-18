"""
Pipeline: Hydrogen Bond Frequency Analysis & Phase Clustering
Simulations: V7, V12 (extensible to more)

Steps:
    1. Load interaction CSVs
    2. Slice transition-state frames
    3. Split into phase 1 / phase 2 by phase_limit
    4. Compute per-phase H-bond frequencies
    5. Cluster phases by frequency profile (dendrogram)
    6. Optional: cluster by delta-frequency (ph2 - ph1) to focus on phase *change*
    7. PCA visualization
    8. Decision Tree for interpretability
"""

# ── Imports ────────────────────────────────────────────────────────────────────
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, export_graphviz, plot_tree
from io import StringIO

# ── Configuration ──────────────────────────────────────────────────────────────
ROOT_DIR = '/content/drive/MyDrive/M1_STAGE/Data/interactions/'   # adapt to your env
OUTPUT_DIR = '/content/drive/MyDrive/M1_STAGE/Manips/Tables/'     # adapt to your env
os.makedirs(OUTPUT_DIR, exist_ok=True)

INTERACTION_FILES = ['res_V12.csv', 'res_V7.csv']
SIMULATION_NAMES  = ['V12', 'V7']

# Frames of interest per simulation (transition-state window)
SELECTED_FRAMES = {
    'V12': range(0,   175 + 1),
    'V7':  range(75,  250 + 1),
    # 'V1':  range(412, 587 + 1),
    # 'V8':  range(62,  237 + 1),
    # 'V21': range(62,  237 + 1),
    # 'V11': range(698, 873 + 1),
}

# Index threshold separating phase 1 from phase 2
PHASE_LIMITS = {
    'V12': 35,
    'V7':  198,
    # 'V1':  517,
    # 'V8':  180,
    # 'V21': 155,
    # 'V11': 1000,
}

# Desired column order in the frequency table
PHASES_ORDERED = ['V12_ph1', 'V12_ph2', 'V7_ph1', 'V7_ph2']

# ── Helpers ────────────────────────────────────────────────────────────────────

def find_hbonds(df: pd.DataFrame) -> list[str]:
    """Return column names that correspond to hydrogen bonds (contain 'hb')."""
    return [col for col in df.columns if 'hb' in col]


def phase_frequency(df_phase: pd.DataFrame, hbond_cols: list[str]) -> dict:
    """
    Compute the occurrence frequency (%) of each H-bond in a phase DataFrame.
    Returns 0.0 for empty phases or missing columns.
    """
    n = len(df_phase)
    freqs = {}
    for col in hbond_cols:
        if col in df_phase.columns and n > 0:
            freqs[col] = (df_phase[col] > 0).sum() / n * 100
        else:
            freqs[col] = 0.0
    return freqs


def slice_transition_state(df_full: pd.DataFrame, sim_name: str) -> pd.DataFrame:
    """Return only the transition-state rows defined in SELECTED_FRAMES."""
    if sim_name not in SELECTED_FRAMES:
        print(f"[WARN] No frame range for {sim_name}. Using full DataFrame.")
        return df_full
    valid = [i for i in df_full.index if i in SELECTED_FRAMES[sim_name]]
    return df_full.loc[valid]


def split_phases(df_ts: pd.DataFrame, sim_name: str):
    """
    Split transition-state DataFrame into phase 1 and phase 2.
    Returns (df_ph1, df_ph2) or raises if no phase limit is defined.
    """
    if sim_name not in PHASE_LIMITS:
        raise ValueError(f"No phase limit defined for {sim_name}.")
    lim = PHASE_LIMITS[sim_name]
    idx = df_ts.index
    ph1 = df_ts.loc[idx[idx <= lim]] if (idx <= lim).any() else pd.DataFrame(columns=df_ts.columns)
    ph2 = df_ts.loc[idx[idx >  lim]] if (idx >  lim).any() else pd.DataFrame(columns=df_ts.columns)
    return ph1, ph2

# ── 1–4: Build frequency table ─────────────────────────────────────────────────
os.chdir(ROOT_DIR)

all_phase_frequencies = {}

for sim_file, sim_name in zip(INTERACTION_FILES, SIMULATION_NAMES):
    df_full   = pd.read_csv(sim_file)
    df_ts     = slice_transition_state(df_full, sim_name)

    if df_ts.empty:
        print(f"[WARN] Empty transition state for {sim_name}. Skipping.")
        continue

    df_ph1, df_ph2 = split_phases(df_ts, sim_name)
    hbonds = find_hbonds(df_ts)          # same bond universe for both phases

    all_phase_frequencies[f"{sim_name}_ph1"] = phase_frequency(df_ph1, hbonds)
    all_phase_frequencies[f"{sim_name}_ph2"] = phase_frequency(df_ph2, hbonds)

    print(f"{sim_name}: {len(df_ph1)} frames in ph1, {len(df_ph2)} frames in ph2, "
          f"{len(hbonds)} H-bonds found.")

# Assemble and order columns
freq_table = pd.DataFrame(all_phase_frequencies)
existing_ordered = [c for c in PHASES_ORDERED if c in freq_table.columns]
freq_table = freq_table[existing_ordered]

# Drop bonds that are zero everywhere
freq_table = freq_table[freq_table.sum(axis=1) > 0]

# Save full table
freq_table.to_csv(os.path.join(OUTPUT_DIR, 'freq_table_all.csv'),
                  index=True, decimal='.', float_format='%.5f')
print(f"\nFrequency table: {freq_table.shape[0]} bonds × {freq_table.shape[1]} phases")

# ── 5: Dendrogram – cluster PHASES by raw frequency profile ───────────────────
# X rows = phases, columns = bond frequencies  →  4 leaves (one per phase)
# This is the correct interpretation: "how similar are the four phases?"

filtered_freq = freq_table[(freq_table > 80).any(axis=1)]   # keep bonds with ≥80% in at least one phase
print(f"Bonds retained after 80% filter: {len(filtered_freq)}")

X_phases = filtered_freq.T           # shape: (n_phases, n_bonds)
y_phases = X_phases.index.tolist()   # phase labels

linked_phases = linkage(X_phases, method='complete')

fig, ax = plt.subplots(figsize=(8, 5))
dendrogram(linked_phases, orientation='top', labels=y_phases,
           distance_sort='descending', ax=ax)
ax.set_title('Phase clustering – raw H-bond frequencies (≥80% filter)')
ax.set_ylabel('Distance')
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'dendrogram_phases_raw.png'), dpi=300)
plt.show()

# ── 6: Dendrogram – cluster PHASES by Δ-frequency (ph2 − ph1) ─────────────────
# This removes simulation-identity signal and focuses on the *change* between phases,
# which is more likely to reveal a common ph1→ph2 transition signature.

delta_rows = {}
for sim_name in SIMULATION_NAMES:
    ph1_col = f"{sim_name}_ph1"
    ph2_col = f"{sim_name}_ph2"
    if ph1_col in freq_table.columns and ph2_col in freq_table.columns:
        delta_rows[sim_name] = freq_table[ph2_col] - freq_table[ph1_col]

delta_df = pd.DataFrame(delta_rows).T   # shape: (n_sims, n_bonds)
delta_df = delta_df.loc[:, (delta_df != 0).any()]   # drop bonds that never change

linked_delta = linkage(delta_df, method='ward')

fig, ax = plt.subplots(figsize=(8, 5))
dendrogram(linked_delta, orientation='top', labels=delta_df.index.tolist(),
           distance_sort='descending', ax=ax)
ax.set_title('Simulation clustering – Δ frequency (ph2 − ph1)\n'
             'Similar = similar transition mechanism')
ax.set_ylabel('Distance')
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'dendrogram_delta_freq.png'), dpi=300)
plt.show()

# ── 6b: Dendrogram – cluster BONDS by frequency profile across phases ──────────
# X rows = bonds, columns = phases  →  many leaves (one per bond)
# Useful to visually group bonds that behave similarly across conditions.

X_bonds = filtered_freq          # shape: (n_bonds, n_phases)   no transpose
linked_bonds = linkage(X_bonds, method='ward')

fig, ax = plt.subplots(figsize=(max(10, len(filtered_freq) // 3), 7))
dendrogram(linked_bonds, orientation='top', labels=X_bonds.index.tolist(),
           distance_sort='descending', leaf_rotation=90, leaf_font_size=7, ax=ax)
ax.set_title('Bond clustering – frequency profile across phases')
ax.set_ylabel('Distance')
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'dendrogram_bonds.png'), dpi=300)
plt.show()

# ── 7: PCA ────────────────────────────────────────────────────────────────────
# Load per-phase frame-level interaction CSVs for frame-wise PCA

phase_csv_map = {
    'V7_ph1':  'V7_ph1_interactions.csv',
    'V7_ph2':  'V7_ph2_interactions.csv',
    'V12_ph1': 'V12_ph1_interactions.csv',
    'V12_ph2': 'V12_ph2_interactions.csv',
}

os.chdir(OUTPUT_DIR)    # phase CSVs are stored in the tables dir

df_pca_input = pd.DataFrame()
for label, fname in phase_csv_map.items():
    if not os.path.exists(fname):
        print(f"[WARN] {fname} not found, skipping for PCA.")
        continue
    df_tmp = pd.read_csv(fname)
    hb_cols = find_hbonds(df_tmp)
    df_tmp  = df_tmp[hb_cols].copy()
    df_tmp['label'] = label
    df_pca_input = pd.concat([df_pca_input, df_tmp], ignore_index=True)

if not df_pca_input.empty:
    # Drop all-zero rows and columns
    df_pca_input = df_pca_input.loc[(df_pca_input.drop('label', axis=1) != 0).any(axis=1)]
    zero_cols = [c for c in df_pca_input.columns
                 if c != 'label' and df_pca_input[c].sum() == 0]
    df_pca_input.drop(columns=zero_cols, inplace=True)

    X_pca = df_pca_input.drop('label', axis=1)
    y_pca = df_pca_input['label']

    scaler     = StandardScaler()
    X_scaled   = scaler.fit_transform(X_pca)

    pca        = PCA(n_components=2)
    components = pca.fit_transform(X_scaled)
    var_expl   = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame(components,
                          columns=['PC1', 'PC2'])
    pca_df['label'] = y_pca.values

    fig = px.scatter(pca_df, x='PC1', y='PC2', color='label',
                     title=f'PCA of frame-level H-bond interactions<br>'
                           f'PC1 {var_expl[0]:.1f}% | PC2 {var_expl[1]:.1f}%',
                     labels={'PC1': f'PC1 ({var_expl[0]:.1f}%)',
                             'PC2': f'PC2 ({var_expl[1]:.1f}%)'})
    fig.update_traces(marker=dict(size=6, opacity=0.6))
    fig.write_html(os.path.join(OUTPUT_DIR, 'pca_phases.html'))
    fig.show()
else:
    print("[INFO] No phase CSVs found; PCA skipped.")

# ── 8: Decision Tree ──────────────────────────────────────────────────────────
# Fit on the aggregated frequency table (phases as samples, bonds as features)
# to find the bonds that best discriminate between sim-phases.

X_dt = X_phases                       # (n_phases, n_bonds)
y_dt = y_phases                       # phase labels

clf = DecisionTreeClassifier(max_depth=4, random_state=42)
clf.fit(X_dt, y_dt)

fig, ax = plt.subplots(figsize=(16, 6))
plot_tree(clf, feature_names=X_phases.columns.tolist(),
          class_names=y_dt, filled=True, rounded=True,
          impurity=False, proportion=False, ax=ax)
ax.set_title('Decision Tree – discriminating bonds between sim-phases')
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'decision_tree.png'), dpi=300, bbox_inches='tight')
plt.show()

print("\nMost discriminating features (top 10 by importance):")
importances = pd.Series(clf.feature_importances_, index=X_phases.columns)
print(importances.sort_values(ascending=False).head(10))
