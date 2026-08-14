"""
RAG Heatmap Generator - Parameterized Version with Polymer Markers

This script creates a 3D visualization showing the best achievable RAG score
for dissolving substrates at any point in Hansen parameter space.

Modified to include markers for specific polymers/substances (In this instance DCM, NMP, PEEK, PES. Can be customised.)
with smart label positioning that adapts to viewing angle.

Adjust the parameters in the CONFIGURATION section below to customize your analysis.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import os
from datetime import datetime

################################################################################
#                              CONFIGURATION
################################################################################

# Input Excel file path (relative to script location or absolute)
EXCEL_FILE = 'Parameter_Set-Up_Pharma.xlsx'

# Grid parameters
GRID_SPACING = 0.5          # Grid spacing in Hansen parameter units (0.5 = high res, 1.0 = medium res, 2.0 = coarse)
DD_RANGE = (14, 22)         # dD range: (min, max)
DP_RANGE = (0, 20)          # dP range: (min, max)
DH_RANGE = (0, 20)          # dH range: (min, max)

# RED calculation parameters
RED_THRESHOLD = 0.3         # Only consider solvents with RED < this value
R_0 = 8.0                   # Interaction radius for RED calculation

# Solvent mixture parameters
MIXING_STEP = 5            # Volume % steps for binary mixtures (5 = 5%, 10%, 15%, 20%, etc.)

# Visualization parameters
POINT_SIZE = 8              # Size of points in scatter plot
POINT_ALPHA = 0.25          # Transparency of points (0.0 - 1.0)
PLOT_DPI = 200              # Resolution of saved image

# Polymer marker parameters
POLYMER_MARKER_SIZE = 30    # Size of polymer markers
POLYMER_MARKER_STYLE = '*'  # Marker style (*, o, s, D, ^, etc.)
POLYMER_LABEL_FONTSIZE = 9  # Font size for polymer labels
SHOW_LABELS = False          # Set to False to hide all labels (legend only)

# Output
OUTPUT_FILENAME = 'RAG_Heatmap.png'  # Will be saved to results/ folder

################################################################################
#                           END CONFIGURATION
################################################################################


# Determine script directory
script_dir = os.path.dirname(os.path.abspath(__file__))
excel_path = os.path.join(script_dir, EXCEL_FILE)
output_dir = os.path.join(script_dir, 'results')

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

# Read the Excel file
print(f"\n{'='*70}")
print(f"RAG HEATMAP GENERATOR WITH POLYMER/SOLVENT MARKERS")
print(f"{'='*70}\n")

# Read polymer data from first sheet
print(f"Reading polymer/substance data from: {EXCEL_FILE} (Polymer sheet)")
polymer_data = pd.read_excel(excel_path, sheet_name='Polymer', header=None)

# Extract polymer parameters (columns 1-5, rows 3-7)
polymers = {}
for col in range(1, 6):  # Columns 1-5 (0-indexed: 1-5)
    polymer_name = polymer_data.iloc[3, col]  # Row 3 has names
    if pd.notna(polymer_name) and polymer_name != '':
        dD = polymer_data.iloc[4, col]  # Row 4 has dD_B
        dP = polymer_data.iloc[5, col]  # Row 5 has dP_B
        dH = polymer_data.iloc[6, col]  # Row 6 has dH_B
        
        if pd.notna(dD) and pd.notna(dP) and pd.notna(dH):
            try:
                polymers[str(polymer_name)] = {
                    'dD': float(dD),
                    'dP': float(dP),
                    'dH': float(dH)
                }
            except (ValueError, TypeError):
                continue

print(f"✓ Loaded {len(polymers)} polymers/substances:")
for name, params in polymers.items():
    print(f"   - {name}: δD={params['dD']}, δP={params['dP']}, δH={params['dH']}")

# Read the solvent data from Excel
print(f"\nReading solvent data from: {EXCEL_FILE} (Solvent sheet)")

sdata = pd.read_excel(excel_path, sheet_name='Solvent', header=None)

# Extract solvent data
allvalues = []
rag_dict = {}

for index, row in sdata.iterrows():
    if index > 0:  # Skip header row
        solvent_name = row[0]
        try:
            dD = float(row[1])
            dP = float(row[2])
            dH = float(row[3])
            
            rag_value = row[4] if len(row) > 4 else None
            if pd.notna(rag_value):
                rag_dict[solvent_name] = str(rag_value).strip()
            else:
                rag_dict[solvent_name] = ''
            
            allvalues.append((solvent_name, dD, dP, dH))
        except (ValueError, TypeError):
            continue

print(f"✓ Loaded {len(allvalues)} solvents")
print(f"✓ RAG scores assigned: {sum(1 for v in rag_dict.values() if v)}")

# Print configuration
print(f"\n📋 Configuration:")
print(f"   Grid spacing: {GRID_SPACING} units")
print(f"   dD range: {DD_RANGE[0]} to {DD_RANGE[1]}")
print(f"   dP range: {DP_RANGE[0]} to {DP_RANGE[1]}")
print(f"   dH range: {DH_RANGE[0]} to {DH_RANGE[1]}")
print(f"   RED threshold: < {RED_THRESHOLD}")
print(f"   Mixing step: {MIXING_STEP}%")

# Convert to numpy arrays for vectorized operations
n_solvents = len(allvalues)
solvent_names = [s[0] for s in allvalues]
solvent_params = np.array([[s[1], s[2], s[3]] for s in allvalues])

# RAG priority (lower = better/safer)
RAG_PRIORITY = {'G': 0, 'A': 1, 'R': 2, 'Grey': 3, '': 4}

print("\n🔧 Pre-computing solvent mixture library...")

# Pre-compute all possible mixtures
mixing_ratios = list(range(MIXING_STEP, 100, MIXING_STEP))
mixture_library = []

# Add pure solvents
for i, (name, dD, dP, dH) in enumerate(allvalues):
    rag = rag_dict.get(name, 'Grey')
    if rag == '':
        rag = 'Grey'
    mixture_library.append({
        'params': np.array([dD, dP, dH]),
        'rag': rag,
        'name': f"{name} (100%)"
    })

# Add binary mixtures
for i in range(n_solvents):
    for j in range(i + 1, n_solvents):
        name1, dD1, dP1, dH1 = allvalues[i]
        name2, dD2, dP2, dH2 = allvalues[j]
        
        rag1 = rag_dict.get(name1, '')
        rag2 = rag_dict.get(name2, '')
        
        # Determine mixture RAG (worst of the two)
        if rag1 == '' and rag2 == '':
            mix_rag = 'Grey'
        elif rag1 == '':
            mix_rag = rag2
        elif rag2 == '':
            mix_rag = rag1
        else:
            priority_order = {'R': 0, 'A': 1, 'Grey': 2, 'G': 3}
            p1 = priority_order.get(rag1, 99)
            p2 = priority_order.get(rag2, 99)
            mix_rag = rag1 if p1 <= p2 else rag2
        
        for ratio in mixing_ratios:
            vol1 = ratio
            vol2 = 100 - ratio
            dD_mix = (vol1 * dD1 + vol2 * dD2) / 100
            dP_mix = (vol1 * dP1 + vol2 * dP2) / 100
            dH_mix = (vol1 * dH1 + vol2 * dH2) / 100
            
            mixture_library.append({
                'params': np.array([dD_mix, dP_mix, dH_mix]),
                'rag': mix_rag,
                'name': f"{name1} ({vol1}%) + {name2} ({vol2}%)"
            })

print(f"   ✓ Created library of {len(mixture_library):,} solvent options")
print(f"     ({n_solvents} pure + {len(mixture_library) - n_solvents:,} mixtures)")

# Convert mixture library to numpy arrays for fast computation
mixture_params = np.array([m['params'] for m in mixture_library])
mixture_rags = [m['rag'] for m in mixture_library]

print("\n🔬 Creating Hansen space grid...")
dD_range = np.arange(DD_RANGE[0], DD_RANGE[1], GRID_SPACING)
dP_range = np.arange(DP_RANGE[0], DP_RANGE[1], GRID_SPACING)
dH_range = np.arange(DH_RANGE[0], DH_RANGE[1], GRID_SPACING)

# Create grid as vectors
dD_grid, dP_grid, dH_grid = np.meshgrid(dD_range, dP_range, dH_range, indexing='ij')
grid_points_flat = np.stack([dD_grid.ravel(), dP_grid.ravel(), dH_grid.ravel()], axis=1)

total_points = len(grid_points_flat)
print(f"   Grid dimensions: {len(dD_range)} × {len(dP_range)} × {len(dH_range)} = {total_points:,} points")

print(f"\n⚙️  Evaluating best achievable RAG at each grid point...")
print(f"   Started at: {datetime.now().strftime('%H:%M:%S')}")

start_time = datetime.now()

# Process in batches for memory efficiency
batch_size = 200
results = []

for batch_start in range(0, total_points, batch_size):
    batch_end = min(batch_start + batch_size, total_points)
    batch_targets = grid_points_flat[batch_start:batch_end]
    
    # Progress reporting
    if batch_start % (batch_size * 10) == 0:
        elapsed = (datetime.now() - start_time).total_seconds()
        progress_pct = 100 * batch_start / total_points
        rate = batch_start / elapsed if elapsed > 0 else 0
        remaining = (total_points - batch_start) / rate / 60 if rate > 0 else 0
        print(f"   Progress: {batch_start:6d}/{total_points} ({progress_pct:5.1f}%) | "
              f"Rate: {rate:.0f} pts/s | ETA: {remaining:.1f} min")
    
    # For each target in the batch
    for target in batch_targets:
        # Calculate RED for ALL mixtures to this target (vectorized)
        diff = mixture_params - target
        diff[:, 0] *= 2  # dD has factor of 4 under sqrt, so *2 outside
        red_values = np.sqrt(np.sum(diff**2, axis=1)) / R_0
        
        # Find mixtures with RED < threshold
        valid_mask = red_values < RED_THRESHOLD
        
        if np.any(valid_mask):
            valid_rags = [mixture_rags[i] for i in np.where(valid_mask)[0]]
            
            # Find best (safest) RAG
            best_rag = min(valid_rags, key=lambda r: RAG_PRIORITY.get(r, 99))
            
            results.append({
                'dD': target[0],
                'dP': target[1],
                'dH': target[2],
                'rag': best_rag
            })

elapsed_total = (datetime.now() - start_time).total_seconds()
print(f"\n✅ Evaluation complete in {elapsed_total:.1f} seconds ({elapsed_total/60:.2f} minutes)")
print(f"   {len(results):,} grid points have achievable solutions (RED < {RED_THRESHOLD})")

# Separate by RAG category
rag_categories = {
    'G': {'color': 'green', 'label': 'Green achievable', 'points': []},
    'A': {'color': 'orange', 'label': 'Amber achievable', 'points': []},
    'R': {'color': 'red', 'label': 'Red achievable', 'points': []},
    'Grey': {'color': 'gray', 'label': 'Grey only', 'points': []}
}

for point in results:
    rag = point['rag']
    if rag in rag_categories:
        rag_categories[rag]['points'].append(point)

# Print summary
print(f"\n📊 Coverage Summary:")
for rag, data in rag_categories.items():
    count = len(data['points'])
    if count > 0:
        pct = 100 * count / total_points
        print(f"   {data['label']:20s}: {count:6,} grid points ({pct:5.1f}% of total)")

no_solution = total_points - len(results)
print(f"   {'No solution':20s}: {no_solution:6,} grid points ({100*no_solution/total_points:5.1f}% of total)")

# Create 3D visualization
print(f"\n🎨 Generating 3D visualization with polymer markers...")

fig = plt.figure(figsize=(20, 14))

views = [
    {'elev': 20, 'azim': 45, 'title': 'View 1 (Front-Right)'},
    {'elev': 20, 'azim': 135, 'title': 'View 2 (Back-Right)'},
    {'elev': 20, 'azim': 225, 'title': 'View 3 (Back-Left)'},
    {'elev': 60, 'azim': 45, 'title': 'View 4 (Top View)'}
]

# Define shapes for different polymers/substances (all black)
polymer_markers = {
    'DCM': 'o',      # Circle
    'NMP': 's',      # Square
    'PEEK': '^',     # Triangle up
    'PES': 'D'       # Diamond
}

# Define custom label offsets for each view to avoid overlaps
# Structure: {view_index: {polymer_name: (dx, dy, dz)}}
label_offsets = {
    0: {  # View 1 (Front-Right) - elev=20, azim=45
        'DCM': (-1.5, 0.5, 0.5),
        'NMP': (0.8, 1.0, 0.3),
        'PEEK': (0.8, -1.2, 0.0),
        'PES': (0.8, 0.8, 0.8)
    },
    1: {  # View 2 (Back-Right) - elev=20, azim=135
        'DCM': (0.8, 0.8, 0.5),
        'NMP': (-1.5, 0.8, 0.3),
        'PEEK': (0.8, -1.2, 0.0),
        'PES': (-1.5, 1.0, 0.8)
    },
    2: {  # View 3 (Back-Left) - elev=20, azim=225
        'DCM': (0.8, -1.2, 0.5),
        'NMP': (0.8, 0.8, 0.3),
        'PEEK': (-1.5, 0.8, 0.0),
        'PES': (0.8, -1.2, 0.8)
    },
    3: {  # View 4 (Top View) - elev=60, azim=45
        'DCM': (-1.2, -1.0, 0.8),
        'NMP': (0.8, 1.2, 0.5),
        'PEEK': (0.8, -1.2, 0.0),
        'PES': (1.0, 1.0, 1.0)
    }
}

for idx, view in enumerate(views, 1):
    ax = fig.add_subplot(2, 2, idx, projection='3d')
    
    # Plot RAG heatmap points
    rag_handles = []
    rag_labels = []
    for rag, data in rag_categories.items():
        if len(data['points']) > 0:
            points = data['points']
            dD_vals = [p['dD'] for p in points]
            dP_vals = [p['dP'] for p in points]
            dH_vals = [p['dH'] for p in points]
            
            scatter = ax.scatter(dD_vals, dP_vals, dH_vals,
                      c=data['color'],
                      alpha=POINT_ALPHA,
                      s=POINT_SIZE,
                      edgecolors='none')
            rag_handles.append(scatter)
            rag_labels.append(f"{data['label']} ({len(points):,})")
    
    # Plot polymer markers
    polymer_handles = []
    polymer_labels = []
    for polymer_name, params in polymers.items():
        dD = params['dD']
        dP = params['dP']
        dH = params['dH']
        
        # Get marker shape for this polymer
        marker_shape = polymer_markers.get(polymer_name, 'X')
        
        # Plot marker (all black)
        scatter = ax.scatter([dD], [dP], [dH],
                  c='black',
                  marker=marker_shape,
                  s=POLYMER_MARKER_SIZE,
                  edgecolors='black',
                  linewidths=0.8,
                  zorder=1000)  # Ensure markers are on top
        
        polymer_handles.append(scatter)
        polymer_labels.append(polymer_name)
        
        # Add label with custom offset for this view (if enabled)
        if SHOW_LABELS:
            view_idx = idx - 1
            offset = label_offsets[view_idx].get(polymer_name, (0.8, 0.8, 0.8))
            
            ax.text(dD + offset[0], 
                    dP + offset[1], 
                    dH + offset[2],
                    polymer_name,
                    fontsize=POLYMER_LABEL_FONTSIZE,
                    fontweight='bold',
                    color='black',
                    bbox=dict(boxstyle='round,pad=0.3', 
                             facecolor='white', 
                             edgecolor='black',
                             alpha=0.8))
    
    ax.set_xlabel('δD (Dispersion)', fontsize=10, fontweight='bold')
    ax.set_ylabel('δP (Polar)', fontsize=10, fontweight='bold')
    ax.set_zlabel('δH (Hydrogen bonding)', fontsize=10, fontweight='bold')
    ax.set_title(view['title'], fontsize=12, fontweight='bold', pad=10)
    ax.view_init(elev=view['elev'], azim=view['azim'])
    
    # Create legend with two separate columns: RAG categories and Polymers
    # First legend: RAG categories
    if rag_handles:
        legend1 = ax.legend(rag_handles, rag_labels, 
                           loc='upper left', fontsize=8, 
                           title='RAG Categories',
                           title_fontsize=9,
                           framealpha=0.9)
        ax.add_artist(legend1)  # Add first legend back
    
    # Second legend: Polymers/Substances
    if polymer_handles:
        ax.legend(polymer_handles, polymer_labels,
                 loc='upper right', fontsize=8,
                 title='Solvents/Polymers',
                 title_fontsize=9,
                 framealpha=0.9)
    
    ax.grid(True, alpha=0.3)
    
    # Set axis limits to include polymer positions
    all_dD = [params['dD'] for params in polymers.values()]
    all_dP = [params['dP'] for params in polymers.values()]
    all_dH = [params['dH'] for params in polymers.values()]
    
    if all_dD:
        ax.set_xlim([min(DD_RANGE[0], min(all_dD)-1), max(DD_RANGE[1], max(all_dD)+1)])
        ax.set_ylim([min(DP_RANGE[0], min(all_dP)-1), max(DP_RANGE[1], max(all_dP)+1)])
        ax.set_zlim([min(DH_RANGE[0], min(all_dH)-1), max(DH_RANGE[1], max(all_dH)+1)])

label_status = "with labels" if SHOW_LABELS else "legend only"
plt.suptitle(f'3D Pharmaceuticals Guide RAG Heatmap - Grid: {GRID_SPACING} units, RED < {RED_THRESHOLD}\n'
             f'Best Achievable RAG Score | {total_points:,} points evaluated',
             fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()

output_path = os.path.join(output_dir, OUTPUT_FILENAME)
plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches='tight', pad_inches=0.4)
print(f"   ✓ Saved to: {output_path}")

plt.show()

print(f"\n{'='*70}")
print(f"✅ RAG HEATMAP WITH POLYMER MARKERS COMPLETE")
print(f"{'='*70}")
print(f"Completed at: {datetime.now().strftime('%H:%M:%S')}")
print(f"\n💡 Interpretation:")
print(f"   🟢 Green regions  = Safe solvents available (RED < {RED_THRESHOLD})")
print(f"   🟠 Amber regions  = Only caution solvents available")
print(f"   🔴 Red regions    = Only hazard solvents available")
print(f"   ⚫ Grey regions   = Only unrated solvents available")
print(f"   (Empty space)    = No suitable solvent mixtures exist")
print(f"\n   ⭐ Black markers show polymer/substance positions:")
print(f"   DCM:  ● (circle)")
print(f"   NMP:  ■ (square)")
print(f"   PEEK: ▲ (triangle)")
print(f"   PES:  ◆ (diamond)")
print()
