"""
Circular solvent pairs as replacement options diagram generator

Creates a circular diagram showing solvent pairs with:
- Each solvent as a coloured point around a circle
- Lines connecting ALL pairs of solvents (no deduplication)
- Crosses (×) with plain black letter labels (a, b, c...) at mixture points
- Compact legend showing letter-to-RED value mapping
- Solvent labels with colored borders and white backgrounds
- Title from Excel row 1

Input: Excel file with columns [RED, Solvent1, Vol(%), Solvent2, Vol(%)]
Output: PNG image of the circular network diagram showing all pairs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

################################################################################
#                              CONFIGURATION
################################################################################

# Input file
INPUT_FILE = 'table1_DCM_P&C.xlsx'   # Path to your Excel file (create after running SolventFinder)
SKIP_ROWS = 2                        # Rows to skip before data starts

# Display settings
N_SOLVENTS = 18                      # Number of most-connected solvents to show

# Visual styling
LINE_WIDTH = 2.5                     # Thickness of connection lines
LINE_ALPHA = 0.6                     # Transparency of lines (0-1)
POINT_SIZE = 450                     # Size of solvent circles

# Cross and letter styling
CROSS_SIZE = 10                      # Size of × markers
CROSS_WIDTH = 2.5                    # Line width of × markers
LETTER_FONTSIZE = 11.25              # Font size for letter labels 
LETTER_OFFSET_X = 0.05               # X offset for letter from cross
LETTER_OFFSET_Y = 0.05               # Y offset for letter from cross

# Legend styling
LEGEND_LINE_SPACING = 0.04           # Vertical spacing between legend entries
LEGEND_FONTSIZE = 9                  # Font size for legend entries

# Output
OUTPUT_FILE = 'circle_diagram.png'
DPI = 200                            # Image resolution

################################################################################
#                           END CONFIGURATION
################################################################################


def clean_solvent_name(name):
    """
    Clean solvent name for display purposes.
    Replaces "1-4-Dioxane" with "1,4-Dioxane" (case-insensitive)
    """
    if pd.isna(name):
        return name
    
    # Replace 2-4-7-9-Tetramethyl-4-7-decanediol with 2,4,7,9-Tetramethyl-4,7-decanediol (case-insensitive)
    cleaned = str(name).replace('1-4-Dioxane', '1,4-Dioxane')
    cleaned = cleaned.replace('1-4-Dioxane', '1,4-Dioxane')
    
    return cleaned


# Read title from row 0
print(f"Reading data from {INPUT_FILE}...")
df_title = pd.read_excel(INPUT_FILE, header=None, nrows=1)
diagram_title = str(df_title.iloc[0, 0]) if pd.notna(df_title.iloc[0, 0]) else 'Solvent Replacement Diagram'
print(f"Title: {diagram_title}")

# Read data
df = pd.read_excel(INPUT_FILE, header=None, skiprows=SKIP_ROWS)
df.columns = ['RED', 'Vol1', 'Solvent1', 'Vol2', 'Solvent2', 'RAG']
df['row_index'] = df.index

print(f"Loaded {len(df)} solvent pairs")

# Count connections per solvent
from collections import Counter
connections = Counter()
for _, row in df.iterrows():
    if pd.notna(row['Solvent1']) and pd.notna(row['Solvent2']):
        connections[row['Solvent1']] += 1
        connections[row['Solvent2']] += 1

# Select top N solvents
top_solvents = [s for s, _ in connections.most_common(N_SOLVENTS)]

# Filter to ALL pairs where both solvents are in top N (NO deduplication)
filtered_pairs = []

for _, row in df.iterrows():
    s1, s2 = row['Solvent1'], row['Solvent2']
    
    if pd.isna(s1) or pd.isna(s2):
        continue
    
    # Include ALL pairs where both solvents are in top N
    if s1 in top_solvents and s2 in top_solvents:
        filtered_pairs.append((s1, s2, row['Vol1'], row['Vol2'], row['RED'], row['row_index']))

print(f"✓ Selected {len(top_solvents)} solvents")
print(f"✓ Showing ALL {len(filtered_pairs)} pairs (no deduplication)")

# Create RED-to-letter mapping
unique_reds = sorted(set([red for _, _, _, _, red, _ in filtered_pairs if pd.notna(red)]))
letters = 'abcdefghijklmnopqrstuvwxyz'
red_to_letter = {red: letters[i] for i, red in enumerate(unique_reds)}

print(f"✓ {len(unique_reds)} unique RED values mapped to letters a-{red_to_letter[unique_reds[-1]]}")

# Function to wrap solvent names (break only at spaces, max 2 words per line)
def wrap_solvent_name(name):
    """Wrap solvent name if longer than 2 words"""
    # First clean the name for display
    name = clean_solvent_name(name)
    
    words = name.split()
    if len(words) <= 2:
        return name
    # Break into lines of 2 words each
    lines = []
    for i in range(0, len(words), 2):
        lines.append(' '.join(words[i:i+2]))
    return '\n'.join(lines)

# Setup circular positions
n = len(top_solvents)
angles = np.linspace(0, 2*np.pi, n, endpoint=False)
radius = 1.0
positions = {solvent: (radius * np.cos(angle), radius * np.sin(angle))
             for solvent, angle in zip(top_solvents, angles)}

# Assign colors
cmap = plt.colormaps.get_cmap('tab20')
solvent_colors = {solvent: cmap(i/n) for i, solvent in enumerate(top_solvents)}

# Create figure with space for legend
fig = plt.figure(figsize=(20, 16))
# Main plot area
ax = fig.add_axes([0.05, 0.05, 0.7, 0.9])
ax.set_facecolor('white')
ax.set_aspect('equal')
ax.axis('off')
ax.set_xlim(-1.7, 1.7)
ax.set_ylim(-1.7, 1.7)

fig.patch.set_facecolor('white')

# Draw tie-lines with gradient colors
for s1, s2, vol1, vol2, red_value, idx in filtered_pairs:
    x1, y1 = positions[s1]
    x2, y2 = positions[s2]
    
    # Calculate cross position based on volumes
    total_vol = vol1 + vol2
    t = vol1 / total_vol
    cross_x = x1 + t * (x2 - x1)
    cross_y = y1 + t * (y2 - y1)
    
    color1 = solvent_colors[s1]
    color2 = solvent_colors[s2]
    
    # Draw line segments
    ax.plot([x1, cross_x], [y1, cross_y], 
            color=color1, alpha=LINE_ALPHA, lw=LINE_WIDTH, zorder=1)
    ax.plot([cross_x, x2], [cross_y, y2],
            color=color2, alpha=LINE_ALPHA, lw=LINE_WIDTH, zorder=1)
    
    # Draw cross marker and letter
    if pd.notna(red_value):
        ax.plot(cross_x, cross_y, 'x', 
               color='black', 
               markersize=CROSS_SIZE,
               markeredgewidth=CROSS_WIDTH,
               alpha=0.9,
               zorder=3)
        
        # Add plain black letter (no background)
        letter = red_to_letter[red_value]
        ax.text(cross_x + LETTER_OFFSET_X, cross_y + LETTER_OFFSET_Y, letter,
               fontsize=LETTER_FONTSIZE,
               fontweight='bold',
               color='black',
               ha='center', va='center',
               zorder=4)

# Draw solvent points with redesigned labels
for solvent, (x, y) in positions.items():
    color = solvent_colors[solvent]
    
    # Draw point
    ax.scatter(x, y, s=POINT_SIZE, color=color, 
               edgecolors='black', linewidths=2.5, zorder=5, alpha=0.95)
    
    # Calculate text position
    angle = np.arctan2(y, x)
    
    # Position text outside the circle
    text_radius = 1.22
    text_x = text_radius * np.cos(angle)
    text_y = text_radius * np.sin(angle)
    
    # Determine horizontal alignment
    if text_x > 0.05:
        ha = 'left'
    elif text_x < -0.05:
        ha = 'right'
    else:
        ha = 'center'
    
    # Wrap the solvent name (this also cleans it)
    wrapped_name = wrap_solvent_name(solvent)
    
    # Label with colored border and white background
    ax.text(text_x, text_y, wrapped_name,
            fontsize=10, ha=ha, va='center',
            rotation=0,
            color='black',  # Black text
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', 
                     facecolor='white',  # White background
                     alpha=1.0,  # Fully opaque
                     edgecolor=color,  # Coloured border
                     linewidth=2.5),  # Thick border
            zorder=6)

# Title
ax.text(0, 1.58, diagram_title,
        ha='center', va='center', fontsize=15, fontweight='bold',
        color='black')
ax.text(0, 1.50, f'{len(top_solvents)} solvents | {len(filtered_pairs)} pairs',
        ha='center', va='center', fontsize=11.5, color='#333333', style='italic')

# Create compact legend panel 
ax_legend = fig.add_axes([0.78, 0.35, 0.18, 0.5])
ax_legend.axis('off')

# Legend title
ax_legend.text(0.5, 0.98, 'RED Values',
              ha='center', va='top',
              fontsize=13, fontweight='bold',
              transform=ax_legend.transAxes)

# Add separator line
ax_legend.plot([0.1, 0.9], [0.94, 0.94], 'k-', lw=1, 
              transform=ax_legend.transAxes)

# Add each letter-RED pair with compact spacing
y_position = 0.90

for red_val in unique_reds:
    letter = red_to_letter[red_val]
    
    # Plain black letter 
    ax_legend.text(0.15, y_position, letter,
                  fontsize=LEGEND_FONTSIZE,
                  fontweight='bold',
                  color='black',
                  ha='center', va='center',
                  transform=ax_legend.transAxes)
    
    # RED value 
    ax_legend.text(0.30, y_position,
                  f'{float(red_val):.3f}',
                  va='center', ha='left',
                  fontsize=LEGEND_FONTSIZE,
                  transform=ax_legend.transAxes)
    
    y_position -= LEGEND_LINE_SPACING

# Add general legend info
y_position -= 0.04
ax_legend.plot([0.1, 0.9], [y_position, y_position], 'k-', lw=0.5,
              alpha=0.3, transform=ax_legend.transAxes)

y_position -= 0.04
legend_info_text = "× = mixture point\nLine colour matches solvent"

ax_legend.text(0.5, y_position,
              legend_info_text,
              ha='center', va='top',
              fontsize=8,
              transform=ax_legend.transAxes)

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved {OUTPUT_FILE}")

