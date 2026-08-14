import csv
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import pandas as pd


# Determine the path of the current file
if getattr(sys, 'frozen', False):
    current_directory = os.path.dirname(sys.executable)
else:
    current_directory = os.path.dirname(os.path.abspath(__file__))


# Paths
Data_LM_parameter = os.path.join(current_directory, 'Parameter_Set-Up_Pharma.xlsx')
Sdata = pd.read_excel(Data_LM_parameter, sheet_name='Solvent', header=None)
Data_Parameter = os.path.join(current_directory, 'Parameter_Set-Up_Pharma.xlsx')
data_df = pd.read_excel(Data_Parameter, sheet_name='Polymer', header=None)
csv_filename = os.path.join(current_directory, 'results/All RED-values.csv')
result_csv_data_name = os.path.join(current_directory, 'results/Solvent list.csv')
output_directory = os.path.join(current_directory, 'results')
Data_Parameter = os.path.join(current_directory, 'Parameter_Set-Up_Pharma.xlsx')
if not os.path.exists(output_directory):
    try:
        os.makedirs(output_directory)
    except OSError as e:
        print(f"Error while creating the folder {output_directory}: {e}")
    else:
        print(f"The folder {output_directory} was successfully created.")
else:
    print(f"The folder {output_directory} already exists.")


# ---------------------------------------------------------------------------
# RAG helpers
# ---------------------------------------------------------------------------
# Priority order used to pick the "worst" score in a solvent mix.
# Lower number = worse 
RAG_PRIORITY = {'Grey': 0, 'R': 1, 'A': 2, 'G': 3}

def worst_rag(rag1, rag2):
    """Return whichever of two RAG scores is worse (lower priority number).
    If either value is missing/unknown, return the other one.
    """
    if rag1 not in RAG_PRIORITY:
        return rag2
    if rag2 not in RAG_PRIORITY:
        return rag1
    return rag1 if RAG_PRIORITY[rag1] <= RAG_PRIORITY[rag2] else rag2


# Creating csv.-files
def writeData(filename, header, data):
    with open(os.path.join(current_directory, filename), mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')
        writer.writerow(header)
        for row in data:
            writer.writerow(row)
    print("Data has been added to the CSV file:", filename)

# header
header1 = ["Vol.-% of S1", "Solvent 1", "Vol.-% of S2", "Solvent 2", "RED Polymer 1", "RED Polymer 2", "RED Polymer 3",
           "RED Polymer 4", "RED Polymer 5"]
header2 = ["Vol.-% of S1", "Solvent 1", "Vol.-% of S2", "Solvent 2", "dD (S1)", "dP (S1)", "dH (S1)", "dD (S2)",
           "dP (S2)", "dH (S2)", "RED"]
header3 = ["Vol.-% of S1", "Solvent 1", "Vol.-% of S2", "Solvent 2", "RED Polymer 1", "RED Polymer 2", "RED Polymer 3",
           "RED Polymer 4", "RED Polymer 5", "Product/Sum"]

# Defining of various functions
def calculate_D_x(S1_x, dD_S1, S2_x, dD_S2):
    dD_x = ((S1_x * dD_S1 + S2_x * dD_S2) / 100)
    return dD_x
def calculate_P_x(S1_x, dP_S1, S2_x, dP_S2):
    dP_x = ((S1_x * dP_S1 + S2_x * dP_S2) / 100)
    return dP_x
def calculate_H_x(S1_x, dH_S1, S2_x, dH_S2):
    dH_x = ((S1_x * dH_S1 + S2_x * dH_S2) / 100)
    return dH_x
def calculate_RED(dD_x, dD_B, dP_x, dP_B, dH_x, dH_B, R_0):
    RED = math.sqrt((4* ( dD_x - dD_B)**2 + (dP_x - dP_B)**2 + (dH_x - dH_B)**2)) /R_0
    return RED
# Function of RED-value calculation and data saving
def calculate_RED_Function(All_RED_values_BM1,dD_B, dP_B, dH_B, R_0,kk):
    MB_list=[]
    pure_solvent_values = []
    # Calculation of solvent mixtures
    for S1_x in start_value_combinations:
        S2_x = 100-S1_x
        for ii in range(len(allvalues)):
            solvent_name, dD_S1_value, dP_S1_value, dH_S1_value = allvalues[ii]
            for jj in range(ii):
                solvent_name2, dD_S2_value, dP_S2_value, dH_S2_value = allvalues[jj]
                if solvent_name != solvent_name2:
                    dD_x_result = calculate_D_x(S1_x, dD_S1_value, S2_x, dD_S2_value)
                    dH_x_result = calculate_H_x(S1_x, dH_S1_value, S2_x, dH_S2_value)
                    dP_x_result = calculate_P_x(S1_x, dP_S1_value, S2_x, dP_S2_value)
                    RED = calculate_RED(dD_x_result, dD_B, dP_x_result, dP_B, dH_x_result, dH_B, R_0)
                    RED = round(RED, 3)
                #Saving
                MB_list.append([S1_x, solvent_name, S2_x, solvent_name2, dD_S1_value, dP_S1_value, dH_S1_value, dD_S2_value, dP_S2_value, dH_S2_value, RED])
                All_RED_values_BM1.append([S1_x, solvent_name, S2_x, solvent_name2, RED])

    # Calculation of pure solvents
    for pure_solvent in allvalues:
        solvent_name, dD_S1_value, dP_S1_value, dH_S1_value = pure_solvent
        dD_S2_value, dP_S2_value, dH_S2_value = 0, 0, 0
        dD_x_result = calculate_D_x(100, dD_S1_value, 0, dD_S2_value)
        dH_x_result = calculate_H_x(100, dH_S1_value, 0, dH_S2_value)
        dP_x_result = calculate_P_x(100, dP_S1_value, 0, dP_S2_value)
        RED = calculate_RED(dD_x_result, dD_B, dP_x_result, dP_B, dH_x_result, dH_B, R_0)
        pure_solvent_values.append([solvent_name, dD_x_result, dP_x_result, dH_x_result, RED])
        RED = round(RED, 3)
        # Saving
        MB_list.append([100, solvent_name, 0, str('-'), dD_S1_value, dP_S1_value, dH_S1_value, dD_S2_value, dP_S2_value, dH_S2_value, RED])
        All_RED_values_BM1.append([100, solvent_name, 0, str('-'), RED])

    print(str(kk) +" Polymer has been calculated")



###############################################  Start of the program  #################################################
# Lists
allvalues = []
List_allRED =[]
All_RED_values_BM1 = []
Ranking = []


# dicts
dict_par = {}
red_dict = {}


# ---------------------------------------------------------------------------
# Import of all Parameters from the Excel file
# ---------------------------------------------------------------------------

# Read in solvent data AND the RAG column (column index 4)
rag_dict = {}   # maps solvent name -> RAG score string

for index, row in Sdata.iterrows():
    solvent_name = row[0]
    # Check whether the row contains data
    if index > 0:
        dD_S1_value = float(row[1])
        dP_S1_value = float(row[2])
        dH_S1_value = float(row[3])
        # Create tuple and add to allvalues
        data_point = (solvent_name, dD_S1_value, dP_S1_value, dH_S1_value)
        allvalues.append(data_point)

        # --- RAG: read from column 4 (E in Excel) ---
        rag_value = row[4] if len(row) > 4 else None
        if pd.notna(rag_value):
            rag_dict[solvent_name] = str(rag_value).strip()
        else:
            rag_dict[solvent_name] = ''   # empty if not provided

# Read the NumerOfPoly value from the Excel file
numer_of_poly_df = pd.read_excel(Data_Parameter, sheet_name='Polymer', header=None, usecols="B", nrows=1)
NumerOfPoly = numer_of_poly_df.iat[0, 0]
# Read the value for STEP_SIZE from the Excel file
step_size_df = pd.read_excel(Data_Parameter, sheet_name='Polymer', header=None, usecols="B", skiprows=1, nrows=1)
STEP_SIZE = step_size_df.iat[0, 0]
# Read the Excel table into a DataFrame
data_df = pd.read_excel(Data_Parameter, sheet_name='Polymer', header=None)
# Read the polymer parameters (dD_B, dP_B, dH_B and R_0) into dict_par from the Excel file
for jj in range(1, int(NumerOfPoly) + 1):
    dict_par[f"dD_B_{jj}"] = data_df.iat[4, jj]
    dict_par[f"dP_B_{jj}"] = data_df.iat[5, jj]
    dict_par[f"dH_B_{jj}"] = data_df.iat[6, jj]
    dict_par[f"R_0_{jj}"] = data_df.iat[7, jj]


# Start value combinations of S1_x and S2_x
start_value_combinations = np.arange(STEP_SIZE, 100, STEP_SIZE)


# Calculation of all RED values. A loop is used for this, which runs through all parameters within the dict and saves the values.
for kk in range(1,NumerOfPoly+1):
    calculate_RED_Function(All_RED_values_BM1,dict_par["dD_B_"+str(kk)], dict_par["dP_B_"+str(kk)], dict_par["dH_B_"+str(kk)], dict_par["R_0_"+str(kk)],kk)


# Calculation number of elements of each polymer
size_single = int(len(All_RED_values_BM1)/NumerOfPoly)


# Save RED values in dict
for pp in range(NumerOfPoly+1):
    red_dict["red_values_{0}".format(pp+1)] = [item[-1] for item in All_RED_values_BM1[size_single*(pp):size_single*(pp+1)]]


# Extract text formats from lists
AmountS1 = [item[0] for item in All_RED_values_BM1[:size_single]]
Solvent1 = [item[1] for item in All_RED_values_BM1[:size_single]]
AmountS2 = [item[2] for item in All_RED_values_BM1[:size_single]]
Solvent2 = [item[3] for item in All_RED_values_BM1[:size_single]]


# Array based on red.dict
Array_RED_values = np.zeros((size_single,NumerOfPoly))


# Write and save data
for tt in range(size_single):
    # loop through polymers
    for hh in range(0, int(NumerOfPoly)):
        Array_RED_values[tt][hh] = (red_dict["red_values_" + str(hh + 1)][tt])

for i in range(size_single):
    expanded_subarray = np.hstack((AmountS1[i], Solvent1[i], AmountS2[i], Solvent2[i], Array_RED_values[i][:]))
    List_allRED.append(expanded_subarray)


# Create All RED-values csv
# Iteration over the lines of List_allRED
for row in List_allRED:
    red_new_values = [float(value) for value in row[4:]]

    # Checking the number of polymers
    if NumerOfPoly == 1:  # If there is only one polymer
        # Check whether all RED values are zero
        if all(value == 0 for value in red_new_values):
            sort = 0
        else:
            # If not all RED values are zero, use the RED values as a product
            sort = sum(value for value in red_new_values)
    else:  # If there is more than one polymer
        # Calculation of the product with zeros as 0.001
        sort = math.prod([value if value != 0 else 0.001 for value in red_new_values]) / sum(red_new_values)

    # Add the calculated values to the line
    row = row.tolist()
    while len(row) < 9:
        row.append('')
    row.append(sort)
    Ranking.append(row)

# Conversion of strings to floats (ignores empty strings)
def safe_float(value):
    try:
        return float(value)
    except ValueError:
        return None  # If conversion fails, None is returned

# Division into two groups with secure conversion
group1 = [row for row in Ranking
    if all(value < 1 for value in [safe_float(v) for v in row[4:-1]] if value is not None)]
group2 = [row for row in Ranking if any(safe_float(value) is not None and safe_float(value) >= 1 for value in row[4:-1])]

# Sort by the last element (sort value)
sorted_group1 = sorted(group1, key=lambda x: x[-1])
sorted_group2 = sorted(group2, key=lambda x: x[-1])

sorted_rowsI = sorted_group1 + sorted_group2

# Function for writing to CSV
def writeData(filename, header, data):
    import csv
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')
        writer.writerow(header)
        for row in data:
            row_without_sort = row[:-1]
            writer.writerow(row_without_sort)

header = ["Vol.-% of S1", "Solvent 1", "Vol.-% of S2", "Solvent 2", "RED Polymer 1", "RED Polymer 2",
          "RED Polymer 3", "RED Polymer 4", "RED Polymer 5"]

writeData(csv_filename, header, sorted_rowsI)


# ---------------------------------------------------------------------------
# Create file "Solvent list.csv"  — now includes a RAG Score column
# ---------------------------------------------------------------------------

# Helper: compute the combined RAG score for a solvent-mix row.
# For a pure solvent (S2 == '-') we simply use S1's RAG.
# For a mixture we take the *worst* of S1 and S2 (R > A > Grey > G).
def get_mix_rag(solvent1_name, solvent2_name):
    rag1 = rag_dict.get(solvent1_name, '')
    rag2 = rag_dict.get(solvent2_name, '') if solvent2_name != '-' else ''
    if rag1 == '' and rag2 == '':
        return ''           # both missing – leave blank
    if rag1 == '':
        return rag2
    if rag2 == '':
        return rag1
    return worst_rag(rag1, rag2)


with open(result_csv_data_name, mode='w', newline='') as result_csv_file:
    csv_writer = csv.writer(result_csv_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    # Write the header — RAG Score is inserted right after Solvent 2 (position 5)
    header1_with_rag = ["Ranking", "Vol.-% of S1", "Solvent 1", "Vol.-% of S2", "Solvent 2", "RAG Score",
                        "RED Polymer 1", "RED Polymer 2", "RED Polymer 3", "RED Polymer 4", "RED Polymer 5"]
    csv_writer.writerow(header1_with_rag)

    lines = []

    # Iteration over the lines of List_allRED
    for row in List_allRED:
        red_new_values = [float(value) for value in row[4:]]

        # Check whether all RED values are less than 1
        if all(value < 1 for value in red_new_values):
            # Checking the number of polymers
            if NumerOfPoly == 1:   # If there is only one polymer
                # Check whether all RED values are zero
                if all(value == 0 for value in red_new_values):
                    sort = 0
                else:
                    # If not all RED values are zero, use the RED values as a product
                    sort = sum(value for value in red_new_values)

            else:  # If there is more than one polymer
                # Calculation of the product with zeros as 0.001
                sort = math.prod([value if value != 0 else 0.001 for value in red_new_values]) / sum(red_new_values)

            # Add the calculated values to the line
            row = row.tolist()
            while len(row) < 9:
                row.append('')
            row.append(sort)
            lines.append(row)

    # Sort
    sorted_rows = sorted(lines, key=lambda x: x[-1])

    # Writing in CSV-file (with RAG column inserted)
    i = 1
    for row in sorted_rows:
        # row layout: [AmountS1, Solvent1, AmountS2, Solvent2, RED1..RED5, sort]
        solvent1_name = row[1]
        solvent2_name = row[3]
        mix_rag = get_mix_rag(solvent1_name, solvent2_name)

        # Build output row: Ranking, AmountS1, Solvent1, AmountS2, Solvent2, RAG, RED1..., (drop sort)
        out_row = [i, row[0], row[1], row[2], row[3], mix_rag] + row[4:-1]
        csv_writer.writerow(out_row)
        i += 1


header = ["Ranking", "Vol.-% of S1", "Solvent 1", "Vol.-% of S2", "Solvent 2", "RAG Score",
          "RED Polymer 1", "RED Polymer 2", "RED Polymer 3", "RED Polymer 4", "RED Polymer 5", "sort"]

# Sort data
sorted_rows = sorted(lines, key=lambda x: (
    float(x[5]) if len(x) >= 6 and x[5] and all(value in ('', '0') for value in x[6:9]) else x[-1]),)


################################### Determining the best mixture to display it in the diagram #####################
with open(result_csv_data_name, mode='r') as result_csv_file:
    csv_reader = csv.reader(result_csv_file, delimiter=';')
    next(csv_reader)  # Skip the header
    # Try to read the first line after the header
    try:
        best_solvent_mixture = next(csv_reader)
    except StopIteration:
        # If there are no lines after the header, best_solvent_mixture is not set
        best_solvent_mixture = None

# If there is a valid line, continue processing
# NOTE: column indices shifted +1 because RAG Score was inserted at position 5
if best_solvent_mixture:
    best_S1_x = float(best_solvent_mixture[1])
    best_S1 = best_solvent_mixture[2]
    best_S2_x = float(best_solvent_mixture[3])
    best_S2 = best_solvent_mixture[4]
    # best_solvent_mixture[5] is now RAG Score — not used for plotting


    # If the first solvent content is 100%, use its Hansen solubility parameter directly
    if best_S1_x == 1:
        best_dD_x, best_dP_x, best_dH_x = dD_S1_value, dP_S1_value, dH_S1_value
    else:
        # Search for Hansen parameters for the solvents, if available
        dD_S1, dP_S1, dH_S1 = next(((dD, dP, dH) for name, dD, dP, dH in allvalues if name == best_S1), (None, None, None))
        dD_S2, dP_S2, dH_S2 = next(((dD, dP, dH) for name, dD, dP, dH in allvalues if name == best_S2), (None, None, None))

        # If both solvents are valid, calculate mean Hansen parameters
        if dD_S1 is not None and dD_S2 is not None:
            best_dD_x = calculate_D_x(best_S1_x, dD_S1, best_S2_x, dD_S2)
            best_dP_x = calculate_P_x(best_S1_x, dP_S1, best_S2_x, dP_S2)
            best_dH_x = calculate_H_x(best_S1_x, dH_S1, best_S2_x, dH_S2)
        else:
            # If only one of the two solvents is valid, adopt its values
            best_dD_x, best_dP_x, best_dH_x = (dD_S1, dP_S1, dH_S1) if dD_S1 is not None else (dD_S2, dP_S2, dH_S2)

##################################################  Diagram  #########################################################

# Define Function
def plot_sphere(R_0, dD_B, dP_B, dH_B, figNmb):
    fig = plt.figure(figsize=(8, 24))

    # View 1
    ax1 = fig.add_subplot(311, projection='3d')  # 3 Subplots, first position
    ax1.view_init(azim=20, elev=20)
    plot_single_view(ax1, R_0, dD_B, dP_B, dH_B, figNmb)
    ax1.set_title('View 1')

    # View 2
    ax2 = fig.add_subplot(312, projection='3d')  # 3 Subplots, second position
    ax2.view_init(azim=80, elev=20)
    plot_single_view(ax2, R_0, dD_B, dP_B, dH_B, figNmb)
    ax2.set_title('View 2')

    # View 3
    ax3 = fig.add_subplot(313, projection='3d')  # 3 Subplots, third position
    ax3.view_init(azim=140, elev=20)
    plot_single_view(ax3, R_0, dD_B, dP_B, dH_B, figNmb)
    ax3.set_title('View 3')

    # Save the combined figure
    screenshot_filename = os.path.join(current_directory, 'results/HSP-Plot.png')
    canvas = FigureCanvas(fig)
    canvas.print_png(screenshot_filename)
    print(f'HSP-Plot has been saved as {screenshot_filename}')

def plot_single_view(ax, R_0, dD_B, dP_B, dH_B, figNmb):
    for jj in range(1, int(NumerOfPoly) + 1):
        R_0 = dict_par["R_0_{0}".format(jj)]
        dD_B = dict_par["dD_B_{0}".format(jj)]
        dP_B = dict_par["dP_B_{0}".format(jj)]
        dH_B = dict_par["dH_B_{0}".format(jj)]
        # draw sphere
        u, v = np.mgrid[0:2 * np.pi:50j, 0:np.pi:50j]
        x = np.cos(u) * np.sin(v) * R_0 + dD_B
        y = np.sin(u) * np.sin(v) * R_0 + dP_B
        z = np.cos(v) * R_0 + dH_B
        color = plt.get_cmap('jet')(jj / int(NumerOfPoly))
        ax.plot_surface(x, y, z, color=color, alpha=0.15, label=f'Polymer {jj}')

    ax.set_xlabel('dD')
    ax.set_ylabel('dP')
    ax.set_zlabel('dH')
    ax.set_xlim((0, 30))
    ax.set_ylim((0, 30))
    ax.set_zlim((0, 30))
    ax.set_box_aspect([1, 1, 1])
    ax.legend()

    # Only plot if the solvent list contains data
    if best_solvent_mixture:
        try:
            best_S1_x = float(best_solvent_mixture[1])
            best_S1 = best_solvent_mixture[2]
            best_S2_x = float(best_solvent_mixture[3])
            best_S2 = best_solvent_mixture[4]

            if best_dD_x is not None and best_dP_x is not None and best_dH_x is not None:
                label_parts = []
                if best_S1_x > 0:
                    label_parts.append(f"{int(best_S1_x)}% {best_S1}")
                if best_S2_x > 0:
                    label_parts.append(f"{int(best_S2_x)}% {best_S2}")

                label_text = " + ".join(label_parts)

                ax.scatter(best_dD_x, best_dP_x, best_dH_x, color='black', s=75, marker='x', label=label_text)
                ax.legend()

        except (ValueError, IndexError) as e:
            print(f"Error when processing CSV data for plotting: {e}")

# Value import and plot call:
jj = 1
R_0 = dict_par["R_0_" + str(jj)]
dD_B = dict_par["dD_B_" + str(jj)]
dP_B = dict_par["dP_B_" + str(jj)]
dH_B = dict_par["dH_B_" + str(jj)]
plot_sphere(R_0, dD_B, dP_B, dH_B, jj)

print("Current working register:", os.getcwd())
