import pandas as pd

# Define the paths
params_path = r'subroutine/TmpDpn/9_Parameters_Fitted_FC_TmpDpn_RD_StrRt0.0001_final.csv'
flow_curves_r_path = r'subroutine/TmpDpn/Tmp298.15K_StrRt0.0001.csv'
output_inp_path = r'subroutine/TmpDpn/Material_DP1000_TxK_SRy.inp'

# Read the CSV file and strip leading/trailing spaces from column names
fc_r_data = pd.read_csv(flow_curves_r_path)
params_data = pd.read_csv(params_path)
fc_r_data.columns = fc_r_data.columns.str.strip()
params_data.columns = params_data.columns.str.strip()

# Print the column names to verify
print(fc_r_data.columns)
print(params_data.columns)

# Extract the necessary columns
strain = fc_r_data['strain']
rd = fc_r_data['RD']
dd = fc_r_data['DD']
td = fc_r_data['TD']
biaxial = fc_r_data['biaxial']
r_0 = fc_r_data['r_RD']
r_45 = fc_r_data['r_DD']
r_90 = fc_r_data['r_TD']

c1 = params_data['C1']
c2 = params_data['C2']
c3 = params_data['C3']
c4 = params_data['C4']
c5 = params_data['C5']
c6 = params_data['C6']
c7 = params_data['C7']
c8 = params_data['C8']
c9 = params_data['C9']

# Create the header for the .inp file
header = f"""** MATERIALS
** 
*Material, name=DP1000
*Depvar
     2,
*User Material, constants=6749
** Young's Modulus, Poisson Ratio
 210000,        0.3,          0.0,         0.0,            0.0,            0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,         0.0,   0.
** Flow curves along 0, 45, 90 and biaxial; r-values along 0, 45 and 90; Temperature Softening and Dynamic Strain Aging parameters
"""


# Create the formatted data string
formatted_data = ''
for i in range(len(strain)):
    formatted_data += f'{strain[i]},      {rd[i]},      {dd[i]},      {td[i]},      {biaxial[i]},      {r_0[i]},      {r_45[i]},      {r_90[i]},      {c1[i]},      {c2[i]},      {c3[i]},      {c4[i]},      {c5[i]},      {c6[i]},      {c7[i]},      {c8[i]},      {c9[i]}\n'

# Add the density at the end
density = "*Density\n 7.85e-09,\n"

# Combine the header, formatted data, and density
inp_content = header + formatted_data + density

# Write the new .inp file
with open(output_inp_path, 'w') as file:
    file.write(inp_content)

print(f'Generated .inp file at {output_inp_path}')