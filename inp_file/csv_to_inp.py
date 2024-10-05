import pandas as pd

# Define the paths
csv_path = r'C:\Users\meian\Desktop\Comp Eng Project\Project_Py\pythonProject1\Simulation\Tmp298.15K_StrRt0.1.csv'
output_inp_path = r'C:\Users\meian\Desktop\Comp Eng Project\Project_Py\pythonProject1\Simulation\Material_DP1000_TxK_SRy.inp'

# Read the CSV file and strip leading/trailing spaces from column names
data = pd.read_csv(csv_path)
data.columns = data.columns.str.strip()

# Print the column names to verify
print(data.columns)

# Extract the necessary columns
strain = data['strain']
rd = data['RD']
dd = data['DD']
td = data['TD']
biaxial = data['biaxial']
r_0 = data['r_RD']
r_45 = data['r_DD']
r_90 = data['r_TD']

# Create the header for the .inp file
header = f"""** MATERIALS
** 
*Material, name=DP1000
*Depvar
     2,
*User Material, constants={len(data.columns) * len(data)}
** Young's Modulus, Poisson Ratio
 210000,        0.3,          0.0,         0.0,            0.0,            0.0,         0.0,   0.
** Flow curves along 0, 45, 90 and biaxial; r-values along 0, 45 and 90
"""


# Create the formatted data string
formatted_data = ''
for i in range(len(strain)):
    formatted_data += f'{strain[i]},      {rd[i]},      {dd[i]},      {td[i]},      {biaxial[i]},      {r_0[i]},      {r_45[i]},      {r_90[i]}\n'

# Add the density at the end
density = "*Density\n 7.85e-09,\n"

# Combine the header, formatted data, and density
inp_content = header + formatted_data + density

# Write the new .inp file
with open(output_inp_path, 'w') as file:
    file.write(inp_content)

print(f'Generated .inp file at {output_inp_path}')