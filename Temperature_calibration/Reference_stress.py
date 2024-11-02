import pandas as pd

# Calculate relative temperature values based on a reference temperature column



file_path = r'Temperature_calibration\Fitted_FC_TmpDpn_RD.xlsx'
sheet_name = 'StrRt0.0001' 
output_file_path = f'Temperature_calibration\Exp_FC_TmpDpn_RD_{sheet_name}_ref.csv'
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Column name for the reference temperature
reference_column = 'Tmp298.15K'

# Check if the reference column exists
if reference_column not in df.columns:
    raise ValueError(f"Reference column '{reference_column}' not found in the CSV file.")

# Iterate over columns to create new columns by dividing with the reference column
for column in df.columns:
    # Skip the reference and non-temperature columns (e.g., "Strain")
    if column.startswith("Tmp"):
        new_column_name = f"{column}_ref"
        df[new_column_name] = df[column] / df[reference_column]

# Save the modified DataFrame to a new CSV file
df.to_csv(output_file_path, index=False)
print(f"Modified file saved as '{output_file_path}'")
