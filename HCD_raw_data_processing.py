

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#==
# HCD data
#==
# Access file location
hcd_sheets_raw = pd.read_excel('C:/Users/Tashi Wischmeyer/Documents/HCD_data/Carbon_dioxide_unattended_CO2.xlsx', index_col = None, sheet_name = 'Raw_data')

# hcd_sheets_raw = []
# for sheet in hcd_xslx_raw.sheet_names:
#     hcd_sheets_raw.append(hcd_xslx_raw.parse(sheet))
#     hcd_raw = pd.concat(hcd_sheets_raw)

hcd = pd.DataFrame(hcd_sheets_raw)
r_count, c_count = hcd.shape
date_type = type(hcd.iloc[1,0])

col = 0;
for name in hcd.columns:
    if(name == 'CO2'):
        break
    col = col+1

in_conc = [0, 20, 10, 5, 4, 3, 2, 1.8, 1.4, 1.2, 1]

conc = [[0],[0],[0],[0],[0],[0],[0],[0],[0]]
span = 5
perc = 0.55
j = 0
# for col in range(1, 10):
row = 30
i = 0
while(row < (r_count-span)) and (j<9):
    # checks that values are non-NULL
    if(type(hcd.iloc[row, 0]) != date_type):
        break
    # experimental difference of 5 pt averages
    exp_diff = np.abs(np.mean(hcd.iloc[(row-span):row, col]) - np.mean(hcd.iloc[(row+span):(row+2*span), col]))
    # threshold for theoretical difference between concentrations
    theo_diff = np.abs((in_conc[i]-in_conc[(i+1)%len(in_conc)])*perc)

    # look for concentration change
    # dependent on expected concentrations
    if exp_diff > theo_diff:
        calc = np.round(np.mean(hcd.iloc[(row-30):(row-14), col]), decimals = 4, out = None)
        val = conc[j]
        val.append(calc)
        if(i == 10):
            j = j + 1
        i = (i + 1)%len(in_conc)
        row = row + 3*span
    else:
        row = row + 1
# if there are not the correct number of entries, the last non-NULL row value is used to calculate the final entry
if(len(conc[col-1]) < len(in_conc)):
    calc = np.round(np.mean(hcd.iloc[(row-30):(row-15), col]), decimals = 4, out = None)
    val = conc[col-1]
    val.append(calc)

for k in range(j+1):
    del conc[k-1][0]
    print(conc[k-1])
raw_data_process.py
