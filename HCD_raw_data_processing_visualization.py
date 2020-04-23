

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sps
#==
# HCD data
#==
# Access file location
hcd_sheets_raw = pd.read_excel('C:/Users/Tashi Wischmeyer/Documents/HCD_data_management/Carbon_dioxide_unattended_CO2.xlsx', index_col = None, sheet_name = 'Raw_data')

target_contaminant = input('What contaminant is the target for this test?')

hcd = pd.DataFrame(hcd_sheets_raw)
r_count, c_count = hcd.shape
date_type = type(hcd.iloc[1,0])

# determine input concentration from col name
if target_contaminant in hcd.columns:
    if(target_contaminant == 'CO2'):
        input_concentration = [20, 10, 5, 4, 3, 2, 1.8, 1.4, 1.2, 1, 0]
else:
    print('Target contaminant not found')

# find index column for target contaminant
col = 0
for name in hcd.columns:
    if name == target_contaminant:
        break
    col = col + 1

concentration_CO2 = [[0],[0],[0],[0],[0],[0],[0],[0],[0]]
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
    theo_diff = np.abs((input_concentration[i]-input_concentration[(i+1)%len(input_concentration)])*perc)

    # look for concentration change
    # dependent on expected concentrations
    if exp_diff > theo_diff:
        calc = np.round(np.mean(hcd.iloc[(row-30):(row-14), col]), decimals = 4, out = None)
        val = concentration_CO2[j]
        val.append(calc)
        if(i == 10):
            j = j + 1
        i = (i + 1)%len(input_concentration)
        row = row + 3*span
    else:
        row = row + 1
# if there are not the correct number of entries, the last non-NULL row value is used to calculate the final entry
if(len(concentration_CO2[col-1]) < len(input_concentration)):
    calc = np.round(np.mean(hcd.iloc[(row-30):(row-15), col]), decimals = 4, out = None)
    val = concentration_CO2[col-1]
    val.append(calc)

for k in range(j+1):
    if concentration_CO2[k-1][0] == 0:
        del concentration_CO2[k-1][0]
    concentration_CO2[k-1].sort(reverse = True)
    print(concentration_CO2[k-1])

# average
avg_response_concentration = np.mean(concentration_CO2, axis = 0)
# std dev
std_response_concentration = np.std(concentration_CO2, axis = 0)
# var
var_response_concentration = np.var(concentration_CO2, axis = 0)
# std error
se_response_concentration = sps.sem(concentration_CO2, axis = 0)

# plotting
area = np.pi*3
figCO2calib = plt.scatter(input_concentration, avg_response_concentration, s = area, alpha = 0.5, c = 'k')
plt.ylim(ymax = 20, ymin = 0)
plt.title('FTIR Calibration curve for CO2')
plt.xlabel('Input concentration [ppm]')
plt.ylabel('Response concentration [ppm]')
plt.errorbar(input_concentration, avg_response_concentration, yerr = se_response_concentration)
plt.show()
# save figure
plt.savefig('C:/Users/Tashi Wischmeyer/Documents/HCD_data_management/Test_curve.pdf')
