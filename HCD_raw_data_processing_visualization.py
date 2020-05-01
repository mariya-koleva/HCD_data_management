

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
        input_concentration = [0, 20, 10, 5, 4, 3, 2, 1.8, 1.4, 1.2, 1]
else:
    print('Target contaminant not found')

# find index column for target contaminant
col = 0
for name in hcd.columns:
    if name == target_contaminant:
        break
    col = col + 1

concentration = []
cycle_vals = []
val_index = []
cycle_index = []
span = 5
perc = 0.5
row = 30
i = 0
while(row < (r_count-span)):
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
        # if we have negative values, zero them because it doesn't make sense to have negatives
        if(calc < 0):
            calc = 0
        cycle_vals.append(calc)
        if(i == 10):
            cycle_vals.sort(reverse = True)
            print(cycle_vals)
            concentration.append(cycle_vals)
            cycle_vals = []
        i = (i + 1)%len(input_concentration)
        row = row + 6*span
    else:
        row = row + 1
# if there are not the correct number of entries, the last non-NULL row value is used to calculate the final entry
if(len(concentration[col-1]) < len(input_concentration)):
    calc = np.round(np.mean(hcd.iloc[(row-30):(row-15), col]), decimals = 4, out = None)
    # if we have negative values, zero them because it doesn't make sense to have negatives
    if(calc < 0):
        calc = 0
    val = concentration[col-1]
    val.append(calc)

input_concentration.sort(reverse = True)
# average
avg_response_concentration = np.mean(concentration, axis = 0)
# std dev
std_response_concentration = np.std(concentration, axis = 0)
# var
var_response_concentration = np.var(concentration, axis = 0)
# std error
se_response_concentration = sps.sem(concentration, axis = 0)

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
