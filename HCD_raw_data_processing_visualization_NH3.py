

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sps
import seaborn
#==
# HCD data
#==
# Access file location
hcd_sheets_raw = pd.read_excel('C:/Users/MKOLEVA/Documents/Masha/Projects/HCD/Experiments/FTIR_data/Multicomponent results/Ammonia_unattended_27.xlsx', index_col = None, sheet_name = 'Raw_data')

while True:
    try:
        target_contaminant = input('What contaminant is the target for this test?')
    except ValueError:
        print("Oopsie, I don't have the capability to deal with that contaminant.")
        continue
    else:
        break
if target_contaminant == "NH3":
    print("This contaminant is in the databse and I'm proceeding with calculations.") 
else:
    print("Oopsie, not able to proceed further as it seems I cannot handle this contaminant.")

hcd = pd.DataFrame(hcd_sheets_raw)
r_count, c_count = hcd.shape
date_type = type(hcd.iloc[1,0])

# determine input concentration from col name
if target_contaminant in hcd.columns:
    if(target_contaminant == 'NH3'):
        input_concentration = [0, 1, 0.7, 0.3, 0.1, 0.08, 0.07, 0.05, 0.03]
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
span = 2
perc = 0.4
row = 10
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
        calc = np.round(np.mean(hcd.iloc[(row-10):(row-4), col]), decimals = 4, out = None)
        # if we have negative values, zero them because it doesn't make sense to have negatives
#        if(calc < 0):
#            calc = 0
        cycle_vals.append(calc)
        if(i == 8):
            cycle_vals.sort(reverse = True)
            print(cycle_vals)
            concentration.append(cycle_vals)
            cycle_vals = []
        i = (i + 1)%len(input_concentration)
        row = row + 6*span
    else:
        row = row + 1
# if there are not the correct number of entries, the last non-NULL row value is used to calculate the final entry
        # the number of concentration points may be different for different contaminant
if(len(concentration[col-1]) < len(input_concentration)):
    calc = np.round(np.mean(hcd.iloc[(row-10):(row-5), col]), decimals = 4, out = None)
    # if we have negative values, zero them because it doesn't make sense to have negatives
#    if(calc < 0):
#        calc = 0
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
# ymax to be generic
plt.ylim(ymax = 1, ymin = -0.5)
plt.title('FTIR Calibration curve for NH3')
plt.xlabel('Input concentration [ppm]')
plt.ylabel('Response concentration [ppm]')
plt.errorbar(input_concentration, avg_response_concentration, yerr = se_response_concentration)
plt.show()

#figCO2calib2 = plt.scatter(input_concentration,concentration)
#
#input_concentration = pd.DataFrame(input_concentration)
#seaborn.pairplot(concentration, vars['input concentration','response concentration'],kind='reg')
# save figure
plt.savefig('C:/Users/MKOLEVA/Documents/Masha/Projects/HCD/Experiments/FTIR_data/Multicomponent results/Test_curve_NH3.png')
