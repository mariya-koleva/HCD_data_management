# File: HCD_log_processing.py
# Brief: Script for reading and processing chemical concentrations from log files and groups of log files
# Author: Tashi Wischmeyer
# Date Created:
# Date Updated: 07-29-2020

# libraries in use
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as sps
import statsmodels.api as sm
import sklearn.metrics as sk
import csv
import os

# determine local directory and navigate to the Data subdirectory
dir = os.path.dirname(__file__)
data_dir = os.path.join(dir, 'Data')
entries = os.listdir(data_dir)      # list of all files/folders in Data subdirectory
print('Data files:')

while True:
    try:
        # list all entries
        for entry in entries:
            print(entry)
        print('\n')
        # get user input to choose file/folder
        file_name = input('Enter one of the above files or folders to analyze:\t')
    except:
        # print error and continue if invalid
        if file_name not in entries:
            print('Error: File or folder not recognized')
            continue
    else:
        break

# determine new file path
file_path = os.path.join(data_dir, file_name)

hcd_raw_log = []
# if the file path point to a file, only one file will be processed
if os.path.isfile(file_path):
    # open file and read in all rows
    hcd_log_file = open(file_path, 'r')
    hcd_raw_rows = hcd_log_file.read().splitlines()
    hcd_log_file.close()
    # parse each row and split by tab delimiter to get elements
    for row in hcd_raw_rows:
        if row:
            row_parse = [col.strip() for col in row.split('\t') if col]
            hcd_raw_log.append(row_parse)
# if the file path is a directory, multiple files must be processed
elif os.path.isdir(file_path):
    data_files = os.listdir(file_path)      # list files in directory
    # process each file and read all rows
    for i in range(len(data_files)):
        # determine each file path
        data_file_path = os.path.join(file_path, data_files[i])
        # open file and read in all rows
        hcd_log_file = open(data_file_path, 'r')
        hcd_raw_rows = hcd_log_file.read().splitlines()
        hcd_log_file.close()


        head = 0    # indicates if the header has been deleted
        # parse each row and split by tab delimiter to get elements
        for row in hcd_raw_rows:
            if row:
                row_parse = [col.strip() for col in row.split('\t') if col]
                hcd_raw_log.append(row_parse)
            # delete header from all files but the first
            if (head == 0) and (i > 0):
                del hcd_raw_log[-1]
                head = 1
else:
    print('Error: Path not recognized.')
    exit()

# convert array to a DataFrame with the first row as the header
hcd = pd.DataFrame(hcd_raw_log[1:], columns = hcd_raw_log[0])
# get number of row and column elements
r_count, c_count = hcd.shape

while True:
    try:
        # get user input for target contaminant
        target_contaminant = input('What contaminant is the target for this test?\t')
    except ValueError:
        # error if invalid input
        print('I don''t have the capability to deal with that contaminant.')
        exit()
    else:
        break
    
while True:
    try:
        # get user input for target contaminant limit from SAE J2719 standards
        contaminant__stdlimit = input('What is the SAE J2719 limit for the target contaminant, e.g. `x` ppm?\t')
    except ValueError:
        # error if invalid input
        print('I don''t have the capability to deal with that limit.')
        exit()
    else:
        break

# determine input concentration and processing variables from col name
contaminant_name = target_contaminant+"(PPM)"
if contaminant_name in hcd.columns:
    if(target_contaminant == "CO2"):
        input_concentration = [0, 20, 10, 5, 4, 3, 2, 1.8, 1.4, 1.2, 1]
        span = 5
        perc = 0.5
    if(target_contaminant == "CH4"):
        input_concentration = [0, 20, 10, 5, 4, 3, 2, 1.8, 1.4, 1.2, 1]
        span = 5
        perc = 0.5
    elif(target_contaminant == "HCHO"):
        input_concentration = [0, 5.7, 3.04, 2.09, 1.045, 0.513, 0.304]
        span = 5
        perc = 0.45
    elif(target_contaminant == 'H2O'):
        input_concentration = [0, 50, 25, 12.5, 10, 7.5, 5, 4.5, 3.5, 3, 2.5]
        span = 5
        perc = 0.5
    elif(target_contaminant == 'NH3'):
        input_concentration = [0, 1, 0.7, 0.3, 0.1, 0.08, 0.07, 0.05, 0.03]
        span = 5
        perc = 0.5
    elif(target_contaminant == 'CO'):
        input_concentration = [0, 2, 1, 0.5, 0.4, 0.3, 0.2, 0.18, 0.14, 0.12, 0.1]
        span = 5
        perc = 0.5
else:
    # error if contaminant is not in the file
    print('Target contaminant not found')
    exit()

# increment through contaminant name to determine data column
col = 0
for name in hcd.columns:
    if name == contaminant_name:
        break
    col = col + 1
# convert column to floating point numerical data
hcd[contaminant_name] = pd.to_numeric(hcd[contaminant_name], errors='coerce')

# initialize iterative objects
concentration = []
cycle_vals = []
val_index = []
cycle_index = []
row = 30
i = 0
while(row < (r_count-span)):
    # determine current calculation
    calc = np.round(np.mean(hcd.iloc[(row-30):(row-14), col]), decimals = 4, out = None)
    next_conc = input_concentration[(i+1)%len(input_concentration)]
    # experimental difference of multi-pt averages
    exp_diff = np.abs(np.mean(hcd.iloc[(row-span):row, col]) - np.mean(hcd.iloc[(row+span):(row+2*span), col]))
    # threshold for theoretical difference between concentrations
    theo_diff = np.abs((input_concentration[i]-next_conc))

    # look for concentration change
    # should be lower than previous value and higher than next value
    if (exp_diff > theo_diff*perc) and ((i < 2) or ((calc < input_concentration[i-1]*0.2 + input_concentration[i]*0.8) and (calc > next_conc*0.8))):
        # if we have negative values, zero them because it doesn't make sense to have negatives
        if(calc < 0):
            calc = 0
        cycle_vals.append(calc)
        # at the end of a full cycle, sort the concentration averages in descending order
        if(i == len(input_concentration)-1):
            print(cycle_vals)
            cycle_vals.sort(reverse = True)
            concentration.append(cycle_vals)
            cycle_vals = []
        i = (i + 1)%len(input_concentration)
        row = row + 4*span
    else:
        row = row + 1
if not concentration:
    print(cycle_vals)
    print("No cycles detected")
    exit()
# if there are not the correct number of entries, the last non-NULL row value is used to calculate the final entry
elif(len(concentration[-1]) < len(input_concentration)):
    calc = np.round(np.mean(hcd.iloc[(row-30):(row-14), col]), decimals = 4, out = None)
    # if we have negative values, zero them because it doesn't make sense to have negatives
    if(calc < 0):
        calc = 0
    val = concentration[-1]

    val.append(calc)
    concentration[-1] = val


input_concentration.sort(reverse = True)
# average
avg_response_concentration = np.mean(concentration, axis = 0)
# std dev
std_response_concentration = np.std(concentration, axis = 0)
# var
var_response_concentration = np.var(concentration, axis = 0)
# std error
se_response_concentration = sps.sem(concentration, axis = 0)

title_string_0 = 'FTIR Raw Data Curve for ' + target_contaminant
# plotting
area = np.pi*3
contaminant__stdlimit_n = float(contaminant__stdlimit)
twos = contaminant__stdlimit_n*np.ones(21, dtype = int)
# fig, figCalib = plt.subplots(2,2)
figCalib = plt.figure(figsize=(10, 10))
ax1 = figCalib.add_subplot(221)
plt.plot(input_concentration, avg_response_concentration, alpha = 0.5, c = 'b', marker = '.')
plt.plot(range(21), twos, linestyle = ':', c = 'k')
plt.ylim(ymax = max(input_concentration)+0.1*max(input_concentration), ymin = 0)
plt.xlim(xmax = max(input_concentration)+0.1*max(input_concentration), xmin = -0.5)
ax1.set_title(title_string_0)
ax1.set_xlabel('Input Concentration [ppm]', fontsize = 8)
ax1.set_ylabel('Response Concentration [ppm]', fontsize = 8)
ax1.errorbar(input_concentration, avg_response_concentration, yerr = se_response_concentration, c = 'b', linestyle = 'none')
plt.legend(['Average Response', contaminant__stdlimit + ' ppm', 'Standard Error'], fontsize = 8)


title_string_1 = 'FTIR Zero-Adjusted Curve for ' + target_contaminant
# plotting
ax2 = figCalib.add_subplot(222)
zeroed_avg_response = np.array(avg_response_concentration)
zeroed_avg_response = zeroed_avg_response - avg_response_concentration[-1]
plt.plot(input_concentration, zeroed_avg_response, alpha = 0.5, c = 'b', marker = '.')
plt.plot(range(21), twos, linestyle = ':', c = 'k')
plt.ylim(ymax = max(input_concentration)+0.1*max(input_concentration), ymin = 0)
plt.xlim(xmax = max(input_concentration)+0.1*max(input_concentration), xmin = -0.5)
ax2.set_title(title_string_1)
ax2.set_xlabel('Input Concentration [ppm]', fontsize = 8)
ax2.set_ylabel('Response Concentration [ppm]', fontsize = 8)
ax2.errorbar(input_concentration, zeroed_avg_response, yerr = se_response_concentration, c = 'b', linestyle = 'none')
plt.legend(['Average Response', contaminant__stdlimit + ' ppm', 'Standard Error'], fontsize = 8)

# linear regression using OLS from statsmodels.OLS
avg_conc = pd.DataFrame(avg_response_concentration)
conc = pd.DataFrame(concentration)
conc = pd.DataFrame.transpose(conc, copy = False)
in_conc = pd.DataFrame(input_concentration)
model_OLS = sm.OLS(in_conc, avg_conc)
results_OLS = model_OLS.fit()
model_attributes = results_OLS.summary();
test_conc = sm.add_constant(conc[0])
predicted_OLS = results_OLS.predict(avg_conc)
r2_OLS = round(results_OLS.rsquared, 5)
params_OLS = np.array(results_OLS.params)

title_string_2 = 'FTIR Ordinary Least Squares Calibration for ' + target_contaminant
# plotting for OLS
ax3 =  figCalib.add_subplot(223)
x = len(concentration)
for i in range(x):
    plt.scatter(conc[i], results_OLS.predict(conc[i]), s = area, alpha = 0.5, c = 'k')
plt.plot(avg_conc, predicted_OLS, c = 'b')
plt.plot(range(21), twos, linestyle = ':', c = 'k')
plt.ylim(ymax = max(input_concentration)+0.1*max(input_concentration), ymin = 0)
plt.xlim(xmax = max(input_concentration)+0.1*max(input_concentration), xmin = -0.5)
ax3.set_title(title_string_2)
ax3.text(3, 1, 'R2 = ' + str(r2_OLS), fontsize = 8)
ax3.text(3, 1.5, 'y=' + str(round(params_OLS[0], 5)) + '*x+' + str(0), fontsize = 8)
ax3.set_ylabel('Expected Concentration [ppm]', fontsize = 8)
ax3.set_xlabel('Response Concentration [ppm]', fontsize = 8)
plt.legend(['Expected Response', contaminant__stdlimit + ' ppm', 'Cycle Values'], fontsize = 8)

# linear regression with OLS from numpy.polyfit
model_polyfit = np.polyfit(avg_response_concentration, input_concentration, 1)
predict_polyfit = np.poly1d(model_polyfit)

r2_polyfit = round(sk.r2_score(input_concentration, predict_polyfit(avg_response_concentration)), 5)

title_string_3 = 'FTIR Polyfit Calibration for ' + target_contaminant
# plotting for polyfit
ax4 =  figCalib.add_subplot(224)
x = len(concentration)
for i in range(x):
    plt.scatter(concentration[i], predict_polyfit(concentration[i]), s = area, alpha = 0.5, c = 'k')
plt.plot(avg_response_concentration, predict_polyfit(avg_response_concentration), c = 'b')
plt.plot(range(21), twos, linestyle = ':', c = 'k')
ax4.set_title(title_string_3)
ax4.text(3, 1, 'R2 = ' + str(r2_polyfit), fontsize = 8)
ax4.text(3, 1.5, 'y=' + str(round(model_polyfit[0], 5)) + '*x+' + str(round(model_polyfit[1], 5)), fontsize = 8)
ax4.set_ylabel('Expected Concentration [ppm]', fontsize = 8)
ax4.set_xlabel('Response Concentration [ppm]', fontsize = 8)
plt.ylim(ymax = max(input_concentration)+0.1*max(input_concentration), ymin = 0)
plt.xlim(xmax = max(input_concentration)+0.1*max(input_concentration), xmin = -0.5)
plt.legend(['Expected Response', contaminant__stdlimit + ' ppm', 'Cycle Values'], fontsize = 8)

plt.tight_layout(pad = 1, h_pad = 0.75, w_pad = 0.75)

# save figure
figCalib.savefig(dir + '//' + 'Results\Calibration_curve_' + target_contaminant + '.pdf')

plt.show()
