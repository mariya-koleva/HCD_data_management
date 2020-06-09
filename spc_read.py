# will read in directly from spc once library problems are resolved
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


hcd_csv_raw = pd.read_csv('C:/Users/Tashi Wischmeyer/Documents/HCD_data_management/Output_spectra.csv', header=None)

r_count, c_count = hcd_csv_raw.shape

contaminant_absorptivity = [[0.25, 1995], [0.43, 2013], [0.43, 2060]]
pathlength = 0.20

hcd_spectra_raw = []
for i in range(r_count):
    temp_str = str(hcd_csv_raw.iloc[i])
    temp_str1 = list(temp_str.split())
    vec = list(map(float, temp_str1[1].split('\\t')))
    hcd_spectra_raw.append(vec)

hcd_spectra = pd.DataFrame(hcd_spectra_raw)
# figure
plt.plot(hcd_spectra.iloc[:,0], hcd_spectra.iloc[:,1])
plt.show()
