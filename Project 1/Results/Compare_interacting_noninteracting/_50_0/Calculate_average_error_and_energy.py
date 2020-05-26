import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import glob, os


import pandas as pd

"Copy into folder containing energies for a range of alpha values to calculate average sigma"
"Core of program made by Øyvind S. Schøyen and Sebastian G. Winther-Larsen"




PERCENTILES = np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272,\
        16.811894, 18.475307, 20.090235, 21.665994, 23.209251,\
        24.724970, 26.216967, 27.688250, 29.141238, 30.577914,\
        31.999927, 33.408664, 34.805306, 36.190869, 37.566235,\
        38.932173, 40.289360, 41.638398, 42.979820, 44.314105,\
        45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

def bootstrap(data, bootstrap_samples):
    boot_vector = np.zeros(bootstrap_samples)
    for i in range(bootstrap_samples):
        boot_vector[i] = np.average(np.random.choice(data, len(data)))
    return np.var(boot_vector), np.std(boot_vector)

def blocking(data, tol=1e-12):
    n = len(data)
    divisor = int(np.log2(n))
    mu = np.mean(data)
    gamma = np.zeros(divisor)
    variances = np.zeros_like(gamma)

    for i in range(divisor):
        n = len(data)
        gamma[i] = np.sum((data[0:n-1] - mu)*(data[1:n] - mu))/n
        variances[i] = np.var(data)
        data_1 = data[0::2]
        data_2 = data[1::2]
        shortest = len(data_2) if len(data_1) > len(data_2) else len(data_1)
        data = 0.5*(data_1[:shortest] + data_2[:shortest])


    if np.sum(variances**2) <= tol:
        return 0
    factor_1 = (gamma/variances)**2
    factor_2 = (2**np.arange(1, divisor+1))[::-1]
    m = np.cumsum((factor_1*factor_2)[::-1])[::-1]

    for i in range(divisor):
        if (i >= len(PERCENTILES)):
            print ("More data is needed!")
            break
        if (m[i] < PERCENTILES[i]):
            break
    variance = variances[i]/(2**(divisor - i))
    return np.sqrt(variance)


if __name__ == "__main__":
    sum_hsig=0
    i=0
    for file in glob.glob("*.dat"):
        x = np.loadtxt(file)
        hatsig=blocking(x)
        print(file)
        print(" %.2E"%(hatsig))
        sum_hsig+=hatsig
        i+=1
    print("Average h_sigma=%.2E"%(sum_hsig/i))
    with open("averages.txt","a+") as f:
        f.write("Sigma (corrected)     = %.2E"%(sum_hsig/i))




filename="IS_high_"


filename+="Average_energies.csv"
#open for reading
filenumber=0
for file in glob.glob("*.dat"):
    print(file)
    f= open(file)
    if f.mode == 'r':
        lines = f.readlines()
    if filenumber==0:
        n=len(lines)
        E_avr=np.zeros(n)

    for i in range(len(lines)):
        E_avr[i]+=float(lines[i])

    filenumber+=1

E_avr=E_avr/filenumber
f = open(filename, "w+")
for i in range(len(E_avr)):
    f.write("%.6f \n"%E_avr[i])
f.close()



with open("averages.txt","a+") as f:
    f.write("\n")
    f.write("E                     = %.5E"%(np.sum(E_avr)/n))
