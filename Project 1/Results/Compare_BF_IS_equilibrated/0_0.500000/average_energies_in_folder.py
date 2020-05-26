import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import glob, os

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
    f.write("\\n")
    f.write("E      = %.2E"%(np.sum(E_avr)/n))
