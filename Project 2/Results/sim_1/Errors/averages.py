import numpy as np
import matplotlib.pyplot as plt
import string as str
from matplotlib.pyplot import figure
import matplotlib as mpl
import glob, os


from matplotlib import rc

file_number=0
for file in glob.glob("*.dat"):
        print(file)
        f= open(file)
        if f.mode == 'r':
            lines = f.readlines()
        if file_number==0:
            n=len(lines)
            E_avr=np.zeros(n)

        for i in range(len(lines)):
            E_avr[i]+=float(lines[i])
        file_number+=1

filename="average_error"
E_avr=E_avr/file_number
f = open(filename, "w+")
for i in range(len(E_avr)):
    f.write("%.6f \n"%E_avr[i])
f.close()
