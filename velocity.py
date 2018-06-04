#for l=6
from pylab import *
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
import numpy as np
x = [20,30,40,50,60,70]
# dataz-> zeroforcing data
# datam-> mmse precoded data
dataz1=[28.9302,30.6724,40.3628,48.6436,50.5029,52.4369]
datam1=[35.8627,43.8871,49.6637,54.0567,57.2512,59.8959]
dataz2=[16.1284,23.518,33.2544,35.7522,38.3449,44.2247]
datam2=[22.1208,29.725,35.6455,39.5561,43.3205,46.191]
dataz3=[4.39814,5.25755,18.974,22.7756,24.4496,26.6526]
datam3=[0.646878,1.35403,2.26649,2.69753,4.09157,4.25345]


line1, = plt.plot(x,dataz1, 'ro', label="$f_DT_s: 0.005-0.05$", linestyle=':',lw=3, ms=9.0)
#line2, = plt.plot(x,datam1, 'ro--', label="$f_DT_s: 0.005-0.05$", lw=3, ms=9.0)
line2, = plt.plot(x,datam1, 'ro', label="$MMSE Precoder$",linestyle='-' ,lw=3, ms=9.0)
line3, = plt.plot(x,dataz2, 'g^', label="$f_DT_s: 0.05-0.1$",linestyle=':', lw=3, ms=9.0)
line4, = plt.plot(x,datam2, 'g^', label="$f_DT_s: 0.05-0.1$", linestyle='-',lw=3, ms=12.0)
#line5, = plt.plot(x,dataz3, 'b*-', label="$f_DT_s: 0.1-0.2$", lw=3, ms=9.0)
line5, = plt.plot(x,dataz3, 'b*', label="$Zero Forcing$",linestyle=':', lw=3, ms=9.0)
line6, = plt.plot(x,datam3, 'b*', label="$f_DT_s: 0.1-0.2$",linestyle='-', lw=3, ms=9.0)

#first_legend = plt.legend(handles=[line1,line3,line6], loc='upper right',bbox_to_anchor=(1.0,0.45))
first_legend = plt.legend(handles=[line1,line3,line6], loc=1)

ax = plt.gca().add_artist(first_legend)

#second_legend = plt.legend(handles=[line2,line5], loc='center')
plt.legend(handles=[line2,line5], loc=4)

plt.xlabel('Number of antennas')
plt.ylabel('Average Sum Rate [bps/Hz]')

plt.grid()
plt.show()






