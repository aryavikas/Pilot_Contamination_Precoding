from pylab import *
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
import numpy as np
x = [20,30,40,50,60,70]
# dataz-> zeroforcing data
# datam-> mmse precoded data
dataz1=[7.23593,29.4365,49.6867,51.1805,27.5845,31.595] # m=3
datam1=[30.5611,38.5675,44.4086,48.4796,52.2539,55.1769]
dataz2=[28.9302,30.6724,40.3628,48.6436,50.5029,52.4369] # m=6
datam2=[35.8627,43.8871,49.6637,54.0567,57.2512,59.8959]
dataz3=[32.123,33.6198,35.7237,42.0311,43.907,52.0216] #m=9
datam3=[36.947,45.1612,50.7759,54.8983,58.3208,61.3721]
dataz4=[12.5506,31.1651,30.4668,36.1479,36.8643,47.013] #m=2
datam4=[27.1591,35.4426,40.5092,45.6732,48.2274,51.4741]

#line1, = plt.plot(x,dataz1, 'ro', label="$m=3$", linestyle=':',lw=3, ms=9.0)
#line2, = plt.plot(x,datam1, 'ro--', label="$f_DT_s: 0.005-0.05$", lw=3, ms=9.0)
#line2, = plt.plot(x,datam1, 'ro', label="$MMSE Precoder$",linestyle='-' ,lw=3, ms=9.0)
line3, = plt.plot(x,dataz2, 'g^', label="$m=6$",linestyle=':', lw=3, ms=9.0)
line4, = plt.plot(x,datam2, 'g^', label="$MMSE Precoder$", linestyle='-',lw=3, ms=12.0)
#line5, = plt.plot(x,dataz3, 'b*-', label="$f_DT_s: 0.1-0.2$", lw=3, ms=9.0)
line5, = plt.plot(x,dataz3, 'b*', label="$Zero Forcing$",linestyle=':', lw=3, ms=9.0)
line6, = plt.plot(x,datam3, 'b*', label="$m=9$",linestyle='-', lw=3, ms=9.0)
line7, = plt.plot(x,dataz4, 'm*', label="$m=2$",linestyle=':', lw=3, ms=9.0)
line8, = plt.plot(x,datam4, 'm*', label="$m=2$",linestyle='-', lw=3, ms=9.0)


#first_legend = plt.legend(handles=[line1,line3,line6], loc='upper right',bbox_to_anchor=(1.0,0.45))
first_legend = plt.legend(handles=[line7,line3,line6], loc=1)

ax = plt.gca().add_artist(first_legend)

#second_legend = plt.legend(handles=[line2,line5], loc='center')
plt.legend(handles=[line4,line5], loc='best')

plt.xlabel('Number of antennas')
plt.ylabel('Average Sum Rate [bps/Hz]')

plt.grid()
plt.show()

