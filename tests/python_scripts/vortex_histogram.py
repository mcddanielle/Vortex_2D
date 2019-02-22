#!/usr/bin/python

"""
.. versionadded:: 1.1.0
   This demo depends on new features added to contourf3d.
"""
import matplotlib
matplotlib.use('Agg')   #for batch jobs!!!

import os, math, sys, fnmatch, csv, glob
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

from pylab import *

import re, numpy


#-----------------------
# CSV reader-----------------------------------------
def import_text(filename, separator):
    for line in csv.reader(open(filename), delimiter=separator, 
                           skipinitialspace=True):
        if line:
            if line[0].startswith("#"):
                print ""
            elif line:
                #print line
                yield line
#end CSV reader------------------------------------------

def return_extrema(data_with_extrema):
    a   =  diff(sign(diff(data_with_extrema))).nonzero()[0] + 1  # local min+max
    min = (diff(sign(diff(data_with_extrema))) > 0).nonzero()[0] + 1 # local min
    max = (diff(sign(diff(data_with_extrema))) < 0).nonzero()[0] + 1 # local max

    return min,max

#end return_extrema-------------------------------------------

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[2]
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

#get one histogram set of data
def get_histogram(file_name):

    histogram_bin=[]; histogram_freq=[];

    #---------------------------Read in Data---------------------------------
    for data in import_text(file+"/histogram_vortex_amplitude", '\t'):  
        pinning_force = float(data[0])
        histogram_bin.append(float(data[1]))
        histogram_freq.append(float(data[2]))
        
    #get minima/maxima of histogram
    ''' #IMPORTANT: extrema are labeled by array indices 
    to get the r values and peak heights, 
    access arrays with these indices 
    i.e. r1 = r[min[0]], g(r1) = g_r[min[0]]'''

    hist_min,hist_max = return_extrema(histogram_freq)
    
    #find first and last transition to nonzero values...
    first_pos = 0;
    last_pos = 0;

    for i in range (0,len(histogram_freq),1):
        if histogram_freq[i] > 0.0:
            first_pos = i
            break
    for i in range (len(histogram_freq)-1,0,-1):
        if histogram_freq[i] > 0.0:
            last_pos = i
            break
    
    total_peak_width = histogram_bin[last_pos]-histogram_bin[first_pos]
    max_to_max_width = histogram_bin[hist_max[-1]] - histogram_bin[hist_max[0]]
    m2m_ratio =  histogram_freq[hist_max[0]] / histogram_freq[hist_max[-1]]

    return [total_peak_width,max_to_max_width, m2m_ratio]

#---------------------------------------------------------------------------
#nominal main.c-------------------------------------------------------------
#---------------------------------------------------------------------------

#instead get many files....
string_files="/home/mcdermott/Codes*/Vortex_2D/tests/N*" 
data_files = glob.glob(string_files)
data_files = sorted(data_files,key=lambda x: float( x.split('Fp')[-1]) )
print data_files
exit()
total_width=[]; m2m_width=[]; m2m_ratio=[]; Fp=[]; Nv=[]; Np=[];

for file in data_files:
    if os.path.isdir(file):

       file_string = file.split('/')[-1]


       num_vortices = file_string.split('_')[0][1:]
       number_pins = file_string.split('_')[1][2:]
       pin_force = file_string.split('_')[2][2:]
       
       total_peak_width,max_to_max_width, max2max_ratio =  get_histogram(file) 
       Fp.append(pin_force)
       Nv.append(num_vortices)
       Np.append(number_pins)
       total_width.append(total_peak_width)
       m2m_width.append(max_to_max_width)
       m2m_ratio.append(max2max_ratio)

#################
# First subplot
#################
fig = plt.figure()

ax = fig.add_subplot(2, 1, 1)
ax.plot(Fp,total_width)
ax.plot(Fp,m2m_width)

ax.set_xlabel('F$_p$')
ax.set_ylabel('Histogram Metrics')

'''
rect = [0.5,0.5,0.5,0.5]  #x0,y0,x_size,y_size (percent)
ax2 = add_subplot_axes(ax,rect)
ax2.plot(histogram_bin,histogram_freq)
ax2.plot(np_histogram_bin[hist_max],np_histogram_freq[hist_max], "o")
ax2.set_yscale('log')
'''
#ax.set_ylim(0, 36.5)

#################
# Second subplot
#################
#
ax = fig.add_subplot(2, 1, 2)
ax.plot(Fp,m2m_ratio)

ax.set_xlabel('F$_p$')
ax.set_ylabel('Histogram Metrics')

################
# Third subplot -
#################

fig.tight_layout()

plt.savefig("Fp_vs_hist_data.eps") #another cmd line opportunity

#plt.show()

