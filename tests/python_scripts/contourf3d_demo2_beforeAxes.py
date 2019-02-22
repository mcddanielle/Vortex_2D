#!/usr/bin/python

"""
.. versionadded:: 1.1.0
   This demo depends on new features added to contourf3d.
"""
import matplotlib
matplotlib.use('Agg')   #for batch jobs!!!

import os, math, sys, fnmatch, csv
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

#---------------------------------------------------------------------------
#nominal main.c-------------------------------------------------------------
#---------------------------------------------------------------------------

#get s(k) image
image_file = 'testplot3.png'
image = plt.imread(image_file)

#standard figure
fig = plt.figure()
#central figure (?)
###ax = fig.gca(projection='3d')

#python test data
#X, Y, Z = axes3d.get_test_data(0.05)

x=[]; 
y=[]; 
z=[];

pos_x = []
pos_y = []

histogram_bin=[]
histogram_freq=[];

r = []
g_r=[]

q6_hist_bin=[]
q6_hist_freq=[]

prof_xbin=[]; prof_x=[]; prof_y=[];
#---------------------------Read in Data---------------------------------
#
#a little hackish, could be done nicely with a hash
#
for data in import_text("xy_output", '\t'):  #call CSV reader here
    x.append( float(data[0]) )
    y.append( float(data[1]) )
    z.append( 100.0*float(data[2]) )

for data in import_text("histogram_vortex_amplitude", '\t'):  
    pinning_force = float(data[0])
    histogram_bin.append(float(data[1]))
    histogram_freq.append(float(data[2]))

for data in import_text("ascii_0.000000", ' '):  
    pos_x.append(float(data[1]))
    pos_y.append(float(data[2]))

for data in import_text('gr_sigma0.020.dat', ' '):  #call CSV reader here
    r.append(    float(data[0]) )
    g_r.append( float(data[1]) )

for data in import_text('Histogram_BondOrder2D.txt', '\t'):  #call CSV reader here
    q6_hist_bin.append( float(data[0]) )
    q6_hist_freq.append(float(data[3]) )

for data in import_text('output', '\t'):
    prof_xbin.append( float(data[0]) )
    prof_x.append( float(data[2]) )
    prof_y.append( float(data[1]) )


#get minima/maxima of histogram
''' #IMPORTANT
extrema are labeled by array indices 
to get the r values and peak heights, 
you must access those arrays with these indices 
i.e. r1 = r[min[0]], g(r1) = g_r[min[0]]
'''
hist_min,hist_max = return_extrema(histogram_freq)

np_histogram_bin = np.array(histogram_bin)
np_histogram_freq = np.array(histogram_freq)

#------------------------------------------------------------------
#for plotting, convert lists to arrays
r_np    = np.array(r)
gr_array  = np.array(g_r)
#

'''
print species_min
print species_max
'''



'''
X_np = numpy.array(X)
Y_np = numpy.array(Y)
Z_np = numpy.array(Z)
'''
xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)


#################
# First subplot - g(r) data
#################
ax = fig.add_subplot(2, 2, 1)
ax.plot(r,g_r)
#ax.plot(np_histogram_bin[hist_max],np_histogram_freq[hist_max], "o")
ax.set_xlabel('particle-particle distance')
ax.set_ylabel('Frequency')
ax.set_xlim(0.5, 4.0)
#add a sub-subplot
rect = [0.4,0.5,0.5,0.5]  #x0,y0,x_size,y_size (percent)
ax1 = add_subplot_axes(ax,rect)
ax1.plot(r,g_r)
ax1.set_yscale('log')
ax1.set_xlim(0.5, 4.0)

#################
# Second subplot
#################
ax = fig.add_subplot(2, 2, 2)
ax.plot(histogram_bin,histogram_freq)
ax.plot(np_histogram_bin[hist_max],np_histogram_freq[hist_max], "o")
ax.set_xlabel('Vortex Amplitude')
ax.set_ylabel('Frequency')

###################
# Inset of Second subplot
##################
#ax = fig.add_subplot(2, 2,)
rect = [0.75,0.55,0.4,0.4]  #x0,y0,x_size,y_size (percent)
ax2 = add_subplot_axes(ax,rect)
ax2.plot(histogram_bin,histogram_freq)
ax2.plot(np_histogram_bin[hist_max],np_histogram_freq[hist_max], "o")
ax2.set_yscale('log')
#ax2.set_xlabel('Vortex Amplitude')
#ax2.set_ylabel('Frequency')

#ax.set_ylim(0, 36.5)

#
################
# Fourth subplot - color map of vortex stamp
#################

#ax = plt.figure(figsize=plt.figaspect(2.))
ax = fig.add_subplot(2, 2, 3)
ax.set_aspect('equal')
cset = ax.contourf(X, Y, Z, zdir='z', offset=-1, cmap=cm.hot)
#cset = ax.contourf(X, Y, Z, zdir='x', offset=-1, cmap=cm.coolwarm)
#cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
###cb = fig.colorbar(ax=ax)#don't work

ax.set_xlabel('X')
ax.set_xlim(0, 36.5)
ax.set_ylabel('Y')
ax.set_ylim(0, 36.5)

####
#Fifth Subplot - real space
####
ax = fig.add_subplot(2, 2, 4)
ax.set_aspect('equal')
#plot(pos_x,pos_y,"o");
ax.scatter(pos_x,pos_y,s=0.5)

ax.set_xlabel('X')
ax.set_xlim(0, 36.5)
ax.set_ylabel('Y')
ax.set_ylim(0, 36.5)


########################
#Sixth Subplot- insert s(k) as a subplot
########################
rect = [-0.38,0.6,0.5,0.5]  #x0,y0,x_size,y_size (percent)
ax2 = add_subplot_axes(ax,rect)
#ax = fig.add_subplot(2, 3, 6)
ax2.set_aspect('equal')
ax2.imshow(image)
ax2.axis('off') # clear x- and y-axes
ax2.set_xlabel('S(k)')
#ax2.set_xlim(0.1,0.9)

########################
#Seventh Subplot- insert q6 as a subplot
########################
rect = [-0.2,0.0,0.3,0.3]  #x0,y0,x_size,y_size (percent)
ax3 = add_subplot_axes(ax,rect)
ax3.plot(q6_hist_bin,q6_hist_freq)
ax3.set_xlabel("q6")
ax3.set_ylim(-0.15,1.0)


fig.tight_layout()

plt.savefig("histogram.eps") #another cmd line opportunity

#plt.show()

