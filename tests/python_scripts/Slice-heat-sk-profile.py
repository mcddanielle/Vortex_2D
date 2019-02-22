#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')   #for batch jobs!!!

import os, math, sys, fnmatch, csv
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from pylab import *

import re, numpy

file="/home/mcdermott/Codes-Scripts/Vortex_2D/tests/" + sys.argv[1]
outfile = "heat" + sys.argv[1].split('/')[0] + ".eps"
#define figure before subroutines
#standard figure
fig = plt.figure(figsize=(9, 9))

#gridded figure
G = gridspec.GridSpec(3,3)#, wspace=0.0, hspace=0.0)

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
def make_single_plot(file):
     #get s(k) image
    image_file = file+'/testplot3.png'
    image = plt.imread(image_file)

    x=[]; y=[]; z=[];
    histogram_bin=[]; histogram_freq=[];
    prof_xbin=[]; prof_x=[]; prof_y=[];
    
    #---------------------------Read in Data---------------------------------
    for data in import_text(file+"/xy_output", '\t'):  #call CSV reader here
        x.append( float(data[0]) )
        y.append( float(data[1]) )
        z.append( 100.0*float(data[2]) )

    for data in import_text(file+"/histogram_vortex_amplitude", '\t'):  
        pinning_force = float(data[0])
        histogram_bin.append(float(data[1]))
        histogram_freq.append(float(data[2]))

    for data in import_text(file+'/output', '\t'):
        prof_xbin.append( float(data[0]) )
        prof_x.append( float(data[2]) )
        prof_y.append( float(data[1]) )

    xi = np.linspace(min(x), max(x))
    yi = np.linspace(min(y), max(y))
    
    X, Y = np.meshgrid(xi, yi)
    Z = griddata(x, y, z, xi, yi)

    #################
    # First subplot - color map of vortex stamp
    #################

    ax = subplot(G[0:2,1:3])
    
    ax.set_aspect('equal')
    cset = ax.contourf(X, Y, Z, zdir='z', offset=-1, cmap=cm.hot)

    ax.set_xlim(0, 36.5)
    ax.set_ylim(0, 36.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ########################
    #Second Subplot- profile in y
    ########################
    '''
    ax2 = subplot(G[0:2,0])
    ax2.plot(prof_x,prof_xbin)
    ax2.set_ylim(0, 36.5)
    ax2.set_ylabel('Y',size="large")
    ax2.set_yticks([0,36.5])
    min_x=min(prof_y)
    max_x=max(prof_y)
    ax2.set_xticks( [min_x,max_x] )
    '''
    ax2 = subplot(G[0:2,0])
    cset = ax2.contour(X, Y, Z, zdir='x', offset=-8)

    ########################
    #Third Subplot- profile in x
    ########################
    ax3 = subplot(G[2,1:3])
    ax3.plot(prof_xbin,prof_y)
    ax3.set_ylabel('X',size="large")
    ax3.set_xlim(0, 36.5)
    ax3.set_xticks([0,36.5])
    min_y=min(prof_x)
    max_y=max(prof_x)
    ax3.set_yticks( [min_y,max_y] )
    
    ########################
    #Fourth Subplot- insert s(k) as a subplot
    ########################
    
    ax2 = subplot(G[-1,0])
    ax2.set_aspect('equal')
    ax2.imshow(image)
    ax2.axis('off') # clear x- and y-axes
    ax2.set_xlabel('S(k)')


make_single_plot(file)

fig.tight_layout()

plt.savefig(outfile) #another cmd line opportunity

#plt.show()

