#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')   #for batch jobs!!!

import os, math, sys, fnmatch, csv, glob
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec

from pylab import *

import re, numpy

#define figure before subroutines
#standard figure
fig = plt.figure()

#gridded figure
G = gridspec.GridSpec(2, 4,wspace=0.0, hspace=0.0)
#G.update(left=0.05, right=1.0, wspace=0.05)
#G.update(hspace=0.05)

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

def get_coordinates(file_name,plot_position,letter):    
    pos_x = []
    pos_y = []

    for data in import_text(file_name+"/ascii_0.000000", ' '):  
        pos_x.append(float(data[1]))
        pos_y.append(float(data[2]))

    ax = subplot(plot_position)        
    ax.set_aspect('equal')
    ax.scatter(pos_x,pos_y,s=0.5)

    ax.tick_params(axis='both', which='major', labelsize='small')
    ax.set_xlabel('X', fontsize='small', labelpad=-5)
    ax.set_xlim(0, 36.5)
    ax.set_xticks([0, 36.5])#, size='small')
    ax.set_ylabel('Y', fontsize='small', labelpad=-15)
    ax.set_ylim(0, 36.5)
    ax.set_yticks([0, 36.5])#, size='small')
    ax.text(1.5,36.5-3,"("+letter+")",backgroundcolor='white')
    #print "plotted " + file_name + "\n"
    return
#---------------------------------------------------------------------------
#nominal main.c-------------------------------------------------------------
#---------------------------------------------------------------------------

#name files:
#string_files="/home/mcdermott/Codes*/Vortex_2D/tests/N160*" 
data_files=[]
for fp in [0.0, 0.6, 1.0, 1.2, 1.5, 2.0, 2.4, 4.0]:
     data_files.append("/pscratch/dmcderm2/LANL/Attract_Repel/N294_rho_0.22"+"/N294_Np1_fp"+fp)
'''
data_files = glob.glob(string_files)
#sort by pinning force
data_files = sorted(data_files,key=lambda x: float( x.split('Fp')[-1]) )
'''
#---------------------------Read in Data---------------------------------

a=0
b=0
letter = ['a','b','c','d','e','f','g','h']
i=0

for file in data_files:
    plot_position = G[a,b]
    #print a, b
    get_coordinates(file,plot_position,letter[i])
    '''
    if a < 1:
        a = a+1
    elif b < 3:
        a = 0
        b = b+1
    '''
    if b < 3:
        b = b+1
    elif a < 1:
        b = 0
        a = a+1
    else:
        break
    i = i+1

G.tight_layout(fig, rect=[0, 0, 1, 1], h_pad=0.0, w_pad=-0.5,pad=0)

#figure(num=None, figsize=(4,2), dpi=80, facecolor='w', edgecolor='k')
fig.set_size_inches(15,7.2)
plt.savefig("ascii.eps") #another cmd line opportunity

#plt.show()

