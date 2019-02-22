#!/usr/bin/python

"""General_Vortex2D_PlotModule.py
Danielle McDermott
May 25, 2014

Collection of subroutines designed for plotting 
2D position data and making figures with Gridspec. 
The data files were generated by code Vortex2D.
Takes as a 'main' function Fig1_400.py, 
so by default it will create Figure 1 of 'Stripes on Stripes' paper
written in May 2014

Some notes:
Currently working to make a more general version of 
ascii_multiplot.py

ascii plot and heat plot have a lot in common,
currently combined into a single subroutine
"""

#######################################################
#Import Python Library Modules
#######################################################
import matplotlib
matplotlib.use('Agg')   #for batch jobs!!!

import os, math, sys, csv, glob, numpy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#contour plots
from matplotlib import cm
#3D plots
#from mpl_toolkits.mplot3d import axes3d


#####################################################
#Define Modules
#####################################################
def import_text(filename, separator):
    '''Simple csv file reader, reads line by line, ignores comments of '#' type

    file_name = absolute name
    separator = ' ', '\t', etc
    '''
    for line in csv.reader(open(filename), delimiter=separator, 
                           skipinitialspace=True):
        if line:
            if line[0].startswith("#"):
                print ""
            elif line:
                #print line
                yield line
    return
#------------------------------------------------------------
#end CSV reader------------------------------------------
#------------------------------------------------------------

def letter_range(start, stop):
    '''populate a list by calling ord(a) to ord(f), or some such.  

    if start > stop, loop will not execute

    start = r
    stop = z
    '''

    for c in xrange(ord(start), ord(stop)):
        yield chr(c)

    return

#------------------------------------------------------------
#end letter_range()-----------------------------------------
#------------------------------------------------------------

def add_subplot_axes(ax,rect,axisbg='w'):
    '''Used to add a subplot to ax (which may itself be a subplot) of size rect

    ax = pointer to subplot
    rect = [x0,y0,x_size,y_size (percent)]
    '''

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

#------------------------------------------------------------
#end add_subplot_axes---------------------------------------
#------------------------------------------------------------

def add_gr_plot(file_name,ax,xlim1=[0.5, 4.0]):
    '''plot g(r) for a given file_name on a given ax

    file_name = absolute file name
    ax = pointer to a subplot
    '''

    ##########################
    #Get Data
    #########################
    r = []
    g_r=[]

    #call CSV reader here
    for data in import_text(file_name, ' '):  
        r.append(    float(data[0]) )
        g_r.append( float(data[1]) )

    ############################
    #Plot data
    ############################
    ###ax = subplot(plot_position)        
    ax.plot(r,g_r)
    ax.set_xlabel('r', fontsize='small', labelpad=0)
    ax.set_ylabel('g(r)', fontsize='small', labelpad=0)
    ax.set_xlim(*xlim1) 
    plt.setp(ax, xticks=[1.0,2.0,3.0,4.0],
             yticks=[0,1,2,3,4,5,6,7])

    return

#------------------------------------------------------------
#end add_gr_plot
#------------------------------------------------------------

def add_sk_image(file_name,ax,column=0,labels=0):
    '''add an image from file_name to a subplot ax

    file_name = '$HOME/all_files/testplot3.png'
    ax = pointer to a subplot
    '''

    #########################################
    #get s(k) image
    #########################################
    image_file = file_name #+'/testplot3.png'
    image = plt.imread(image_file)

    ax.set_aspect('equal')
    ax.imshow(image)       #plot image
    if labels:
        ax.set_xlabel('$k_x$')  #turned off
        ax.set_ylabel('$k_y$')  #turned off
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.axis('off')         # clear x- and y-axes

    return

#------------------------------------------------------------
#end add_sk_image()------------------------------------------
#------------------------------------------------------------

def get_coordinates(file_name,ax,labels=1,Sx=[0,36.5],Sy=[0,36.5],num_dens=0, heat=0):    
    '''
    Plot coordinates of given file, default particle positions, 
    optional: heat map of vortex density

    file_name is an absolute path, in earlier versions it was a directory
    file_name = "$HOME/all_files/ascii_0.00000" or 
    file_name = "$HOME/all_files/xy_output"

    ax is pointer to a subplot of the figure, i.e
    ax = fig.add_subplot(grid_position) 

    heat is a flag to ignore particle
    ascii plot and heat plot have a lot in common,
    currently combined into a single subroutine
    heat = 0 
    heat = 1
    '''
    ###################################################
    #Import Data
    ###################################################
    if heat == 0:        #pure particle positions
        pos_x = []
        pos_y = []

        for data in import_text(file_name, ' '):  
            pos_x.append(float(data[1]))
            pos_y.append(float(data[2]))

        #################################################
        #Calculate effective number density (for one pin)
        #################################################
        if num_dens == 1:
            #find the x values centered about the single pin minima
            x_c = []
            for x in pos_x: 
                if abs(x-0.75*Sx[1]) > Sx[1]/2.0:
                    x += Sx[1]
                    x_c.append(x)

            x_min = min(x_c)
            x_max = max(x_c)

            effective_area = (x_max - x_min)*Sy[1]

            eff_num_density = len(pos_x)/effective_area
            print "effective density of", file_name.split("/")[-1], ":", eff_num_density
        #end num_dens if
        
    else:
        #####################################
        #Heat Data
        #####################################
        x=[]; 
        y=[]; 
        z=[];
        data_file = file_name 

        for data in import_text(data_file, '\t'):  #call CSV reader here
            x.append( float(data[0]) )
            y.append( float(data[1]) )
            z.append( 100.0*float(data[2]) )
            
        ############################################
        #Contour Plotting
        ############################################
        #create regular x-y grid
        xi = numpy.linspace(min(x), max(x))
        yi = numpy.linspace(min(y), max(y))
        
        #create a mesh x-y grid from the numpy.linspace
        #these are necessary for the plot
        X, Y = numpy.meshgrid(xi, yi)
    
        #take the Heat (x,y,z) data and grid it on the xi, yi (not the mesh)
        Z = matplotlib.mlab.griddata(x, y, z, xi, yi)
        

    ###########################################################
    #Data has been imported, set up the plot
    ###########################################################
    #ax is declared in main function
    ax.set_aspect('equal')

    if heat == 0:
        #############################################
        #simple x-y plot
        #############################################
        ax.scatter(pos_x,pos_y,s=0.5) 
    else:
        #############################################
        # color map of vortex stamp
        #############################################
        cset = ax.contourf(X, Y, Z, zdir='z', offset=-1, cmap=cm.hot)
    
    ax.set_xlim(*Sx)
    ax.set_ylim(*Sy)

    if labels:
        ax.tick_params(axis='both', which='major', labelsize='small')
        ax.set_xlabel('X', fontsize='small', labelpad=0)
        #ax.set_xticks(*Sx) #, size='small')
        ax.set_ylabel('Y', fontsize='small', labelpad=0)
        #ax.set_yticks(*Sy) #, size='small')
    else:
        plt.tick_params(\
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # labels along the bottom edge are off
        plt.tick_params(\
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',      # ticks along the bottom edge are off
            right='off',         # ticks along the top edge are off
            labelleft='off') # labels along the bottom edge are off


    
    ################################################################
    #Structure Factor -- broken
    #should add a subplot to the ascii file, 
    #try such things in the main function
    ################################################################
    '''
    if sk == 1:
        add_sk_image(file_name,ax,column = 0)
    '''
    #if gr:
    #    add_gr_plot(file_name,ax)

    return
#------------------------------------------------------------
#end plot_particle_positions
#------------------------------------------------------------

    
#---------------------------------------------------------------------------
#nominal main.c-------------------------------------------------------------
#---------------------------------------------------------------------------
if __name__ == "__main__":
    ###################################################
    #Get Parameters 
    #(1) single value of Nv, Fp, Np
    #2x2 figure
    ###################################################
    Sy = [0,36.5]
    Sx = [0,36.5]
    verbose = 0

    ##########################
    #directory to plot assume we are in it
    ##########################
    directory = ""

    ############################################################### 
    #several of these parameters are hardwired latter in the 'loop'
    ###############################################################
    gr = 1             #add g(r), pair correlation data (not used)
    sk = 2             #add s(k), structure factor, data (not used)
    heat = 1           #plot 'vortex density' (used!!!)
    pin_array = [0]    #not used
    force_array = [0]  #not used

    #######################################################
    #define number of row/columns to make a gridded figure
    #######################################################
    rows = 2
    columns = 2

    #########################################################
    ####################parameters got!######################
    #########################################################

    #turn on/off the "X" "Y" and axis ticks
    labels = 1

    #########################################################
    #identify data files
    #########################################################
        
    position_file = directory + "ascii_0.000000"
    heat_file = directory + "xy_output"
    sk_file = directory + 'testplot3.png'

    out_file = "Fig1.eps"
    gr_file ='gr_sigma0.001.dat'
    data_files = [gr_file, position_file, heat_file, sk_file]

    #########################################################
    ####################Data Files ID'd!#####################
    #########################################################

    ###########################
    #Define standard figure
    ###########################
    fig = plt.figure( figsize=(columns*2,rows*2) )
        
    ############################################################
    #created grid for figure (subplots connect them, IMPORTANT)
    ############################################################
    G = gridspec.GridSpec(rows, columns, wspace=0.1, hspace=0.1)
    #G.update(hspace=0.05)
        
    #######################################
    #populate a list for annotations
    #######################################
    letter = []
    for x in letter_range('a', 'z'):
        letter.append( x )  
    
    #########################################################
    #add letter annotation at (xt,yt)
    #########################################################
    xt = Sx[0]+1.5
    #y -position of text depends on whether or not label exists
        
    if labels:
        yt = Sy[1]-5.0        
    else:
        yt = Sy[1]-4.0

    #######################################
    #Iterate through grid adding subplots
    
    #by rights, should be a loop, 
    #but the function called is different every time
    #######################################
    i = 0
    verbose=1
    for a in range(rows):
        if verbose:
            print "a=",a
        for b in range(columns):
            if verbose:
                print "a=",a,"b=",b

            #grab a grid pointer
            plot_position = G[a,b]    

            #make a subplot at that location
            ax = fig.add_subplot(plot_position)        

            #grab file from the declared list
            file = data_files[i]
            if verbose:
                print "i=", i, "file=", file

            ##########################################
            #Case/Switch (ish) to plot (a-d)
            ##########################################
            if i == 0:    #(a) GR plot
                add_gr_plot(file,ax)

            elif i == 1:  #(b) x-y positions plot
                get_coordinates(file,ax)

            elif i == 2:  #(c) vortex density map
                get_coordinates(file,ax,heat=1,labels=1)

            elif i == 3:  #(d) structure factor
                add_sk_image(file,ax,labels=1)
            else:
                print "some error in grid, plots, etc"
                exit(-1)

            ###########################################
            ##Add subplot letter (a), etc
            ##text position is hardwired, needs work
            ###########################################

            xt0 = xt
            yt0 = yt

            ax.text(xt0,yt0,"("+letter[i]+")",backgroundcolor='white',size=16)
        
            i += 1

            ###########################################
            #for loop ends here
            ###########################################


    ########################################
    #Tweak the figure
    ########################################
        
    if labels:
        G.tight_layout(fig, rect=[0, 0, 1, 1], h_pad=0.2, w_pad=0.2,pad=0.2)     
    else:
        p_val = -1.4
        G.tight_layout(fig, rect=[0, 0, 1, 1], h_pad=p_val,w_pad=p_val,pad=0.0)
    ###################################################
    #Save the Figure
    ###################################################
    plt.savefig(out_file) #another cmd line opportunity
    
    exit()