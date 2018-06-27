#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 18:30:32 2018

@author: u1490431
"""

## Working on NAPH analysis for FENS poster

# Distracted or not peaks

TDTfileslist = ['NAPH02_distraction', 'NAPH03_distraction', 'NAPH04_distraction',\
                'NAPH05_distraction', 'NAPH07_distraction', 'NAPH08_distraction',\
                'NAPH09_distraction', 'NAPH10_distraction']

#TDTfileslist = ['NAPH09_distraction']
TDTfilepath = '/Volumes/KP_HARD_DRI/kp259/NAPH1/MATLAB/'

# NAPH09 is 100% distracted (can't analyse non distracted trials for him)


# Licking days 
#NAPH02_licktrain.mat
#NAPH03_licktrain.mat
#NAPH04_licktrain.mat
#NAPH05_licktrain.mat
#NAPH07_licktrain.mat
#NAPH08_licktrain.mat
#NAPH09_licktrain.mat
#NAPH10_licktrain.mat


# Assign empty lists for storing arrays of burst/run lengths
allDistracted = []
allNotDistracted = []
allDistractors = []
allRatLicks = []
allRatBlue = []
allRatUV = []
allRatFS = []
# Loop through files and extracts all info needed for snips 
for filename in TDTfileslist:
    
    file = TDTfilepath + filename
    ratdata = loadmatfile(file)
    distracted = ratdata['distracted']
    notdistracted = ratdata['notdistracted']
    distractors = ratdata['distractors']
        
    allRatBlue.append(ratdata['blue'])
    allRatUV.append(ratdata['uv'])
    allRatFS.append(ratdata['fs'])
    
    allRatLicks.append(ratdata['licks'])
    allDistracted.append(distracted)
    allNotDistracted.append(notdistracted)
    allDistractors.append(distractors)

# ===============================================================================
# ===============================================================================
# ===============================================================================
# ===============================================================================

# Calculating different peaks, and means based on on DISTRACTED and NOT DISTRACTED
# To find peaks and averaged activity across time 

# Called means because give 14 arrays which are the MEAN for each rat 
allblueMeans = []
alluvMeans = []
for i, val in enumerate(allDistracted):

        # make a blue and uv snip for all 14
        blueSnips, ppsBlue = snipper(allRatBlue[i], allDistracted[i], fs=allRatFS[i], bins=300)
        uvSnips, ppsUV = snipper(allRatUV[i], allDistracted[i], fs=allRatFS[i], bins=300)
# # these four lines used later to define means plot (made after runs) 
        # Makes a mean for each rat's snips
        blueMean = np.mean(blueSnips, axis=0)
        allblueMeans.append(blueMean)
        uvMean = np.mean(uvSnips, axis=0)
        alluvMeans.append(uvMean)

        
#values for averaged activity for the 2 seconds preceeding (14 values, 1 per rat) =

#TwoSecBEFOREactivity = (np.mean(allblueMeans, axis=1)) # axis 1 gives 14 values, axis 0 gives 1 list of 20 (the averaged average snip)

# Want to produce slices of the lists already here
# I have 14 lists of 300. 0 to 100 are the 10 seconds preding distraction
# 100 to 300 are 20 seconds after 

# Split into (1) 80 - 100 (the 2 seconds before distraction)
           # (2) 110 - 300 (the 19 seconds after - to see if supression)

# JUST BLUE ACTIVITY IGNORING THE UV (NO SUBTRACTION AS PRETTY MUCH ZERO)



all2SecBefore = []
all20SecAfter = []
for eachlist in allblueMeans:
    slice1 = eachlist[80:100]
    #print(len(slice1))
    
    slice2 = eachlist[100:300]
    #print(len(slice2))
    
    all2SecBefore.append(slice1)
    all20SecAfter.append(slice2) 
    
    #then out of the loop mean these (axis 1 not 0) = gives 14 points for the barscatter

      
MeanOf2SecDISTRACTED = (np.mean(all2SecBefore, axis=1))
MeanOf20SecDISTRACTED = (np.mean(all20SecAfter, axis=1))

# Repeat for not distracted - reassignes all values 

allblueMeans = []
alluvMeans = []
for i, val in enumerate(allNotDistracted):

        # make a blue and uv snip for all 14
        try:
            
            blueSnips, ppsBlue = snipper(allRatBlue[i], allNotDistracted[i], fs=allRatFS[i], bins=300)
            uvSnips, ppsUV = snipper(allRatUV[i], allNotDistracted[i], fs=allRatFS[i], bins=300)
        except:
            pass
# # these four lines used later to define means plot (made after runs) 
        # Makes a mean for each rat's snips
        blueMean = np.mean(blueSnips, axis=0)
        allblueMeans.append(blueMean)
        uvMean = np.mean(uvSnips, axis=0)
        alluvMeans.append(uvMean)


all2SecBefore = []
all20SecAfter = []
for eachlist in allblueMeans:
    slice1 = eachlist[80:100]
    #print(len(slice1))
    
    slice2 = eachlist[100:300]
    #print(len(slice2))
    
    all2SecBefore.append(slice1)
    all20SecAfter.append(slice2) 
    
MeanOf2SecNOT_DISTRACTED = (np.mean(all2SecBefore, axis=1))
MeanOf20SecNOT_DISTRACTED = (np.mean(all20SecAfter, axis=1))

PEAK2SEC = np.empty((2,), dtype=np.object)
PEAK2SEC[0] = MeanOf2SecDISTRACTED 
PEAK2SEC[1] = MeanOf2SecNOT_DISTRACTED

peak20sec = np.empty((2,), dtype=np.object)
peak20sec[0] = MeanOf20SecDISTRACTED
peak20sec[1] = MeanOf20SecNOT_DISTRACTED

# BAR SCATTER FUNCTION = HERE BECAUSE DOESNT WORK IF IMPORTED FROM ALL FUNCS. will fix later
#
colors = ['thistle', 'gold']
colors2 = ['k','k']
colors3 = ['white', 'white']
def barscatter(data, transpose = False,
                groupwidth = .75,
                barwidth = .9,
                paired = False,
                barfacecoloroption = 'same', # other options 'between' or 'individual'
                barfacecolor = ['white'],
                baredgecoloroption = 'same',
                baredgecolor = ['black'],
                baralpha = 1,
                scatterfacecoloroption = 'same',
                scatterfacecolor = ['white'],
                scatteredgecoloroption = 'same',
                scatteredgecolor = ['grey'],
                scatterlinecolor = 'grey', # Don't put this value in a list
                scattersize = 80,
                scatteralpha = 1,
                linewidth=1,
                ylabel = 'none',
                xlabel = 'none',
                title = 'none',
                grouplabel = 'auto',
                itemlabel = 'none',
                yaxisparams = 'auto',
                show_legend = 'none',
                legendloc='upper right',
                ax=[]):
#
#    if type(data) == float
    # Check if transpose = True
    if transpose == True:
        data = np.transpose(data)
        
    # Initialize arrays and calculate number of groups, bars, items, and means
    
    barMeans = np.zeros((np.shape(data)))
    items = np.zeros((np.shape(data)))
    
    nGroups = np.shape(data)[0]
    groupx = np.arange(1,nGroups+1)

    if len(np.shape(data)) > 1:
        grouped = True
        barspergroup = np.shape(data)[1]
        barwidth = (barwidth * groupwidth) / barspergroup
        
        for i in range(np.shape(data)[0]):
            for j in range(np.shape(data)[1]):
                barMeans[i][j] = np.mean(data[i][j])
                items[i][j] = len(data[i][j])
        
    else:
        grouped = False
        paired = True
        barspergroup = 1
        
        for i in range(np.shape(data)[0]):
            barMeans[i] = np.mean(data[i])
            items[i] = len(data[i])
    
    # Calculate x values for bars and scatters    
    xvals = np.zeros((np.shape(data)))
    barallocation = groupwidth / barspergroup
    k = (groupwidth/2) - (barallocation/2)
    
    if grouped == True:
        
        for i in range(np.shape(data)[0]):
            xrange = np.linspace(i+1-k, i+1+k, barspergroup)
            for j in range(barspergroup):
                xvals[i][j] = xrange[j]
    else:
        xvals = groupx
    
# Set colors for bars and scatters  
    colors = ['thistle', 'gold']
    colors2 = ['k','k']
    colors3 = ['white', 'white']
    
    barfacecolorArray = setcolors("between", colors, 1, 2, data, paired_scatter = True)
    baredgecolorArray = setcolors("between", colors, 1, 2, data, paired_scatter = True)
     
    scfacecolorArray = setcolors("between", colors3, 1, 2, data, paired_scatter = True)
    scedgecolorArray = setcolors("between", colors2, 1, 2, data, paired_scatter = True)
 #   scfacecolorArray = setcolors("between", colors3, nGroups=nGroups, barspergroup=barspergroup, data=dataX, paired_scatter = True)
    
# Initialize figure
    if ax == []:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.tick_params(axis='both', which='major', labelsize=14)
 
    
    # Make bars
    barlist = []
    barx = []
    for x, y, bfc, bec in zip(xvals.flatten(), barMeans.flatten(),
                              barfacecolorArray, baredgecolorArray):
        barx.append(x)
        barlist.append(ax.bar(x, y, barwidth,
                         facecolor = bfc, edgecolor = bec,
                         zorder=-1))
    
    # Make scatters
    sclist = []
    if paired == False:
        for x, Yarray, scf, sce  in zip(xvals.flatten(), data.flatten(),
                                        scfacecolorArray, scedgecolorArray):
            for y in Yarray:
                sclist.append(ax.scatter(x, y, s = scattersize,
                         c = scf,
                         edgecolors = sce,
                         zorder=1))

    else:
        try:
            print('hey')
            np.shape(data)[1]
            for x, Yarray, scf, sce in zip(xvals, data, scfacecolorArray, scedgecolorArray):
                for y in np.transpose(Yarray.tolist()):
                    sclist.append(ax.plot(x, y, '-o', markersize = scattersize/10,
                             color = scatterlinecolor,
                             linewidth=linewidth,
                             markerfacecolor = scf,
                             markeredgecolor = sce))

# Explicitly added color here, issue with assignment of scf and sce 
        except IndexError:
                    
            print(len(data[0]))
            for n,_ in enumerate(data[0]):
                y = [y[n-1] for y in data]
                sclist.append(ax.plot(xvals, y, '-o', markersize = scattersize/10,
                             color = 'grey',
                             linewidth=linewidth,
                             markerfacecolor = 'white',
                             markeredgecolor = 'k'))
                

    # Label axes
    if ylabel != 'none':
        plt.ylabel(ylabel, fontsize=14)
    
    if xlabel != 'none':
        plt.xlabel(xlabel)
        
    if title != 'none':
        plt.title(title, fontsize=14)
    
    # Set range and tick values for Y axis
    if yaxisparams != 'auto':
        ax.set_ylim(yaxisparams[0])
        plt.yticks(yaxisparams[1])
       
    # X ticks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off') # labels along the bottom edge are off
    
    if grouplabel == 'auto':
        plt.tick_params(labelbottom='off')
    else:
        plt.xticks(range(1,nGroups+1), grouplabel)
    
    # Hide the right and top spines and set bottom to zero
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    
    
    if show_legend == 'within':
        if len(itemlabel) != barspergroup:
            print('Not enough item labels for legend!')
        else:
            legendbar = []
            legendtext = []
            for i in range(barspergroup):
                legendbar.append(barlist[i])
                legendtext.append(itemlabel[i])
            plt.legend(legendbar, legendtext, loc=legendloc)

    ax.set(ylabel='Peak ∆F')
    ax.yaxis.label.set_size(14)      
  #  fig.savefig('/Volumes/KPMSB352/PHOTOMETRY MMIN18/PDF figures/BarScat_20secAfter.pdf', bbox_inches="tight")        
    
    return ax, barx, barlist, sclist
      
def setcolors(coloroption, colors, barspergroup, nGroups, data, paired_scatter = False):
            
    nColors = len(colors)
    
    if (paired_scatter == True) & (coloroption == 'within'):
        print('Not possible to make a Paired scatter plot with Within setting.')
        coloroption = 'same'
        
    if coloroption == 'within':
        if nColors < barspergroup:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > barspergroup:
            colors = colors[:barspergroup]
        coloroutput = [colors for i in data]
        coloroutput = list(chain(*coloroutput))
        
    if coloroption == 'between':
        if nColors < nGroups:
            print('Not enough colors for this option! Reverting to one color.')
            coloroption = 'same'
        elif nColors > nGroups:
            colors = colors[:nGroups]
        if paired_scatter == False:
            coloroutput = [[c]*barspergroup for c in colors]
            coloroutput = list(chain(*coloroutput))
        else:
            coloroutput = colors
            
    if coloroption == 'individual':
        if nColors < nGroups*barspergroup:
            print('Not enough colors for this color option')
            coloroption = 'same'
        elif nColors > nGroups*barspergroup:
            coloroutput = colors[:nGroups*barspergroup]
        else: 
            coloroutput = colors
    
    if coloroption == 'same':
        coloroutput = [colors[0] for x in range(len(data.flatten()))]

    return coloroutput


colors = ['thistle', 'gold']
colors2 = ['k','k']
colors3 = ['white', 'white']
mpl.rc('lines', markeredgewidth=1)
ax = barscatter(PEAK2SEC, paired=True, scatterlinecolor='k', ylabel='Distracted trials (%)')
ax2 = barscatter(peak20sec, paired=True, scatterlinecolor='k', ylabel='Distracted trials (%)')


# Not sure about accuracy here but not significant. will run in R and cross ref (is the same)
from scipy import stats as stt
t2sec = stt.ttest_rel(MeanOf2SecNOT_DISTRACTED,MeanOf2SecDISTRACTED)
t20sec = stt.ttest_rel(MeanOf20SecNOT_DISTRACTED,MeanOf20SecDISTRACTED)
