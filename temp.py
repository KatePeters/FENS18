#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 09:29:50 2018

@author: u1490431
"""

import scipy.io as sio


TDTfileslist = ['NAPH02_distraction', 'NAPH03_distraction', 'NAPH04_distraction',\
                'NAPH05_distraction', 'NAPH07_distraction', 'NAPH08_distraction',\
                'NAPH09_distraction', 'NAPH10_distraction']

#TDTfileslist = ['NAPH09_distraction']
TDTfilepath = '/Volumes/KP_HARD_DRI/kp259/NAPH1/MATLAB/'

#file = '/Volumes/KP_HARD_DRI/kp259/NAPH1/MATLAB/NAPH07_licktrain.mat'

for x in TDTfileslist:
    a = sio.loadmat(TDTfilepath + x, squeeze_me=True, struct_as_record=False)
    print(type(a))
    
    
    
# LiA_   La2_    


def mapTTLs(matdict):
    for x in ['LiA_', 'La2_']:
        try:
            licks = getattr(matdict['output'], x)
        except:
            print('File has no ' + x)
            
    lickson = licks.onset
    licksoff = licks.offset 
            
    return lickson, licksoff 

testlickson, testlicksoff = mapTTLs(a)
        
# if want distractionr add argument of list of ttls (x) make more general - licks - ttls 


        
#getattr(b,'LiA_')