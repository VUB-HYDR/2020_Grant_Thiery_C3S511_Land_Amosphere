#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:52:49 2020

@author: luke
"""


# =============================================================================
# SUMMARY
# =============================================================================


# This script runs correlation analyses between ERA5/ERA5-Land lake mixed layer  
# temperatures and evaporation or snowfall to identify 
# lake-atmosphere interactions.

# Options for analysis type can be adjusted under heading "OPTIONS - ANALYSIS".
    # Only adapt flags/thresholds with comment labels "<< SELECT >>" 

# Options for plotting can be adjusted under "OPTIONS - PLOTTING".


# =============================================================================
# IMPORT
# =============================================================================


import sys
import os


#==============================================================================
# PATH
#==============================================================================


# change current working directory to pathway of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# identify current directory for relative paths
curDIR = os.getcwd()

# data input directories
lmltDIR = os.path.join(curDIR, 'lakes/mixlayertemp/daily')
licdDIR = os.path.join(curDIR, 'lakes/icedepth')
mskDIR = os.path.join(curDIR, 'lakes/lakecover')
prDIR = os.path.join(curDIR, 'precipitation/daily')
etDIR = os.path.join(curDIR, 'evaporation_runoff/daily')
tasDIR = os.path.join(curDIR, 'temperature_pressure/daily')
uvDIR = os.path.join(curDIR, 'wind/monthly')


#==============================================================================
# OPTIONS - ANALYSIS
#==============================================================================


# adjust these flag settings for analysis choices only change "<< SELECT >>" lines

# << SELECT >>
flag_proc=1;      # 0; process and plot 
                  # 1; plot only (after processing has already been completed)

flag_svplt=1;     # 0: do not save plot
                  # 1: save plot in picDIR

# << SELECT >>
flag_ds=1;        # 0: ERA5
                  # 1: ERA5-Land

# << SELECT >>             
flag_coup=0;      # 0: process corr(lmlt,et) 
                  # 1: process corr(lmlt,tp) 
                  # 2: process corr(lmlt,tp) masked by corr(lmlt,et)
                  # 3: process lake effect snow; corr(lmlt,sf)
                  # 4: process lake effect snow; corr([lmlt-tas],sf)

# << SELECT >>
flag_corrtst=0;   # 0: pearson
                  # 1: spearman (rank)

# << SELECT >>
flag_subset=0;    # 0: do not subset part of globe
                  # 1: subset part of the globe

# << SELECT >> 
flag_wind=1;      # 0: do not overlay North America with wind vectors
                  # 1: zoom on great lakes; overlay with wind vectors

# << SELECT >>                                
cl_thrsh=0.15;    # minimum fraction of grid-scale lake cover to consider
et_thrsh=0.2;     # minimum corr(lmlt,et) to mask corr(lmlt,tp) given flag_coup == 2
ice_thrsh = 0.001 # ice depth threshold for masking open water (vals under this are ice-free)

# projections for global plotting
if flag_subset == 0:                
    
    if flag_coup == 0 or flag_coup == 1 or flag_coup == 2:    # lmlt -> et -> tp
        
        flag_proj=0;      # 0: robinson
                          
    elif flag_coup == 3 or flag_coup == 4:    # lake effect snow
        
        flag_proj=1       # 1: azimuthal equidistant

# projections for regional plotting
elif flag_subset == 1:  # choose region
    
    # << SELECT >>
    flag_proj=2;          # 2: lambert conformal (Canada/US)
                          # 3: platecaree (Africa)
    
# string for ERA5 vs ERA5-Land supplementary file reading
if flag_ds == 0:
    
    ds = 'era5'
    y1 = '1979'

elif flag_ds == 1:
    
    ds = 'era5-land'
    y1 = '1981'

variables = ['lmlt',    # lake mixed layer temp
             'licd',    # lake ice cover depth
             'sf',      # snowfall
             'tp',      # total precipitation
             'et',      # evapotranspiration
             'tas']     # 2m temperature
correlation_tests = ['pearson', 
                     'spearman']
projections = ['Robinson',
               'AzimuthalEquidistant',
               'LambertConformal',
               'PlateCarree']

corr_test = correlation_tests[flag_corrtst]
projec = projections[flag_proj]

# lake cover and ice cover masks
msk_file = ds + '_lakes_cl_mask.nc'
licd_oct_msk_file = ds + '_lakes_icedepth_timmean_OCT_' + y1 + '_2019.nc'
licd_nov_msk_file = ds + '_lakes_icedepth_timmean_NOV_' + y1 + '_2019.nc'
licd_dec_msk_file = ds + '_lakes_icedepth_timmean_DEC_' + y1 + '_2019.nc'

# fall wind vector maps
uv_oct_vec = 'oct_' + ds + '_u_v_10m_' + y1 + '_2019.nc'
uv_nov_vec = 'nov_' + ds + '_u_v_10m_' + y1 + '_2019.nc'
uv_dec_vec = 'dec_' + ds + '_u_v_10m_' + y1 + '_2019.nc'
uv_ond_vec = 'ond_' + ds + '_u_v_10m_' + y1 + '_2019.nc'

# output for image
if flag_subset == 0:
    
    picDIR = '/theia/scratch/projects/climate/users/lgrant/' + ds + '/coupling/globe'
    ncDIR = '/theia/scratch/projects/climate/users/lgrant/' + ds + '/coupling/output_nc'
    
elif flag_subset == 1:
    
    picDIR = '/theia/scratch/projects/climate/users/lgrant/' + ds + '/coupling/subset'
    ncDIR = '/theia/scratch/projects/climate/users/lgrant/' + ds + '/coupling/output_nc'

# input directory for plotting only
if flag_proc == 1:
    
    corrDIR = '/theia/scratch/projects/climate/users/lgrant/' + ds + '/coupling/output_nc'


#==============================================================================
# OPTIONS - PLOTTING
#==============================================================================


# font settings << SELECT >>   
title_font = 12
cbtitle_font = 12
tick_font = 9

# colorbar features << SELECT >>
col_cbticlbl = '0'   # colorbar color of tick labels
col_cbtic = '0.5'   # colorbar color of ticks
col_cbedg = '0.9'   # colorbar color of edge
cb_ticlen = 3.5   # colorbar length of ticks
cb_ticwid = 0.4   # colorbar thickness of ticks
cb_edgthic = 0   # colorbar thickness of edges between colors
cblabel = 'corr'  # colorbar label
col_zero = 'gray'   # zero change color
sbplt_lw = 0.1   # linewidth on projection panels
cstlin_lw = 0.2   # linewidth for coastlines

# figsize, panel spacing, colorbar loc, plot boundaries << SELECT >>
if projec == 'Robinson':
    
    x = 12   # width of fig (inches)
    y = 8   # height of fig (inches)
    pad = 10   # distance from panel to plot border
    w_pad = 5   # horiz distance between panels
    h_pad = 5   # vert distance between panels
    cb_x0 = 0.235   # x starting point
    cb_y0 = 0.10   # y starting point
    cb_xlen = 0.55   # x length
    cb_ylen = 0.015   # y length
    east = 50 # ignore east, west, north and south for this projection
    west = -15
    north = 30
    south = -30
    
elif projec == 'AzimuthalEquidistant':
    
    x = 12
    y = 10
    pad = 0
    w_pad = 5
    h_pad = 2.5
    cb_x0 = 0.235   
    cb_y0 = 0.05   
    cb_xlen = 0.55   
    cb_ylen = 0.015
    east = 50 # ignore east, west, north and south for this projection
    west = -15
    north = 30
    south = -30
    
elif projec == 'LambertConformal':
    
    x = 12
    y = 10
    pad = 0
    w_pad = 5
    h_pad = 2.5
    cb_x0 = 0.235   
    cb_y0 = 0.05   
    cb_xlen = 0.55   
    cb_ylen = 0.015
    
    if flag_wind == 0:
        
        east = -63
        west = -123
        north = 75
        south = 37
        
    elif flag_wind == 1:
        
        east = -75
        west = -93
        north = 49
        south = 41
        

elif projec == 'PlateCarree':
    
    x = 10
    y = 12
    pad = 0
    w_pad = 5
    h_pad = 2.5
    cb_x0 = 0.235   
    cb_y0 = 0.05   
    cb_xlen = 0.55   
    cb_ylen = 0.015
    east = 50
    west = -15
    north = 30
    south = -30


#==============================================================================
# proc LMLT -> ET -> corr(lmlt,et)
#==============================================================================


if flag_coup == 0:
    
    lak_var = variables[0]
    pr_var = variables[4]
    
    from funcs import *
    
    if flag_proc == 0:

        # lists of lmlt and et files
        lak_files = []
        et_files = []
        
        # grab available files
        for file in [file for file in sorted(os.listdir(lmltDIR))\
                     if lak_var in file]:
            lak_files.append(file)
        for file in [file for file in sorted(os.listdir(etDIR))\
                     if pr_var in file]:
            et_files.append(file)
    
        # pick warmest 3 months
        os.chdir(lmltDIR)
        monthly_lmlt,weights = months_mean(lak_files)
        months = warm_months(monthly_lmlt, weights)
        
        # retrieve maps of correlations
        corr_maps_avg,labels,corr_maps = corr(lmltDIR,etDIR,
                                lak_files,et_files,
                                months,corr_test)
        
    elif flag_proc == 1:
        
        labels = {}
        labels[7] = 'July'
        labels[8] = 'August'
        labels[9] = 'September'
        labels['all'] = 'July - August - September'
        
        corr_maps_avg = {}
        os.chdir(corrDIR)
        for label in labels.values():
            # file; eg couple_era5-land_4_December.nc
            if label != 'July - August - September':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_' + label + '.nc')
            elif label == 'July - August - September':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_July_August_September.nc')
        
    
#==============================================================================
# proc LMLT -> RAINFALL -> corr(lmlt,tp)
#==============================================================================


if flag_coup == 1:
    
    lak_var = variables[0]
    pr_var = variables[3]
    
    from funcs import *
    
    if flag_proc == 0:

        # lists of lmlt and tp files
        lak_files = []
        tp_files = []
        
        # grab available files
        for file in [file for file in sorted(os.listdir(lmltDIR))\
                     if lak_var in file]:
            lak_files.append(file)
        for file in [file for file in sorted(os.listdir(prDIR))\
                     if pr_var in file]:
            tp_files.append(file)
            
        # pick warmest 3 months
        os.chdir(lmltDIR)
        monthly_lmlt,weights = months_mean(lak_files)
        months = warm_months(monthly_lmlt, weights)
        
        # retrieve maps of correlations
        corr_maps_avg,labels,corr_maps = corr(lmltDIR,prDIR,
                                lak_files,tp_files,
                                months,corr_test)
        
    elif flag_proc == 1:
        
        labels = {}
        labels[7] = 'July'
        labels[8] = 'August'
        labels[9] = 'September'
        labels['all'] = 'July - August - September'
        
        corr_maps_avg = {}
        os.chdir(corrDIR)
        for label in labels.values():
            # file; eg couple_era5-land_4_December.nc
            if label != 'July - August - September':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_' + label + '.nc')
            elif label == 'July - August - September':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_July_August_September.nc')
        
    
#==============================================================================
# proc LMLT -> RAINFALL where corr(LMLT,ET) > et_thrsh
#==============================================================================


if flag_coup == 2:
    
    lak_var = variables[0]
    pr_var = variables[4]
    
    from funcs import *
    
    if flag_proc == 0:

        # lists of lmlt and et files
        lak_files = []
        et_files = []
        
        # grab available files
        for file in [file for file in sorted(os.listdir(lmltDIR))\
                     if lak_var in file]:
            lak_files.append(file)
        for file in [file for file in sorted(os.listdir(etDIR))\
                     if pr_var in file]:
            et_files.append(file)
    
        # pick warmest 3 months
        os.chdir(lmltDIR)
        monthly_lmlt,weights = months_mean(lak_files)
        months = warm_months(monthly_lmlt, weights)
        
        # retrieve maps of correlations for lmlt <-> et
        corr_maps_avg_et,labels,corr_maps_et = corr(lmltDIR,etDIR,
                                lak_files,et_files,
                                months,corr_test)
        
        # redefine pr_var for corr(lmlt,tp)
        pr_var = variables[3]
    
        # list for tp files
        tp_files = []
        
        # grab available files
        for file in [file for file in sorted(os.listdir(prDIR))\
                     if pr_var in file]:
            tp_files.append(file)
            
        # fix months
        months = months[:-1]
        
        # retrieve maps of correlations
        corr_maps_avg,labels,corr_maps = corr(lmltDIR,prDIR,
                                lak_files,tp_files,
                                months,corr_test)
        
        # isolate to regions with "strong" lmlt -> et coupling
        for lab in labels.values():
            corr_maps_avg[lab] = corr_maps_avg[lab].where(corr_maps_avg_et[lab]>et_thrsh)
            corr_maps[lab] = corr_maps[lab].where(corr_maps_et[lab]>et_thrsh)
            
    
    elif flag_proc == 1:
        
        labels = {}
        labels[7] = 'July'
        labels[8] = 'August'
        labels[9] = 'September'
        labels['all'] = 'July - August - September'
        
        corr_maps_avg = {}
        os.chdir(corrDIR)
        for label in labels.values():
            # file; eg couple_era5-land_4_December.nc
            if label != 'July - August - September':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_' + label + '.nc')
            elif label == 'July - August - September':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_July_August_September.nc')


#==============================================================================
# proc LAKE EFFECT SNOW -> corr(lmlt,sf)
#==============================================================================
    
    
elif flag_coup == 3:

    lak_var = variables[0]
    pr_var = variables[2]
    
    from funcs import *
    
    if flag_proc == 0:

        # lists of lmlt and sf files
        lak_files = []
        sf_files = []
        
        
        # grab available files
        for file in [file for file in sorted(os.listdir(lmltDIR))\
                     if lak_var in file]:
            lak_files.append(file)
        for file in [file for file in sorted(os.listdir(prDIR))\
                     if pr_var in file]:
            sf_files.append(file)

        # October, November and December chosen based on icestart timmean plot from spqb
        months = [10,11,12]
        
        # retrieve maps of correlations
        corr_maps_avg,labels,corr_maps = corr(lmltDIR,prDIR,
                                lak_files,sf_files,
                                months,corr_test)
    
    elif flag_proc == 1:
        
        labels = {}
        labels[10] = 'October'
        labels[11] = 'November'
        labels[12] = 'December'
        labels['all'] = 'October - November - December'
        
        corr_maps_avg = {}
        os.chdir(corrDIR)
        for label in labels.values():
            # file; eg couple_era5-land_4_December.nc
            if label != 'October - November - December':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_' + label + '.nc')
            elif label == 'October - November - December':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_October_November_December.nc')


#==============================================================================
# proc LAKE EFFECT SNOW -> corr([lmlt-tas],sf)
#==============================================================================
 
    
elif flag_coup == 4:
    
    lak_var = variables[0]
    pr_var = variables[2]
    tas_var = variables[-1]
    
    from funcs import *
    
    if flag_proc == 0:
    
        # lists of lmlt, sf and tas files
        lak_files = []
        sf_files = []
        tas_files = []
        
        # grab available files
        for file in [file for file in sorted(os.listdir(lmltDIR))\
                     if lak_var in file]:
            lak_files.append(file)
        for file in [file for file in sorted(os.listdir(prDIR))\
                     if pr_var in file]:
            sf_files.append(file)
        for file in [file for file in sorted(os.listdir(tasDIR))\
                     if tas_var in file]:
            tas_files.append(file)
        
        # October, November and December chosen based on icestart timmean plot from spqb
        months = [10,11,12]
        
        # retrieve maps of correlations
        corr_maps_avg,labels,corr_maps = corr_les2(lmltDIR,prDIR,tasDIR,
                                lak_files,sf_files,tas_files,
                                months,corr_test)
    
    
    elif flag_proc == 1:
        
        labels = {}
        labels[10] = 'October'
        labels[11] = 'November'
        labels[12] = 'December'
        labels['all'] = 'October - November - December'
        
        corr_maps_avg = {}
        os.chdir(corrDIR)
        for label in labels.values():
            # file; eg couple_era5-land_4_December.nc
            if label != 'October - November - December':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_' + label + '.nc')
            elif label == 'October - November - December':
                corr_maps_avg[label] = reader('couple_' + ds + '_' + str(flag_coup) + '_October_November_December.nc')
    
    
#==============================================================================
# MASK + PLOT RESULTS
#==============================================================================   


# wind vectors
os.chdir(uvDIR)
uv = {}

if flag_subset == 1 and flag_wind == 1 and projec == 'LambertConformal':
    for lab in labels.values():
        if lab == 'October':
            uv[lab] = reader(uv_oct_vec)
        elif lab == 'November':
            uv[lab] = reader(uv_nov_vec)
        elif lab == 'December':
            uv[lab] = reader(uv_dec_vec)
        elif lab == 'October - November - December':
            uv[lab] = reader(uv_ond_vec)

# mask data based on lake area and ice presence
corr_maps_avg = masker(corr_maps_avg,
                       labels,
                       pr_var,
                       mskDIR,
                       msk_file,
                       cl_thrsh,
                       ice_thrsh,
                       licdDIR,
                       licd_oct_msk_file,
                       licd_nov_msk_file,
                       licd_dec_msk_file)

# plot maps of correlations
corr_plot = plot_corr(flag_svplt,ds,picDIR,corr_maps_avg,labels,
                      x,
                      y,
                      east,
                      west,
                      north,
                      south,
                      projec,
                      title_font,
                      cbtitle_font,
                      tick_font,
                      col_cbticlbl,
                      col_cbtic,
                      col_cbedg,
                      cb_ticlen,
                      cb_ticwid,
                      cb_edgthic,
                      cb_x0,
                      cb_y0,
                      cb_xlen,
                      cb_ylen,
                      cblabel,
                      col_zero,
                      sbplt_lw,
                      cstlin_lw,
                      pad,
                      w_pad,
                      h_pad,
                      lak_var,
                      pr_var,
                      flag_coup,
                      flag_wind,
                      cl_thrsh,
                      uv)


#==============================================================================
# SAVE TO NC
#==============================================================================


if flag_proc == 0:
    
    for lab in labels.values():
        if not '-' in lab:
            corr_maps[lab].to_netcdf(path=ncDIR + '/couple_' + ds +'_' + str(flag_coup) + '_' + lab + '.nc')
        elif '-' in lab:
            lab = lab.replace(" - ", "_")
            corr_maps[lab].to_netcdf(path=ncDIR + '/couple_' + ds +'_' + str(flag_coup) + '_' + lab + '.nc')
        
elif flag_proc == 1:
    
    None
    
        

