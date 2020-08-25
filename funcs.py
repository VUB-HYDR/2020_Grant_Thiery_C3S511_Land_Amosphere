#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 08:37:09 2020

@author: luke
"""

# =============================================================================
# IMPORT
# =============================================================================

import os
import sys
import xarray as xr
import numpy as np
from scipy import stats as sts
from scipy import signal as sgl
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.stats import spearmanr
import bottleneck


# =============================================================================
# FUNCTIONS
# =============================================================================

def reader(file):
    
    """ Read in netcdfs based on variable.
    
    Parameters
    ----------
    file : files in data directory
    
    Returns
    ------- 
    Xarray data array
    """
    
    if not 'icedepth_timmean' in file:
        ds = xr.open_dataset(file)
        ds = ds.squeeze(drop=True)
        
    elif 'icedepth_timmean' in file:
        ds = xr.open_dataset(file, decode_times=False)
        ds = ds.squeeze(drop=True)
        
    if 'lmlt' in file:
        da = ds.lmlt
        
    elif 'licd' in file or 'icedepth_timmean' in file:
        da = ds.licd
        
    elif 'tp' in file:
        da = ds.tp
        
    elif 'sf' in file:
        da = ds.sf
        
    elif 'cl' in file:
        da = ds.cl
        
    elif 'era5_et' in file:
        da = ds.e*-1   # invert values for evap as (+) and condensation as (-)
    
    elif 'era5-land_et' in file:
        da = ds.evaow*-1   # invert values for evap as (+) and condensation as (-)
        
    elif 'tas' in file:
        da = ds.t2m
        
    elif '_u_v_' in file:
        da = ds
        
    elif 'couple' in file:
        da = ds.__xarray_dataarray_variable__
        
    return da

def date_range(concat_dim,m):
        
    """"Give time var for correlation map output in nc files.
    
    Parameters
    ---------- 
    concat_dim : np array of len = len(number corr maps)
    m : integer for month
    
    Returns
    -------
    xarray var for time
    """
    
    y1 = 1979
    y2 = y1+concat_dim[-1]
    periods = y2 - y1 + 1
    y1 = str(y1)
    y2 = str(y2)
    m1 = m
    if len(str(m1)) == 1:
        m1 = f'0{str(m1)}'
    elif len(str(m1)) == 2:
        m1 = str(m1)
    d1 = '15'
    start = f'{y1}-{m1}-{d1}'
    end = f'{y2}-{m1}-{d1}'
    date_range = pd.date_range(start=start,end=end,periods=periods)
    time = xr.Variable('time', date_range)
    
    return time
    
def months_mean(lak_files):
    
    """"Convert to monthly files and grab weights.
    
    Parameters
    ---------- 
    lake_files : all available lmlt files
    
    Returns
    -------
    monthly concatenated series for all available daily files
    """
    
    monthly_arrays = []
    for l in lak_files:
        
        daily_data = reader(l)
        monthly_data = daily_data.resample(time="1MS").mean(dim="time")
        monthly_arrays.append(monthly_data)
        
    monthly_series = xr.concat(monthly_arrays,dim="time")
    weights = np.cos(np.deg2rad(monthly_series.latitude))
    return monthly_series,weights

def weighted_mean(da, weights, dim):
    """ Reduce da by a weighted mean along some dimension(s).
    
    Parameters
    ----------
    da : DataArray
        Object over which the weighted reduction operation is applied.
    weights : DataArray
        An array of weights associated with the values in this Dataset.
    dim : str or sequence of str, optional
        Dimension(s) over which to apply the weighted `mean`.

    Returns
    -------
    weighted_mean : DataArray
        New DataArray with weighted mean applied to its data and
        the indicated dimension(s) removed.
    """
    
    weighted_sum = (da * weights).sum(dim=dim, skipna=True)
    masked_weights = weights.where(da.notnull())
    sum_of_weights = masked_weights.sum(dim=dim, skipna=True)
    valid_weights = sum_of_weights != 0
    sum_of_weights = sum_of_weights.where(valid_weights)
    
    return weighted_sum / sum_of_weights

def warm_months(monthly_series,weights):
    
    """ Find warmest three months across monthly dataset with weighted global mean.
    
    Parameters
    ----------
    monthly_series : result of "months" function
        Object over which the weighted reduction operation is applied.
    weights : DataArray (also result of "months" function)
        An array of weights associated with the values in this Dataset.

    Returns
    -------
    Mode of top 3 warmest months
    """    

    glob_means = weighted_mean(monthly_series,weights,('latitude', 'longitude'))
    
    years = monthly_series.resample(time="YS").mean(dim="time").time.values
    years = pd.to_datetime(years)
    
    hot_dict = {}
    hot_list = []
    
    for y in years:
        
        glob_means_y = glob_means.sel(time=str(y.year))
    
        max1 = glob_means_y.where(glob_means_y==glob_means_y.max(), drop=True).squeeze()
    
        glob_means_y2 = glob_means_y.where(glob_means_y.time.values != max1.time.values)
        max2 = glob_means_y2.where(glob_means_y2==glob_means_y2.max(skipna=True), drop=True).squeeze()
    
        glob_means_y3 = glob_means_y2.where(glob_means_y2.time.values != max2.time.values)
        max3 = glob_means_y3.where(glob_means_y3==glob_means_y3.max(skipna=True), drop=True).squeeze()
        
        max1_month = pd.to_datetime(max1.time.values).month
        max2_month = pd.to_datetime(max2.time.values).month
        max3_month = pd.to_datetime(max3.time.values).month
        
        hot_dict[str(y.year)] = [max1_month, max2_month, max3_month]
        hot_list.append(hot_dict[str(y.year)])
        
    hot_months = np.vstack(hot_list)
    m1 = sts.mode(hot_months[:,0])[0][0]
    m2 = sts.mode(hot_months[:,1])[0][0]
    m3 = sts.mode(hot_months[:,2])[0][0]
        
    return sorted([m1,m2,m3])

def covariance_gufunc(x, y):
    return ((x - x.mean(axis=-1, keepdims=True))
            * (y - y.mean(axis=-1, keepdims=True))).mean(axis=-1)

def pearson_correlation_gufunc(x, y):
    return covariance_gufunc(x, y) / (x.std(axis=-1) * y.std(axis=-1))

def pearson_correlation(x,y,dim):
    return xr.apply_ufunc(
        pearson_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        output_dtypes=[float])

def spearman_correlation_gufunc(x, y):
    x_ranks = bottleneck.rankdata(x, axis=-1)
    y_ranks = bottleneck.rankdata(y, axis=-1)
    return pearson_correlation_gufunc(x_ranks, y_ranks)

def spearman_correlation(x, y, dim):
    return xr.apply_ufunc(
        spearman_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        output_dtypes=[float])

def month_labels(months):
    
    """ Retrieve month labels for keying corr maps and plotting .
    
    Parameters
    ----------
    months : result of warm_months func

    Returns
    -------
    4-keyed dictionary with labels for 3 individual months and 'month1 - month2 - month3'
    """    
    
    month_dict = {1: 'January',
                  2: 'February',
                  3: 'March',
                  4: 'April',
                  5: 'May',
                  6: 'June',
                  7: 'July',
                  8: 'August',
                  9: 'September',
                  10: 'October',
                  11: 'November',
                  12: 'December'}
    labels = {}
    
    for month in months:
        
        labels[month] = month_dict[month]
        
    labels['all'] = month_dict[months[0]] + ' - ' +month_dict[months[1]] + ' - ' + month_dict[months[2]]
    return labels

def corr(lmltDIR,prDIR,
         lak_files,tp_files,
         months,corr_test):
    
    """ Take corr(lmlt,tp) from given three months of lake lmlt.
    Can work on lake effect snow corr(lmlt,sf) or any corr(lmlt,var) proc
    
    Parameters
    ----------
    files : list of annual files at daily t-step and individual mask file
    months : result of warm_months func
    corr_test : pearson or spearman (rank) correlation
    
        pearson uses funcs; "covariance_gufunc" and "pearson_correlation_gufunc"
    
        spearman uses funcs; "covariance_gufunc","pearson_correlation_gufunc",
            "spearman_correlation_gufunc" and "spearman_correlation"

    Returns
    -------
    4 corr maps as avg for three months and season
    """      
    
    labels = month_labels(months)

    corr_maps_list = {}
    corr_maps_avg = {}
    
    for m in months:
        
        corr_maps_list[labels[m]] = []
    
    corr_maps_list[labels['all']] = []
    
    
    for l,t in zip(lak_files,tp_files):
        
        os.chdir(lmltDIR)
        da_l = reader(l)
        
        os.chdir(prDIR)
        da_t = reader(t)
        
        corr_maps = {}
        
        count = 0
        # select months, correlate
        for m in months:

            count += 1
            da_l_m = da_l.sel(time=np.in1d(da_l['time.month'],m))
            da_t_m = da_t.sel(time=np.in1d(da_t['time.month'],m)).where(~np.isnan(da_l_m.values))
            
            if count == 3:
                da_l_all = da_l = da_l.sel(time=np.in1d(da_l['time.month'],months))
                da_t_all = da_t = da_t.sel(time=np.in1d(da_t['time.month'],months))
            
            # setup for corr
            if corr_test == 'pearson':
                corr_maps[labels[m]] = pearson_correlation(da_l_m,da_t_m,dim='time')
        
                if count == 3:
                    corr_maps[labels['all']] = pearson_correlation(da_l_all,da_t_all,dim='time')
                
            elif corr_test == 'spearman':
                corr_maps[labels[m]] = spearman_correlation(da_l_m,da_t_m,dim='time')
                
                if count == 3:
                    corr_maps[labels['all']] = spearman_correlation(da_l_all,da_t_all,dim='time')
                    
            
            corr_maps_list[labels[m]].append(corr_maps[labels[m]])
            if count == 3:
                corr_maps_list[labels['all']].append(corr_maps[labels['all']])
    
    months.append(months[1])
    for label,m in zip(labels.values(),months):
        
        concat_dim = np.arange(len(corr_maps_list[label]))
        time = date_range(concat_dim,m)
        corr_maps[label] = xr.concat(corr_maps_list[label],dim=time)
        
        corr_maps_avg[label] = xr.concat(corr_maps_list[label],dim=concat_dim).mean(dim='concat_dim')
        
    return corr_maps_avg,labels,corr_maps

def corr_les2(lmltDIR,prDIR,tasDIR,
                lak_files,tp_files,tas_files,
                months,corr_test):
    
    """ Take corr([lmlt-tas,sf) from warmest three months of lake lmlt or latest ice free months.
    "les2" for lake effect snow proc v2 for looking to lmlt-tas
    
    Parameters
    ----------
    files : list of annual files at daily t-step and individual mask file
    months : result of warm_months func or chosen winter months
    corr_test : pearson or spearman (rank) correlation
    
        pearson uses funcs; "covariance_gufunc" and "pearson_correlation_gufunc"
    
        spearman uses funcs; "covariance_gufunc","pearson_correlation_gufunc",
            "spearman_correlation_gufunc" and "spearman_correlation"

    Returns
    -------
    4 corr maps as avg for three months and season
    """      
    
    labels = month_labels(months)

    corr_maps_list = {}
    corr_maps_avg = {}
    
    for m in months:
        
        corr_maps_list[labels[m]] = []
    
    corr_maps_list[labels['all']] = []
    
    
    for l,t,tas in zip(lak_files,tp_files,tas_files):
        
        os.chdir(lmltDIR)
        da_l = reader(l)
        
        os.chdir(prDIR)
        da_t = reader(t)
        
        os.chdir(tasDIR)
        da_tas = reader(tas)
        
        # take difference between air and water temp for correlation
        da_l = da_l - da_tas
        
        corr_maps = {}
        
        count = 0
        # select months, detrend data arrays, correlate
        for m in months:

            count += 1
            da_l_m = da_l.sel(time=np.in1d(da_l['time.month'],m))
            da_t_m = da_t.sel(time=np.in1d(da_t['time.month'],m)).where(~np.isnan(da_l_m.values))
            
            if count == 3:
                da_l_all = da_l = da_l.sel(time=np.in1d(da_l['time.month'],months))
                da_t_all = da_t = da_t.sel(time=np.in1d(da_t['time.month'],months))
            
            # setup for corr
            if corr_test == 'pearson':
                corr_maps[labels[m]] = pearson_correlation(da_l_m,da_t_m,dim='time')
        
                if count == 3:
                    corr_maps[labels['all']] = pearson_correlation(da_l_all,da_t_all,dim='time')
                
            elif corr_test == 'spearman':
                corr_maps[labels[m]] = spearman_correlation(da_l_m,da_t_m,dim='time')
                
                if count == 3:
                    corr_maps[labels['all']] = spearman_correlation(da_l_all,da_t_all,dim='time')
                    
            
            corr_maps_list[labels[m]].append(corr_maps[labels[m]])
            if count == 3:
                corr_maps_list[labels['all']].append(corr_maps[labels['all']])
    
    months.append(months[1])
    for label,m in zip(labels.values(),months):
        
        concat_dim = np.arange(len(corr_maps_list[label]))
        time = date_range(concat_dim,m)
        corr_maps[label] = xr.concat(corr_maps_list[label],dim=time)
        
        corr_maps_avg[label] = xr.concat(corr_maps_list[label],dim=concat_dim).mean(dim='concat_dim')
        
    return corr_maps_avg,labels,corr_maps

def masker(corr_maps_avg,
           labels,
           pr_var,
           mskDIR,
           msk_file,
           cl_thrsh,
           ice_thrsh,
           licdDIR,
           licd_oct_msk_file,
           licd_nov_msk_file,
           licd_dec_msk_file):
    
    """ First masking all correlation maps for lakes covering a minimum grid cell fraction.
    Then, mask correlation maps for lake effect snow by lakes which do not have ice cover (defined
    by a minimum threshold for ice depth)
    
    Parameters
    ----------
    corr_maps_avg : output from correlation processing
    labels : months from analysis
    pr_var : snowfall (sf) or total precipitation (tp)
    mskDIR : directory with lake cover mask
    cl_thrsh : gridscale fraction of lake cover threshold for consideration
    ice_thrsh : threshold for ice depth under which lakes are considered ice free for 
        the lake effect snow processing.
    licdDIR : directory with ice depth files used for masking months in lake effect snow correlation

    Returns
    -------
    4 corr maps as avg for three months and season
    """     
    
    os.chdir(mskDIR)
    msk = reader(msk_file)

    for lab in labels.values():
            
        corr_maps_avg[lab] = corr_maps_avg[lab].where(msk > cl_thrsh)
        
        if pr_var == 'sf':
            
            os.chdir(licdDIR)
            if lab == 'October':
                licd_msk = reader(licd_oct_msk_file)
                
            elif lab == 'November':
                licd_msk = reader(licd_nov_msk_file)
                
            elif lab == 'December' or lab == 'October - November - December':
                licd_msk = reader(licd_dec_msk_file)                
            
            licd_msk_vals = licd_msk.values
            licd_msk_align = xr.DataArray(licd_msk_vals,coords=msk.coords,dims=msk.dims)
            corr_maps_avg[lab] = corr_maps_avg[lab].where(licd_msk_align < ice_thrsh)
    
    return corr_maps_avg

def projector(projec,east,west,north,south):
    
    """ Provide extent details based on projection type.
    
    Parameters
    ----------
    All defined on main.py

    Returns
    -------
    'new_extent' under given projection and object for projection ('proj')
    """    
    
        # defining extents based on projection
    if projec == 'Robinson':
        proj = ccrs.Robinson()
        new_extent = 'filler'
        
    elif projec == 'AzimuthalEquidistant':
        radius=6371228
        ex = 9024309
        adj = 4000000
        cornerx=ex+adj 
        newcornerx=cornerx/2.
        newcornery=newcornerx
        globe = ccrs.Globe(ellipse=None, semimajor_axis=radius, 
                   semiminor_axis=radius)
        proj = ccrs.AzimuthalEquidistant(central_latitude=90,
                                          central_longitude=-30,
                                          globe=globe)
        new_extent = [-newcornerx*1.25,newcornerx*1.25,-newcornery,newcornery]
        
    elif projec == 'LambertConformal':
        canada_east = east
        canada_west = west
        canada_north = north
        canada_south = south
        
        standard_parallels = (49, 77)
        central_longitude = -(91 + 52 / 60)
        
        proj = ccrs.LambertConformal(central_longitude=central_longitude,
                                     standard_parallels=standard_parallels)
        new_extent = [canada_west,canada_east,canada_south,canada_north]
        
    elif projec == 'PlateCarree':
        africa_east = east
        africa_west = west
        africa_north = north
        africa_south = south
        
        proj = ccrs.PlateCarree()
        new_extent = [africa_west,africa_east,africa_south,africa_north]
        
    return proj,new_extent
            
def plot_corr(flag_svplt,ds,picDIR,corr_maps_avg,labels,
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
              uv):
    
    """ Plot given correlation maps.
    
    Parameters
    ----------
    All defined on main.py

    Returns
    -------
    4 corr maps in .png
    """    
    
    letters = ['a','b','c','d']
    
    # identify colors
    cmap_whole = plt.cm.get_cmap('RdBu_r')
    cmap55 = cmap_whole(0.01)
    cmap50 = cmap_whole(0.05)   #blue
    cmap45 = cmap_whole(0.1)
    cmap40 = cmap_whole(0.15)
    cmap35 = cmap_whole(0.2)
    cmap30 = cmap_whole(0.25)
    cmap25 = cmap_whole(0.3)
    cmap20 = cmap_whole(0.325)
    cmap10 = cmap_whole(0.4)
    cmap5 = cmap_whole(0.475)
    cmap0 = col_zero
    cmap_5 = cmap_whole(0.525)
    cmap_10 = cmap_whole(0.6)
    cmap_20 = cmap_whole(0.625)
    cmap_25 = cmap_whole(0.7)
    cmap_30 = cmap_whole(0.75)
    cmap_35 = cmap_whole(0.8)
    cmap_40 = cmap_whole(0.85)
    cmap_45 = cmap_whole(0.9)
    cmap_50 = cmap_whole(0.95)  #red
    cmap_55 = cmap_whole(0.99)
    
    colors = [cmap_45,cmap_35,cmap_25,cmap_10,
              cmap0,
              cmap10,cmap25,cmap35,cmap45]
    
    # declare list of colors for discrete colormap of colorbar
    cmap = mpl.colors.ListedColormap(colors,N=len(colors))

    # colorbar args
    values = [-1,-.75,-.50,-.25,-0.01,0.01,.25,.5,.75,1]
    tick_locs = [-1,-.75,-.50,-.25,0,.25,.5,.75,1]
    tick_labels = ['-1','-0.75','-0.50','-0.25','0','0.25','0.50','0.75','1']
    norm = mpl.colors.BoundaryNorm(values,cmap.N)

    # defining extents based on projection
    proj,new_extent = projector(projec,east,west,north,south)
    
    # initiate plotting
    f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,
                   figsize=(x,y),
                   subplot_kw=dict(projection=proj))
    
    # add colorbar axis
    cbax = f.add_axes([cb_x0, cb_y0, cb_xlen, cb_ylen])
    axes = [ax1,ax2,ax3,ax4]
    
    # loop through paired axes/labels
    count = 0
    for ax,lab in zip(axes,labels.values()):
        count += 1
        corr_maps_avg[lab].plot.pcolormesh(
            ax=axes[count-1], 
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            cbar_ax=cbax,
            center=0,
            norm=norm,
            vmin=-1,
            vmax=1,
            add_labels=False)
        ax.set_title(lab,fontsize=title_font,loc='center',pad=5)
        ax.set_title(letters[count-1],fontsize=title_font,fontweight='bold',loc='left')
        ax.coastlines(linewidth=cstlin_lw);
        
        if projec == 'LambertConformal':
            # add lakes at higher resolution for Canada
            lake_50m = cfeature.NaturalEarthFeature('physical','lakes','50m',
                                                    edgecolor='k',
                                                    linewidth=0.5,
                                                    facecolor="none",
                                                    alpha=1)
            ax.add_feature(lake_50m)
            
            if flag_wind == 1:
                # plot wind vectors
                skip_vec = (slice(None,None,6),slice(None,None,6))
                skip_latlon = slice(None,None,6)
                u = uv[lab].u10
                v = uv[lab].v10
                uv_lat = u.latitude.values
                uv_lon = u.longitude.values
                ax.quiver(uv_lon[skip_latlon],uv_lat[skip_latlon],
                          u.values[skip_vec],v.values[skip_vec],
                          transform=ccrs.PlateCarree(),
                          scale = 75,
                          headwidth=7,
                          headlength=15,
                          color='k')
                
            elif flag_wind == 0:
                None
            
                
        elif projec == 'AzimuthalEquidistant':
            ax.add_feature(cfeature.LAKES,
                           edgecolor='k',
                           alpha=1,
                           linewidth=0.25,
                           facecolor="none")
            
        ax.add_feature(cfeature.OCEAN,
                       edgecolor='k',
                       alpha=1,
                       linewidth=0,
                       facecolor="lightgray")
        
        ax.outline_patch.set_linewidth(sbplt_lw)
        
        if projec == 'AzimuthalEquidistant':
            ax.set_extent(new_extent,proj)
            
        elif projec == 'LambertConformal':
            ax.set_extent(new_extent)
            
        elif projec == 'PlateCarree':
            ax.set_extent(new_extent)

    cb = mpl.colorbar.ColorbarBase(ax=cbax, cmap=cmap,
                                   norm=norm,
                                   spacing='uniform',
                                   orientation='horizontal',
                                   extend='neither',
                                   ticks=tick_locs,
                                   drawedges=False)
    new_lak_var = 'LMLT'
    if pr_var == 'tp':
        new_pr_var = 'TP'
    elif pr_var == 'et':
        new_pr_var = 'E'
    elif pr_var == 'sf':
        new_pr_var = 'SF'
    
    if flag_coup != 4:
        cb.set_label(cblabel+'('+str(new_lak_var)+','+str(new_pr_var)+'), cl'+'$\geq$'+str(cl_thrsh),
                     size=cbtitle_font)
    elif flag_coup == 4:
        cb.set_label(cblabel+'(['+str(new_lak_var)+'-TAS],'+str(new_pr_var)+'), cl' +'$\geq$'+str(cl_thrsh),
                     size=cbtitle_font)
        
    cb.ax.xaxis.set_label_position('top');
    cb.ax.tick_params(labelcolor=col_cbticlbl,
                      labelsize=tick_font,
                      color=col_cbtic,
                      length=cb_ticlen,
                      width=cb_ticwid,
                      direction='out'); 
    cb.ax.set_xticklabels(tick_labels)
    cb.outline.set_edgecolor(col_cbedg)
    cb.outline.set_linewidth(cb_edgthic)
    
    plt.tight_layout(pad=pad,w_pad=w_pad,h_pad=h_pad)   
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
        
    elif flag_svplt == 1:
        
        if flag_coup != 4:
            
            if flag_coup != 2:
                f.savefig(picDIR+'/' + ds + '_corr_' + lak_var + '_' + pr_var +'.png',bbox_inches='tight',dpi=500)
            
            elif flag_coup == 2:
                f.savefig(picDIR+'/' + ds + '_corr_' + lak_var + '_' + pr_var + '_et_msk.png',bbox_inches='tight',dpi=500)
            
        elif flag_coup == 4:
            
            if flag_wind == 0:
                f.savefig(picDIR+'/' + ds + '_corr_' + lak_var + '-tas_' + pr_var +'.png',bbox_inches='tight',dpi=500)
                
            elif flag_wind == 1:
                f.savefig(picDIR+'/' + ds + '_corr_' + lak_var + '-tas_' + pr_var +'_uv.png',bbox_inches='tight',dpi=500)