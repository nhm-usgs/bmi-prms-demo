import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import geopandas as gpd
import pandas as pd
import numpy as np


def bmi_prms6_value_plot(data, n_index, val, label, start, end, tax = None):
    tax = tax or plt.gca()
    #test if val exists in both and get nhru or nsegment
    dim_type = None
    try:
        dim_type = data[val].dims[1]

        if dim_type == 'nhru':
            data_val = data[val].sel(nhru=n_index, time=slice(start, end)).to_pandas()
            # dprms_val = dprms[val].sel(nhru=n_index, time=slice(start, end))
            data_val.plot.line(ax=tax, label=label)
            tax.legend()
            # line1, = dprms_val.plot.line(x='time', ax=tax, add_legend=True)

        elif dim_type == 'nsegment':
            data_val = data[val].sel(nsegment=n_index, time=slice(start, end)).to_pandas()
            # dprms_val = dprms[val].sel(nsegment=n_index, time=slice(start, end)).to_pandas()

            data_val.plot(ax=tax, label=label)
            tax.legend()
            # line1, = dprms_val.plot(label='PRMS6')

        tax.set_title(f'{val} {n_index}')

    except Exception as err:
        print('Error', {err})


def bmi_prms6_residual_plot(dbmi, dprms, n_index, val, label, start, end, tax = None):
    tax = tax or plt.gca()
    dim_type = dbmi[val].dims[1]
    try:
        if dim_type == 'nhru':
            data_val = dbmi[val] - dprms[val]
            data = data_val.sel(nhru=n_index, time=slice(start, end)).to_pandas()
            # bmi = dbmi[val]
            # prms = dprms.sel(nhru=n_index, time=slice(start, end))[val]
        elif dim_type == 'nsegment':
            data_val = dbmi[val] - dprms[val]
            data = data_val.sel(nsegment=n_index, time=slice(start, end)).to_pandas()
            # bmi = dbmi.sel[val]
            # prms = dprms.sel(nsegment=n_index, time=slice(start, end))[val]
        # res = prms-bmi
        data.plot(ax=tax, label=label)
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        tax.legend()
        tax.set_title('Residual (prms-bmi)')
    except Exception as err:
        print('Error', {err})


def get_feat_coord(feat, data_set, feat_id):
    lat_da = data_set[feat + '_lat']
    lat = lat_da[feat_id-1].values
    lon_da = data_set[feat + '_lon']
    lon = lon_da[feat_id-1].values
    return lat,lon


def get_hrus_for_box(ds, lat_min, lat_max, lon_min, lon_max):
    sel = ds.hru_lat.sel(hruid=((ds.hru_lat.values >= lat_min)
                            & (ds.hru_lat.values <= lat_max)))
    ids_1 = sel.hruid.values
    sel_1 = ds.hru_lon.sel(hruid=ids_1)
    sel_2 = sel_1.sel(hruid=((sel_1.values >= lon_min) & (sel_1.values <= lon_max)))
    ids_2 = sel_2.hruid.values
    return ids_2


def get_segs_for_box(ds, lat_min, lat_max, lon_min, lon_max):
    sel = ds.seg_lat.sel(segid=((ds.seg_lat.values >= lat_min)
                            & (ds.seg_lat.values <= lat_max)))
    ids_1 = sel.segid.values
    sel_1 = ds.seg_lon.sel(segid=ids_1)
    sel_2 = sel_1.sel(segid=((sel_1.values >= lon_min) & (sel_1.values <= lon_max)))
    ids_2 = sel_2.segid.values
    return ids_2


def get_values_for_DOY(ds, timestamp, hru_ids, var_name):
    if (timestamp < pd.Timestamp('1979-10-01') or timestamp > pd.Timestamp('1980-09-30')):
        print("The date you provided is outside of range 1979-10-01 to 1980-09-30")
        return None
        
    time_range = pd.date_range(timestamp, freq='1Y', periods=40)
    dif = timestamp - time_range[0]
    time_range = time_range + dif
    # print(time_range)

    date_list = []
    val_list = []
    for ts in time_range:
        try:
            date_str = str(ts.year).zfill(4) + '-' + str(ts.month).zfill(2) + '-' + str(ts.day).zfill(2)
            ds_sel = ds[var_name].sel(hruid=hru_ids, time=date_str)
            val = ds_sel.values[0][0]
            date_list.append(date_str + 'T05:00:00')
            val_list.append(val)
        except:
            pass
        
    val_np = np.asarray(val_list, dtype=np.float64)
    val_np = val_np.reshape((1, val_np.shape[0]))
    hru_ids_np = np.asarray(hru_ids, dtype=np.int32)
    date_np = np.asarray(date_list, dtype='datetime64[ns]')
    
    attrs = ds[var_name].attrs
    da_new = xr.DataArray(data=val_np, dims=['hruid','time'],
                          coords={'hruid':hru_ids_np,'time':date_np},
                          attrs=attrs)

    return da_new
    
    
def plot_climate(c_xarray, hru_index, val, start, end, tax=None):
    tax = tax or plt.gca()
    hru_ids = c_xarray.hru.values
    simclimate = c_xarray.sel(time=slice(start, end))

    line, = simclimate.sel(hru=hru_ids[hru_index])[val].plot(ax=tax)
    tax.set_title(val)
    

def bmi_prms6_value_splot(gdf, mbmi, value, tvmin, tvmax, index, timesel, pax = None):
    tax = pax or plt.gca()
    gdf[value] = mbmi.get_value(value)
    divider = make_axes_locatable(tax)
    tcax = divider.append_axes(position='right', size='5%', pad=0.1)
    gdf.plot(column=value, vmin=tvmin, vmax=tvmax, ax=tax, legend=True, cax=tcax)
    tax.set_title(value)


def get_DataSet_prms6(summary, myparam):
    # merge spatial locations of hru and segments into summary file
    ds = xr.open_dataset(summary)
    param = xr.open_dataset(myparam)
    hru_lat = param.get("hru_lat")
    ds['hru_lat'] = hru_lat
    hru_lon = param.get("hru_lon")
    ds['hru_lon'] = hru_lon
    seg_lat = param.get("seg_lat")
    ds['seg_lat'] = seg_lat
    seg_lon = param.get("seg_lon")
    ds['seg_lon'] = seg_lon
    return ds


def get_gdf(file, msurf):
    gdf =gpd.read_file(file)
    pd.set_option('mode.chained_assignment', None)
    gdf_ps = gdf[gdf['hru_id_nat'].isin(msurf.var['nhm_id'].data)]
    dindex = np.zeros(np.shape(gdf_ps.hru_id_nat.values), dtype=np.int8)
    for index, val in np.ndenumerate(msurf.var['nhm_id'].data):
        tind = np.int(np.where(gdf_ps['hru_id_nat'].values == msurf.var['nhm_id'].data[index])[0])
    #     print(type(tind), tind)
        dindex[tind] = np.array([index])
    # print(dindex)
    gdf_ps.loc[:,'tindex'] = dindex
    gdf_ps.sort_values(by=['tindex'], inplace=True)
    return gdf_ps


def get_gdf_streams(file, msurf):
    gdf =gpd.read_file(file)
    pd.set_option('mode.chained_assignment', None)
    gdf_ps = gdf[gdf['seg_id_nat'].isin(msurf.var['nhm_seg'].data)]
    dindex = np.zeros(np.shape(gdf_ps.seg_id_nat.values), dtype=np.int8)
    for index, val in np.ndenumerate(msurf.var['nhm_seg'].data):
        tind = np.int(np.where(gdf_ps['seg_id_nat'].values == msurf.var['nhm_seg'].data[index])[0])
    #     print(type(tind), tind)
        dindex[tind] = np.array([index])
    # print(dindex)
    gdf_ps.loc[:,'tindex'] = dindex
    gdf_ps.sort_values(by=['tindex'], inplace=True)
    return gdf_ps


def plot_climate2(clim_file, gdf_ps, msurf):
    
    clim = xr.open_dataset(clim_file)
    ptime = msurf.var['nowtime'].data
    timesel = dt.datetime(ptime[0], ptime[1], ptime[2])
    start_date = timesel
    gdf_ps['tmax'] = clim.tmax.sel(time=timesel)
    gdf_ps['tmin'] = clim.tmin.sel(time=timesel)
    gdf_ps['prcp'] = clim.prcp.sel(time=timesel)
    fig, ax = plt.subplots(ncols=3)
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    h_tmax = gdf_ps.tmax.max()
    l_tmax = gdf_ps.tmax.min()
    h_tmin= gdf_ps.tmin.max()
    l_tmin= gdf_ps.tmin.min()
    h_tmax = gdf_ps.tmax.max()
    l_tmax = gdf_ps.tmax.min()
    h_ppt= gdf_ps.prcp.max()
    l_ppt= gdf_ps.prcp.min()

    gdf_ps.plot(column='tmax', ax=ax[0], vmin=l_tmax, vmax=h_tmax, legend=True,
                label='tmax',cax=cax0)
    gdf_ps.plot(column='tmin', ax=ax[1], vmin=l_tmin, vmax=h_tmin, legend=True, 
                label='tmin',cax=cax1)
    gdf_ps.plot(column='prcp', ax=ax[2], vmin=l_ppt, vmax=h_ppt, legend=True,
                label='prcp',cax=cax2)
    for i in range(3):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if i == 0:
            ax[i].set_title('tmax')
        elif i == 1:
            ax[i].set_title('tmin')
        elif i == 2:
            ax[i].set_title('prcp')
    plt.tight_layout()
    
    return clim


def example_plot(clim, gdf_ps, msurf, msoil, j, timesel):
    gdf_ps['tmax'] = clim.tmax.sel(time=timesel)
    gdf_ps['tmin'] = clim.tmin.sel(time=timesel)
    gdf_ps['prcp'] = clim.prcp.sel(time=timesel)

    gdf_ps['infil'] = msurf.var['infil'].data
    gdf_ps['sroff'] = msurf.var['sroff'].data
    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data

    fig, ax = plt.subplots(ncols=6, figsize = (12,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=20.0, vmax=65.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='tmin', vmin=20.0, vmax=65.0, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=0.7, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='infil', vmin=0.0, vmax=0.7, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='sroff', vmin=0.0, vmax=0.25, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=1.75, ax=ax[5], legend=True, cax=cax5)
    for i in range(6):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('tmin')
            elif i == 2:
                ax[i].set_title('prcp')
            elif i == 3:
                ax[i].set_title('infil')
            elif i == 4:
                ax[i].set_title('sroff')
            elif i == 5:
                ax[i].set_title('soil_moist_tot')
    plt.tight_layout()


def example_plot_strm(clim, gdf_ps, gdf_stream, msurf, msoil, mgw, mstrm, j, timesel):
    gdf_ps['tmax'] = (msurf.get_value('tmax')*(9.0/5.0)) + 32.0
    gdf_ps['tmin'] = (msurf.get_value('tmin')*(9.0/5.0)) + 32.0
    gdf_ps['prcp'] = msurf.get_value('hru_ppt')

    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data
    gdf_ps['sroff'] = msurf.var['sroff'].data
    gdf_ps['ssres_flow'] = msoil.var['ssres_flow'].data
    gdf_ps['gwres_flow'] = mgw.var['gwres_flow'].data
    gdf_stream['seg_outflow'] = mstrm.var['seg_outflow'].data

    fig, ax = plt.subplots(ncols=7, figsize = (14,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
    divider6 = make_axes_locatable(ax[6])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
    cax6 = divider6.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=40.0, vmax=85.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=0.7, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=3.0, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='sroff', vmin=0.0, vmax=0.15, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='ssres_flow', vmin=0.0, vmax=0.05, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='gwres_flow', vmin=0.0, vmax=0.02, ax=ax[5], legend=True, cax=cax5)
    gdf_stream.plot(column='seg_outflow', vmin=0.0, vmax=200, ax=ax[6], legend=True, cax=cax6)

        
    
    for i in range(7):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('prcp')
            elif i == 2:
                ax[i].set_title('soil_moist_tot')
            elif i == 3:
                ax[i].set_title('sroff')
            elif i == 4:
                ax[i].set_title('ssres_flow')
            elif i == 5:
                ax[i].set_title('gwres_flow')
            elif i == 6:
                ax[i].set_title('seg_outflow')
    plt.tight_layout()
    
    


def gm_example_plot(gdf_ps, gmdata, msurf, msoil, j, timesel):
    gdf_ps['tmax'] = (gmdata.tmax.data[j,:]*(9/5))+32.0
    gdf_ps['tmin'] = (gmdata.tmin.data[j,:]*(9/5))+32.0
    gdf_ps['prcp'] = gmdata.precip.data[j,:]*.0393701

    gdf_ps['infil'] = msurf.var['infil'].data
    gdf_ps['sroff'] = msurf.var['sroff'].data
    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data

    fig, ax = plt.subplots(ncols=6, figsize = (12,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
    # divider6 = make_axes_locatable(ax[6])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
    # cax6 = divider6.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=50.0, vmax=70.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='tmin', vmin=20.0, vmax=45.0, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=.75, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='infil', vmin=0.0, vmax=0.7, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='sroff', vmin=0.0, vmax=0.25, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=1.75, ax=ax[5], legend=True, cax=cax5)
    # gdf_ps.plot(column='soil_moist_tot', vmin=0.0, vmax=1.5, ax=ax[6], legend=True, cax=cax6)
    for i in range(6):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('tmin')
            elif i == 2:
                ax[i].set_title('prcp')
            elif i == 3:
                ax[i].set_title('soil_to_gw')
            elif i == 4:
                ax[i].set_title('ssr_to_gw')
            elif i == 5:
                ax[i].set_title('soil_moist_tot')
            # elif i == 6:
            #     ax[i].set_title('soil_moist_tot')
    plt.tight_layout()


def example_plot2(gdf_ps, msurf, msoil, j, timesel):
    gdf_ps['tmax'] = (msurf.get_value('tmax')*(9.0/5.0)) + 32.0
    gdf_ps['tmin'] = (msurf.get_value('tmin')*(9.0/5.0)) + 32.0
    gdf_ps['prcp'] = msurf.get_value('hru_ppt')

    gdf_ps['soil_to_gw'] = msoil.var['soil_to_gw'].data
    gdf_ps['ssr_to_gw'] = msoil.var['ssr_to_gw'].data
    gdf_ps['ssres_flow'] = msoil.var['ssres_flow'].data
    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data

    fig, ax = plt.subplots(ncols=7, figsize = (14,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
    divider6 = make_axes_locatable(ax[6])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
    cax6 = divider6.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=20.0, vmax=75.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='tmin', vmin=20.0, vmax=75.0, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=0.7, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='soil_to_gw', vmin=0.0, vmax=0.1, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='ssr_to_gw', vmin=0.0, vmax=0.15, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='ssres_flow', vmin=0.0, vmax=0.1, ax=ax[5], legend=True, cax=cax5)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=3.0, ax=ax[6], legend=True, cax=cax6)
    for i in range(7):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('tmin')
            elif i == 2:
                ax[i].set_title('prcp')
            elif i == 3:
                ax[i].set_title('soil_to_gw')
            elif i == 4:
                ax[i].set_title('ssr_to_gw')
            elif i == 5:
                ax[i].set_title('ssres_flow')
            elif i == 6:
                ax[i].set_title('soil_moist_tot')
    plt.tight_layout()