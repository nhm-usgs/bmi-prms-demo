import numpy as np
from numpy.ma import masked


def getaverage(data, wghts):
    try:
        v_ave = np.average(data, weights=wghts)
    except ZeroDivisionError:
        v_ave = np.nan
    return v_ave


def np_get_wval(ndata, wghts, hru_id=None):
    mdata = np.ma.masked_array(ndata, np.isnan(ndata))
    tmp = np.ma.average(mdata, weights=wghts)
    if tmp is masked:
        # print(f'returning masked value: {hru_id}', ndata)
        return np.nan
    else:
        return tmp
