# -*- coding: utf-8 -*-
import datetime
import re
import urllib
from pathlib import Path
import pandas as pd

import requests
import xarray as xr
import numpy as np
import yaml

from helpers import np_get_wval, getaverage


class Gridmet:
    SCHEME = "http"
    NETLOC = "thredds.northwestknowledge.net:8080"
    PATH = {
        "daily_maximum_temperature": "thredds/ncss/agg_met_tmmx_1979_CurrentYear_CONUS.nc",
        "daily_minimum_temperature": "thredds/ncss/agg_met_tmmn_1979_CurrentYear_CONUS.nc",
        "precipitation_amount": "thredds/ncss/agg_met_pr_1979_CurrentYear_CONUS.nc",
    }

    def __init__(self, start_date=None, end_date=None, hrumap=None, hru_id=None, wght_file=None, lazy=True,
                 cache_dir=None, config_file=None):
        self._wghts = None
        self._wghts_id = None
        lazy = True
        cache_dir = None
        if config_file is not None:
            with open(config_file, 'r') as fp:
                parameters = yaml.safe_load(fp)
            for key, value in parameters.items():
                setattr(self, key, value)
        else:
            self._start_date = Gridmet.datetime_or_yesterday(start_date)
            self._end_date = Gridmet.datetime_or_yesterday(end_date)
            self.hru_id = hru_id
            self.wght_file = wght_file
            self.return_map = hrumap


        if self._start_date > self._end_date:
            raise ValueError(
                "start date ({0}) must be before end date ({1})".format(
                    self.start_date, self.end_date
                )
            )
        if self._end_date > datetime.date.today():
            raise ValueError(
                "end date cannot be a future date ({0} > {1}".format(
                    self.end_date, datetime.date.today()
                )
            )
        self._delta = (self._end_date - self._start_date)
        self._m_tmin_data = None
        self._m_tmax_data = None
        self._m_prcp_data = None
        if self.return_map:
            if self.wght_file:
                self._wghts = pd.read_csv(self.wght_file)
                self._wghts_id = self._wghts.columns[1]
                self._unique_hru_ids = self._wghts.groupby(self._wghts_id)
                self._m_tmin_data = np.zeros(shape=(self._delta.days + 1, len(self.hru_id)))
                self._m_tmax_data = np.zeros(shape=(self._delta.days + 1, len(self.hru_id)))
                self._m_prcp_data = np.zeros(shape=(self._delta.days + 1, len(self.hru_id)))
            else:
                self.return_map = False
                print('mapping to hru ids requires weights file')

        if cache_dir is None:
            cache_dir = Path("~/.gridmet")
        self._cache_dir = Path(cache_dir).expanduser().resolve()

        self._dataset = None
        if not lazy:
            for name in self.PATH:
                self._lazy_load(name)

    @staticmethod
    def clear_cache(cache_dir=None):
        for fname in Gridmet.list_cache(cache_dir=cache_dir):
            fname.unlink()

    @staticmethod
    def list_cache(cache_dir=None):
        if cache_dir is None:
            cache_dir = Path("~/.gridmet")
        cache_dir = Path(cache_dir).expanduser().resolve()

        pattern = r"(?P<var>[a-z_]*)_(?P<start>[0-9\-]*)_(?P<end>[0-9\-]*)\.nc"

        cached_files = []
        for fname in [p.name for p in cache_dir.glob("*.nc")]:
            match = re.match(pattern, fname)
            if match and match.group("var") in Gridmet.PATH:
                try:
                    datetime.date.fromisoformat(match.group("start"))
                    datetime.date.fromisoformat(match.group("end"))
                except ValueError:
                    pass
                else:
                    cached_files.append(cache_dir / fname)

        return cached_files

    @staticmethod
    def datetime_or_yesterday(val):
        if val is None:
            return datetime.date.today() - datetime.timedelta(days=1)
        elif isinstance(val, str):
#             return datetime.date.fromisoformat(val)
            format_str = '%Y-%m-%d'
            return datetime.datetime.strptime(val, format_str).date()
        else:
            return val

    @classmethod
    def from_today(cls, days, lazy=True):
        if days <= 0:
            raise ValueError("number of days must be positive ({0})".format(days))

        end_date = datetime.date.today()
        start_date = end_date - datetime.timedelta(days=days)

        return cls(start_date, end_date, lazy=lazy)

    @property
    def cache_dir(self):
        return self._cache_dir

    @property
    def start_date(self):
        return str(self._start_date)

    @property
    def end_date(self):
        return str(self._end_date)

    @property
    def dataset(self):
        return self._dataset

    def _fetch_and_open(self, name):
        self._cache_dir.mkdir(exist_ok=True)
        return xr.open_dataset(
            Gridmet.fetch_var(
                name, self.start_date, self.end_date, cache_dir=self._cache_dir
            )
        )

    def _lazy_load(self, name):
        if self._dataset is None:
            self._dataset = self._fetch_and_open(name)

        try:
            self._dataset[name]
        except KeyError:
            self._dataset = self._dataset.merge(self._fetch_and_open(name))

        return self._dataset[name]

    @property
    def tmax(self):
        tname = "daily_maximum_temperature"
        ds = self._lazy_load(tname)

        if self.return_map:
            flt_val = ds.values.flatten(order='K')
            for i in np.arange(ds.coords['day'].size):
                for j in np.arange(len(self.hru_id)):
                    weight_id_rows = self._unique_hru_ids.get_group(self.hru_id[j])
                    tw = weight_id_rows.w.values
                    tgid = weight_id_rows.grid_ids.values
                    if np.isnan(getaverage(flt_val[tgid], tw)):
                        self._m_tmax_data[i, j] = getaverage(ds.values[i, :, :].flatten(order='K')[tgid], tw) - 273.15
                    else:
                        self._m_tmax_data[i, j] = np_get_wval(ds.values[i, :, :].flatten(order='K')[tgid], tw) - 273.15

            return xr.DataArray(self._m_tmax_data, dims=['day', 'hru_id'],
                                coords={'day': list(ds.day.coords['day'].values),
                                        'hru_id': list(self.hru_id)}, )
        else:
            return ds

    @property
    def tmin(self):
        tname = "daily_minimum_temperature"
        ds = self._lazy_load(tname)

        if self.return_map:
            flt_val = ds.values.flatten(order='K')
            for i in np.arange(ds.coords['day'].size):
                for j in np.arange(len(self.hru_id)):
                    weight_id_rows = self._unique_hru_ids.get_group(self.hru_id[j])
                    tw = weight_id_rows.w.values
                    tgid = weight_id_rows.grid_ids.values
                    if np.isnan(getaverage(flt_val[tgid], tw)):
                        self._m_tmin_data[i, j] = getaverage(ds.values[i, :, :].flatten(order='K')[tgid], tw) - 273.15
                    else:
                        self._m_tmin_data[i, j] = np_get_wval(ds.values[i, :, :].flatten(order='K')[tgid], tw) - 273.15

            return xr.DataArray(self._m_tmin_data, dims=['day', 'hru_id'],
                                coords={'day': list(ds.day.coords['day'].values),
                                        'hru_id': list(self.hru_id)}, )
        else:
            return ds

    @property
    def precip(self):
        tname = "precipitation_amount"
        ds = self._lazy_load(tname)
        if self.return_map:
            flt_val = ds.values.flatten(order='K')
            for i in np.arange(ds.coords['day'].size):
                for j in np.arange(len(self.hru_id)):
                    weight_id_rows = self._unique_hru_ids.get_group(self.hru_id[j])
                    tw = weight_id_rows.w.values
                    tgid = weight_id_rows.grid_ids.values
                    if np.isnan(getaverage(flt_val[tgid], tw)):
                        self._m_prcp_data[i, j] = getaverage(ds.values[i, :, :].flatten(order='K')[tgid], tw)
                    else:
                        self._m_prcp_data[i, j] = np_get_wval(ds.values[i, :, :].flatten(order='K')[tgid], tw)

            return xr.DataArray(self._m_prcp_data, dims=['day', 'hru_id'],
                                coords={'day': list(ds.day.coords['day'].values),
                                        'hru_id': list(self.hru_id)}, )
        else:
            return ds

    @classmethod
    def fetch_var(cls, name, start_date, end_date=None, cache_dir="."):
        if name not in cls.PATH:
            raise ValueError(
                "name not understood ({0} not in {1})".format(name, ", ".join(cls.PATH))
            )

        end_date = end_date or datetime.date.today()
        fname = Path(cache_dir) / "{var}_{start}_{end}.nc".format(
            var=name, start=start_date, end=end_date
        )

        if not fname.is_file():
            params = {
                "var": name,
                "north": "49.4000",
                "west": "-124.7666",
                "east": "-67.0583",
                "south": "25.0666",
                "disableLLSubset": "on",
                "disableProjSubset": "on",
                "horizStride": "1",
                "time_start": start_date + "T00:00:00Z",
                "time_end": end_date + "T00:00:00Z",
                "timeStride": "1",
                "accept": "netcdf",
            }

            response = requests.get(cls.data_url(name), params=params)
            response.raise_for_status()

            with fname.open("wb") as fp:
                fp.write(response.content)

        return fname.absolute()

    @classmethod
    def data_url(cls, name):
        return urllib.parse.urlunparse(
            (cls.SCHEME, cls.NETLOC, cls.PATH[name], "", "", "")
        )
