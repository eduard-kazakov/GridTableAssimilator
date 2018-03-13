import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime


def CreateTable(spatial_mask_npy, dates, gridded_params, single_params_path=None):
    spatial_mask = np.load(spatial_mask_npy)

    general_sub_tables = []

    for current_gridded_param in gridded_params:
        print current_gridded_param['name']
        all_param_files = glob.glob(os.path.join(current_gridded_param['path'], '*.*'))

        valid_param_files = []

        # DATE FILTER
        for param_file in all_param_files:
            try:
                current_date = datetime.strptime(os.path.basename(param_file), current_gridded_param['naming_mask'])
            except:
                continue

            if dates['mode'] == 'from_to':
                if dates['start_date'] <= current_date <= dates['end_date']:
                    valid_param_files.append(param_file)

            elif dates['mode'] == 'list':
                if current_date in dates['list_of_dates']:
                    valid_param_files.append(param_file)

            elif dates['mode'] == 'set':

                active_years = dates['years']
                active_months = dates['months']
                active_days = dates['days']

                if not active_years:
                    active_years = range(1900, 2100, 1)
                if not active_months:
                    active_months = range(1, 13, 1)
                if not active_days:
                    active_days = range(1, 32, 1)

                if (current_date.year in active_years) and (current_date.month in active_months) and (
                    current_date.day in active_days):
                    valid_param_files.append(param_file)

            else:
                continue

        # MASKING AND PUT TO TABLE
        sub_tables = []
        for valid_param_file in valid_param_files:
            current_date = datetime.strptime(os.path.basename(valid_param_file), current_gridded_param['naming_mask'])
            current_grid = np.load(valid_param_file)
            current_indices = np.indices(current_grid.shape)

            current_grid = current_grid[spatial_mask == 1]
            current_rows = current_indices[0][spatial_mask == 1]
            current_cols = current_indices[1][spatial_mask == 1]

            current_sub_table = pd.DataFrame(
                {current_gridded_param['name']: current_grid, 'ROW': current_rows, 'COL': current_cols})
            current_sub_table['PERIOD'] = current_date

            sub_tables.append(current_sub_table)

        if sub_tables:
            general_param_sub_table = pd.concat(sub_tables)
            general_sub_tables.append(general_param_sub_table)

    i = 0

    final_table = general_sub_tables[0]
    while i < len(general_sub_tables):
        if i == 0:
            i += 1
            continue

        final_table = pd.merge(final_table, general_sub_tables[i], on=['PERIOD', 'ROW', 'COL'], how='outer')
        i += 1

    if single_params_path:
        single_params = pd.read_csv(single_params_path)
        single_params['PERIOD'] = pd.to_datetime(single_params['PERIOD'], format='%Y%m%d')  # .astype(pd.datetime)
        single_params = single_params[(single_params['PERIOD'] >= start_date) & (single_params['PERIOD'] <= end_date)]
        final_table = pd.merge(final_table, single_params, on=['PERIOD'], how='outer')

    return final_table

'''
### HOW To use
### First step: define spatial region as binary npy mask (1 -> use location, 0 -> not)
spatial_mask_npy = '/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/Sea_Masks/npy/barentsSea.npy'

### Second step: define input grid locations

gridded_params = [{'name':'BloomMask','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/spectra_mask/','naming_mask':'ESACCI-OC-L3S-RRS-MERGED-8D_DAILY_4km_GEO_PML_RRS-%Y%m%d-fv3.1.npy'},
                  {'name':'CoccoConcentraton','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/cct/','naming_mask':'%Y%m%d_cct.npy'},
                  {'name':'NO3_WOD','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/NO3/reprojected/WOD13/','naming_mask':'%Y%m%d.npy'},
                  {'name':'NO3_GLODAP','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/NO3/reprojected/GLODAP/','naming_mask':'%Y%m%d.npy'},
                  {'name':'SST','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/sst/reprojected/pathfinder5.3/8days/gapfilling/','naming_mask':'%Y%m%d.npy'},
                  {'name':'PAR','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/par/reprojected/renamed/','naming_mask':'%Y%m%d_par_mean.npy'},
                  {'name':'Salinity','path':'/mnt/nfs0/RSF2017_PozdnyakovDV/data/gridded/salinity/reprojected/SMOS/8days/','naming_mask':'%Y%m%d.npy'}]

### Third step: define dates. Three ways to do it: 'from_to', 'list', 'set'

dates_v1 = {'mode':'from_to','start_date':datetime(2011,8,1),'end_date':datetime(2011,8,30)}
dates_v2 = {'mode':'list','list_of_dates':[datetime(2011,1,1),datetime(2011,1,9)]}
dates_v3 = {'mode':'set','years':[2002],'months':[8],'days':[]}

### Fourth step: define path to single data, if needed

single_params = '/mnt/nfs0/data_sonarc/data/heap/INDICIES_Processing/INDICIES_Processed_v2.csv'

### Fifth step: RUN!
data = CreateTable(spatial_mask_npy,dates_v3,gridded_params,None)
'''