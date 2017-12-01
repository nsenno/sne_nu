# Code that will read in the data from ../data/raw_sne_data.csv and manipulate it into 
# a form for the coincidence analysis. See accompanying jupyter notebook for detailed description 

import pandas as pd 
import numpy as np

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import Angle

# Define function to calculate the average values of R.A. and Dec. 
# ang_unit is a kwarg which describes how the angular data is presented 
#
#   R.A. -- u.hourangle
#   Dec. -- u.degree


def get_average_angle(ang, ang_unit = u.degree) : 
    return np.mean([y for y in map(lambda x : Angle(x,ang_unit).rad,ang.split(','))])

# Define function to calculate the average value of redshift (z)
# Most values are strings that are separated by ','
# However, NaNs are floats. A -1 is returned to signal NA data 

def get_average_z(z) : 
    return np.mean(list(map(float,z.split(',')))) if type(z) is str else -1 

# read in the SNe data 

sne_data = pd.read_csv('data/raw_sne_data.csv',index_col='Name')

# convert the Gregorian dates to MJD using Time from astropy.time package 
# add a new column to the df with MJD dates 

sne_data['Max Date (MJD)'] = Time([date for date in sne_data['Max Date'].str.replace('/','-')]).mjd 

# Clean the R.A., Dec., and z values with one list comprehension 
# Create new columns for R.A. and Dec. in rad
# Replace the z column with the new values 

sne_data['R.A. (rad)'], sne_data['Dec. (rad)'], sne_data['z'] = list(zip(*[(get_average_angle(row['R.A.'],\
        u.hourangle), get_average_angle(row['Dec.']), get_average_z(row['z'])) for _, row in sne_data.iterrows()]))

bad_types = ['SLSN-II','II']
right_type = [sne_type not in bad_types for sne_type in sne_data.Type]


#relevant time window set by neutrino data 

min_mjd = 55750
max_mjd = 56068 

in_time_window = (sne_data['Max Date (MJD)'] >= min_mjd) & (max_mjd >= sne_data['Max Date (MJD)'])

# Write the data to the file 

sne_data.loc[in_time_window*right_type, ['Max Date (MJD)','R.A. (rad)','Dec. (rad)','z','Type']].sort_values( \
        'Max Date (MJD)').to_csv('cleaned_sne_data.csv')
