import pandas as pd 
import numpy as np

nu_data = pd.read_csv('data/upgoing_events.txt',delim_whitespace=True)

#convert data from degrees to radians 

nu_data['RA (rad)'] = nu_data.RA.apply(np.radians)
nu_data['Dec (rad)'] = nu_data.Dec.apply(np.radians)
nu_data['ang_err (rad)'] = nu_data.ang_err.apply(np.radians)

#write new data to file with correct column ordering

nu_data[['MJD','ang_err (rad)','RA (rad)','Dec (rad)']].to_csv('data/cleaned_nu_data.csv',index=False)
