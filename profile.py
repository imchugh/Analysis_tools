# Standard modules
import Tkinter, tkFileDialog
from configobj import ConfigObj
import netCDF4
import os
import numpy as np
import pandas as pd
import datetime as dt
from scipy import stats
import re
import matplotlib.pyplot as plt
import pdb

# Custom modules
import DataIO as io

def profile_run():
    
    print 'Select control file containing run parameters'
    cf_path_name=io.file_select_dialog()
    cf=ConfigObj(cf_path_name)
    work_file_path=os.path.join(cf['files']['input_path'],cf['files']['input_file'])
    d=io.xlsx_to_pandas(work_file_path)
    print work_file_path
    return d

def truncate_profile_data():
    
    frequency_mins_out=30
    lag=0 #(calculated lag in number of time steps - system displacement / flow rate)
    bin_size=6 # Number of samples over which to average (defaults to 2 if 0 entered)

    print 'Truncating data to '+str(frequency_mins_out)+' frequency:'
		
    # Import data
    df=pd.read_pickle(os.path.join(path,profile_name))
	
    # Find the time window to be used given the lag and bin size settings above
    temp_df=pd.DataFrame(index=df.index)
    temp_df['mins']=temp_df.index
    temp_df['mins']=temp_df['mins'].apply(lambda x: x.timetuple().tm_min)
    temp_df['mins']=temp_df['mins']%frequency_mins_out # Rewrite minutes to half hourly repetition if desired frequency is 30min
    arr=np.unique(temp_df['mins']) # Array containing all unique values of mins
    select_arr=np.arange(bin_size)-((bin_size-2)/2)+lag # Create window of indices to acceptable values
    valid_arr=arr[select_arr] # Create window of acceptable values
	
    # Create boolean to retrieve central timestamp for each averaged observation
    temp_df['valid_tmstmp']=temp_df['mins']==0
      
    # Create boolean to retrieve appropriate data for averaging interval (bin size)
    if select_arr[0]<0:
        temp_df['valid_cases']=(temp_df['mins']>=valid_arr[0])|(temp_df['mins']<=valid_arr[-1])
    else:
        temp_df['valid_cases']=(temp_df['mins']>=valid_arr[0])&(temp_df['mins']<=valid_arr[-1])
    
    print '    Removing incomplete intervals from beginning and end of time series'
    
    # Trim the ends of the df to remove incomplete intervals
    init_count_list=[0,-1]
    increment_list=[1,-1]
    trim_list=[]
    for i in xrange(2):
        count=init_count_list[i]
        valid_count=0
        while valid_count<bin_size:
            cond=temp_df['valid_cases'].iloc[count]
            count=count+increment_list[i]
            if cond:
                valid_count=valid_count+1
            else:
                valid_count=0
        count=count-(bin_size*increment_list[i])
        trim_list.append(count)	
    temp_df=temp_df[trim_list[0]:trim_list[1]+1]
    df=df.reindex(temp_df.index)
			
    # Create a sequentially numbered grouping variable using the valid case boolean
    df['grp']=np.nan
    step=(frequency_mins_out/2)
    for i in xrange(0,len(df),step):
        df['grp'].iloc[i:i+bin_size]=i

    # Create the output df and do the averaging
    out_df_index=temp_df['Timestamp'][temp_df['valid_tmstmp']]
    out_df=df.groupby('grp').mean()
    out_df.index=out_df_index
    
    print 'Returning data'    
    
    return out_df
	
def process_data():

    # Set parameters
    path='/home/imchugh/Processing/Whroo/Profile data/'
    name_ancillary='slow_profile_corrected.df'
    r=8.3143 # Universal gas constant
    K_con=273.15 # Kelvin conversion
    frequency_mins_out=30
    		
    # Import truncated profile data
    df=pd.read_pickle(os.path.join(path,name))
    	
    # Import slow profile data file that contains temperatures
    ancillary_df=pd.read_pickle(os.path.join(path,name_ancillary))
    ancillary_df=ancillary_df.reindex(df.index)
    	
    # Get the temperature and CO2 column names
    CO2_cols_list=[i for i in df.columns if 'Cc' in i]
    Ta_cols_list=[i for i in ancillary_df.columns if 'Ta' in i]
    
    # Get the actual heights from the column names and calculate layer depths
    heights=[map(float,re.findall(r'[-+]?\d*\.\d+|\d+',i)) for i in CO2_cols_list]
    heights=[i[0] for i in heights]
    layer_thickness=[heights[i+1]-heights[i] for i in xrange(len(heights)-1)]
    layer_thickness=[heights[0]]+layer_thickness
    			
    # Create layer names
    CO2_layer_list=['Cc_'+str(i)+'m' for i in layer_thickness]
    Ta_layer_list=['Ta_'+str(i)+'m' for i in layer_thickness]
    
    # Calculate CO2 averages for each layer
    # (assume value at lowest measurement height is average for that layer (i.e. 0-0.5m))
    for i in xrange(len(CO2_cols_list)):
    	if i==0:
    		df[CO2_layer_list[i]]=df[CO2_cols_list[i]]
    	else:
    		df[CO2_layer_list[i]]=(df[CO2_cols_list[i]]+df[CO2_cols_list[i-1]])/2
    			
    # Calculate temperature averages for each layer
    # (assume value at lowest measurement height is average for that layer (i.e. 0-0.5m))
    for i in xrange(len(Ta_cols_list)):
    	if i==0:
                df[Ta_layer_list[i]]=ancillary_df[Ta_cols_list[i]]
    	else:
    		df[Ta_layer_list[i]]=(ancillary_df[Ta_cols_list[i]]+ancillary_df[Ta_cols_list[i-1]])/2
    				
    # Create the output df
    CO2_dels_list=['del_'+str(i) for i in CO2_cols_list]
    CO2_stor_list=['CO2_stor_'+str(i)+'m' for i in layer_thickness]
    outdf_cols_list=CO2_dels_list+CO2_stor_list+Ta_layer_list+['ps']
    out_df=pd.DataFrame(index=df.index[1:],columns=outdf_cols_list)
    out_df[Ta_layer_list]=df[Ta_layer_list][1:]
    out_df['ps']=ancillary_df['ps'][1:]
    	
    # Calculate CO2 time deltas for each layer
    for i in xrange(len(CO2_cols_list)):
    	out_df[CO2_dels_list[i]]=df[CO2_layer_list[i]]-df[CO2_layer_list[i]].shift()
    		
    # Do the storage calculation for each layer, scaling to umol m-2 s-1 over layer thickness
    for i in xrange(len(CO2_cols_list)):	
    	out_df[CO2_stor_list[i]]=(out_df['ps'][1:]*10**2/(r*(K_con+out_df[Ta_layer_list[i]][1:]))*
    							  out_df[CO2_dels_list[i]]/(frequency_mins_out*60)*layer_thickness[i])
    	
    # Sum all heights and remove sums where any height was nan
    out_df['CO2_stor_tot']=out_df[CO2_stor_list].sum(axis=1)
    out_df['CO2_stor_tot']=out_df['CO2_stor_tot'][~np.isnan(out_df[CO2_stor_list]).any(axis=1)]
    				
    # Output
    out_df.to_pickle(os.path.join(path,'profile_storage.df'))
    out_df.to_csv(os.path.join(path,'profile_storage.csv'),sep=',')
    	
    # Plotting
    CO2_stor_list.append('CO2_stor_tot')
    for i in CO2_stor_list:
        plt.plot(out_df.index,out_df[i],label=i)
        plt.legend(loc='upper left')
        plt.show()
	
