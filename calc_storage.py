import netCDF4
import os
import numpy as np
import pandas as pd
import datetime as dt
from scipy import stats
import re
import matplotlib.pyplot as plt
import pdb

def stitch_data():

	# Set parameters
	path='C:\Offline_data_processing\Data\Whroo\Profile data'
	name1='Whroo_profile_IRGA_avg.dat.backup'
	name2='Whroo_profile_IRGA_avg.dat'
	name_ancillary1='Whroo_profile_Slow_avg.dat.backup'
	name_ancillary2='Whroo_profile_Slow_avg.dat'
	
	### Primary data ###
	
	# Import and clean up first dataset
	data1_df=pd.read_csv(os.path.join(path,name1),skiprows=[0,2,3])
	data1_df.index=pd.to_datetime(data1_df['TIMESTAMP'])
	data1_df.drop(['TIMESTAMP','RECORD'],axis=1,inplace=True)
	data1_df.drop_duplicates(inplace=True)
	data1_df=data1_df.astype(float)

	# Change the earliest data to two minute frequency to match later data
	sub_df=data1_df[:'2012-02-28 12:02:00']
	new_index=pd.date_range(start=sub_df.index[0],end=sub_df.index[-1],freq='1Min')
	sub_df=sub_df.reindex(new_index,inplace=True)
	sub_df['grp']=np.arange(len(sub_df))/2
	new_df=sub_df.groupby('grp').mean()
	new_index=pd.date_range(start=sub_df.index[1],end=sub_df.index[-1],freq='2Min')
	new_df.index=new_index

	# Import and clean up second dataset
	data2_df=pd.read_csv(os.path.join(path,name2),skiprows=[0,2,3])
	data2_df.index=pd.to_datetime(data2_df['TIMESTAMP'])
	data2_df.drop(['TIMESTAMP','RECORD'],axis=1,inplace=True)
	data2_df.drop_duplicates(inplace=True)
	data2_df=data2_df.astype(float)
				
	# Put all of the data together and check no gaps
	full_df=pd.concat([new_df,data1_df.ix['2012-02-28 12:10:00':],data2_df],axis=0)
	full_df.drop_duplicates(inplace=True)
	new_index=pd.date_range(start=full_df.index[1],end=full_df.index[-1],freq='2Min')
	full_df=full_df.reindex(new_index,inplace=True)

	### Ancillary data ###
	
	# Import first slow profile dataset to get temperature
	ancillary1_df=pd.read_csv(os.path.join(path,name_ancillary1),skiprows=[0,2,3])
	ancillary1_df.index=pd.to_datetime(ancillary1_df['TIMESTAMP'])
	ancillary1_df.drop(['TIMESTAMP','RECORD'],axis=1,inplace=True)
	ancillary1_df=ancillary1_df.astype(float)
	
	# Import second slow profile dataset to get temperature
	ancillary2_df=pd.read_csv(os.path.join(path,name_ancillary2),skiprows=[0,2,3])
	ancillary2_df.index=pd.to_datetime(ancillary2_df['TIMESTAMP'])
	ancillary2_df.drop(['TIMESTAMP','RECORD'],axis=1,inplace=True)
	ancillary2_df=ancillary2_df.astype(float)
	
	# Join them and dump duplicates
	ancillary_df=pd.concat([ancillary1_df,ancillary2_df])
	ancillary_df.drop_duplicates(inplace=True)
	
	### Output all ###
	full_df.to_pickle(os.path.join(path,'profile_uncorrected.df'))
	ancillary_df.to_pickle(os.path.join(path,'slow_profile_uncorrected.df'))
			
# to do:	
# 	- Use Cape Grimm data as southern hemisphere baseline for CO2
#	- include daytime maximum and nocturnal minimum data for estimating CO2 baseline
	
def correctCO2_data():
	
	# Set parameters
	path='C:\Offline_data_processing\Data\Whroo\Profile data'
	name='profile_uncorrected.df'
	baddata_dates=['2013-08-24','2013-10-29']
	badcoeff_dates=['2012-06-28 11:00:00','2012-10-17 12:50:00']
	instrument_dates=[['2011-12-02','2012-06-28'],
					  ['2012-06-28','2012-10-13'],
					  ['2012-10-13','2013-08-23'],
					  ['2013-10-29','2014-06-02']]
	instrument_times=[['12:00:00','10:58:00'],
					  ['11:00:00','12:00:00'],
					  ['12:02:00','23:58:00'],
					  ['12:00:00','23:58:00']]
	coeff_correct=2.5
	true_heights=[0.5,2,4,8,16,32]
	CO2_range=[300,600]
	CO2_base=390
		
	# Import data
	df=pd.read_pickle(os.path.join(path,name))
	
	# Get the CO2 columns
	CO2_cols_list=[i for i in df.columns if 'Cc' in i]

	# Correct the data for 1) no data; 2) wrong instrument scaling coefficients; 3) range checks; 4) drift; 5) reversed label assignment of CO2
		
	# 1 and 2 above
	for i in CO2_cols_list:
		df.ix[baddata_dates[0]:baddata_dates[1]]=np.nan
		df[i].ix[badcoeff_dates[0]:badcoeff_dates[1]]=df[i].ix[badcoeff_dates[0]:badcoeff_dates[1]]*coeff_correct
		
	# 3 above
	for i in CO2_cols_list:
		df[i]=np.where(df[i]<CO2_range[0],np.nan,np.where(df[i]>CO2_range[1],np.nan,df[i]))
	
	# 4 above (use approximate 2011 CO2 baseline to anchor beginning of series (~390ppm))
	
	# Use nocturnal minima and daytime maxima to estimate the drift for each instrument residence (only need one level)
	noct_df=df[CO2_cols_list].ix[df.index.indexer_between_time(dt.time(20),dt.time(6),include_start=True,include_end=False)]
	daily_df=pd.DataFrame({'Count':noct_df['Cc_LI840_1m'].groupby([lambda x: x.year,lambda y: y.dayofyear]).count(),
						   'Min':noct_df['Cc_LI840_1m'].groupby([lambda x: x.year,lambda y: y.dayofyear]).min()}).reset_index()
	daily_df['Min']=np.where(daily_df['Count']<300,np.nan,daily_df['Min'])
	daily_df.index=(daily_df['level_0'].apply(lambda x: dt.datetime(x,1,1))+
                    daily_df['level_1'].apply(lambda x: dt.timedelta(int(x))))
	daily_df['Drift_correct']=np.zeros(len(daily_df))
	params=[]
	for i,d in enumerate(instrument_dates):
		temp_df=pd.DataFrame({'x':np.arange(len(daily_df.ix[d[0]:d[1]])),
							  'y':daily_df['Min'].ix[d[0]:d[1]]})
		temp_df.dropna(how='any',axis=0,inplace=True)
		params.append(np.polyfit(temp_df['x'],temp_df['y'],1))
	
	# Correct the data
	df['Drift_correct']=np.nan
	for i in xrange(len(instrument_dates)):
		dt_tm_start=instrument_dates[i][0]+' '+instrument_times[i][0]
		dt_tm_end=instrument_dates[i][1]+' '+instrument_times[i][1]
		df['Drift_correct'].ix[dt_tm_start:dt_tm_end]=(np.arange(len(df.ix[dt_tm_start:dt_tm_end]))
													   *(-params[i][0]/720)-(params[i][1]-CO2_base))	 
	for j in CO2_cols_list:
		df[j]=df[j]+df['Drift_correct']
	
	# Drop the correction variable
	df.drop('Drift_correct',axis=1,inplace=True)
	
	# 5 above (reverse the column names and reassign the heights)
	
	# Rename the columns using the true heights
	CO2col_indices_list=[df.columns.get_loc(CO2_cols_list[0]),df.columns.get_loc(CO2_cols_list[-1])]
	true_heights.reverse()
	new_CO2_cols_list=['Cc_'+str(i)+'m' for i in true_heights]
	new_labels=dict([tuple((CO2_cols_list[i],new_CO2_cols_list[i])) for i in xrange(len(CO2_cols_list))])
	df=df.rename(columns=new_labels)
	
	# Reorder the columns so that they are ascending
	new_CO2_cols_list.reverse()
	all_cols_list=df.columns.tolist()
	all_cols_list[CO2col_indices_list[0]:CO2col_indices_list[1]+1]=new_CO2_cols_list
	df=df.reindex_axis(all_cols_list,axis=1)
	
	# Output the data	
	df.to_pickle(os.path.join(path,'profile_corrected.df'))

def correctTa_data():
	
	# Set parameters
	path='C:\Offline_data_processing\Data\Whroo\Profile data'
	name='slow_profile_uncorrected.df'
	baddata1m_dates=['2012-10-20','2013-04-10']
	baddata8m_dates=['2013-07-18','2014-05-08']
	baddata32m_dates=[['2011-12-02','2012-06-28'],
					  ['2013-02-01','2013-04-10']]
	true_heights=[0.5,2,4,8,16,32]
	range_limits=[-20,50]
	
	#Import the data
	df=pd.read_pickle(os.path.join(path,name))
	
	# Find the temperature column names
	Ta_cols_list=[i for i in df.columns if 'Ta' in i]
	
	# Remove known bad data
	df[Ta_cols_list[0]].ix[baddata1m_dates[0]:baddata1m_dates[1]]=np.nan
	df[Ta_cols_list[3]].ix[baddata8m_dates[0]:baddata8m_dates[1]]=np.nan
	for i in xrange(len(baddata32m_dates)):
		df[Ta_cols_list[5]].ix[baddata32m_dates[i][0]:baddata32m_dates[i][1]]=np.nan
	
	# remove data outside range limits
	for i in Ta_cols_list:
		df[i]=np.where(df[i]<range_limits[0],np.nan,np.where(df[i]>range_limits[1],np.nan,df[i]))
			
	# Rename the columns using the true heights
	Tacol_indices_list=[df.columns.get_loc(Ta_cols_list[0]),df.columns.get_loc(Ta_cols_list[-1])]
	new_Ta_cols_list=['Ta'+str(i)+'m' for i in true_heights]
	new_labels=dict([tuple((Ta_cols_list[i],new_Ta_cols_list[i])) for i in xrange(len(Ta_cols_list))])
	df=df.rename(columns=new_labels)
	
	# Fill temperature data for each level from level with best correlation
	for i in new_Ta_cols_list:
		vars_list=list(new_Ta_cols_list)
		vars_list.remove(i)
		results_df=pd.DataFrame(columns=['slope','int','r2','p_val','se'],index=vars_list)
		for j in vars_list:
			temp_df=df[[i,j]].dropna(how='any',axis=0)
			results_df.ix[j]=stats.linregress(temp_df[j],temp_df[i])
		results_df.sort('r2',ascending=False,inplace=True)
		vars_list=results_df.index
		for j in vars_list:
			if len(df[j].dropna())==len(df):
				break
			temp_S=df[j]*results_df['slope'].ix[j]+results_df['int'].ix[j]
			df[i]=np.where(pd.isnull(df[i]),temp_S,df[i])
	
	# Fix the obviously wrong temperature near ground level (due to bug in CSI AM25T instruction)
	df[new_Ta_cols_list[0]]=df[new_Ta_cols_list[0]]+(df[new_Ta_cols_list[1]].mean()-df[new_Ta_cols_list[0]].mean())
	
	# Output data
	df.to_pickle(os.path.join(path,'slow_profile_corrected.df'))
		
def truncate_data():
	
	# Set parameters
	path='C:\Offline_data_processing\Data\Whroo\Profile data'
	name='profile_corrected.df'
	frequency_mins_out=30
	lag=0 #(calculated lag in number of time steps - system displacement / flow rate)
	bin_size=6 # Number of samples over which to average (defaults to 2 if 0 entered)
		
	# Import data
	df=pd.read_pickle(os.path.join(path,name))
	
	# Find the time window to be used given the lag and bin size settings above
	temp_df=pd.DataFrame({'Timestamp':df.index})
	temp_df.index=temp_df['Timestamp']
	temp_df['mins']=temp_df['Timestamp'].apply(lambda x: x.timetuple().tm_min)
	temp_df['mins']=temp_df['mins']%frequency_mins_out # Rewrite minutes to half hourly repetition if desired frequency is 30min
	arr=np.unique(temp_df['mins']) # Array containing all unique values of mins
	select_arr=np.arange(bin_size)-((bin_size-2)/2)+lag # Create window of indices to acceptable values
	valid_arr=arr[select_arr] # Create window of acceptable values
		
	# Create boolean to retrieve central timestamp for each averaged observation
	temp_df['valid_tstmp']=temp_df['mins']==0
	
	# Create boolean to retrieve appropriate data for averaging interval (bin size)
	if select_arr[0]<0:
		temp_df['valid_cases']=(temp_df['mins']>=valid_arr[0])|(temp_df['mins']<=valid_arr[-1])
	else:
		temp_df['valid_cases']=(temp_df['mins']>=valid_arr[0])&(temp_df['mins']<=valid_arr[-1])
	
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
	out_df_index=temp_df['Timestamp'][temp_df['valid_tstmp']]
	out_df=df.groupby('grp').mean()
	out_df.index=out_df_index
	
	out_df.to_pickle(os.path.join(path,'profile_truncated.df'))
	
def process_data():

	# Set parameters
	path='C:\Offline_data_processing\Data\Whroo\Profile data'
	name='profile_truncated.df'
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
	