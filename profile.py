# Standard modules
import Tkinter, tkFileDialog
from configobj import ConfigObj
import os
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
def get_configs():
    root = Tkinter.Tk(); root.withdraw()
    cfName = tkFileDialog.askopenfilename(initialdir='')
    root.destroy()
    cf=ConfigObj(cfName)

    # Build dictionaries of config file contents
    paths_dict={}        
    variables_dict=dict(cf['variables'])
    options_dict=dict(cf['options'])

    # Set input file and output path and create directories for plots and results
    paths_dict['file_in']=os.path.join(cf['files']['input_path'],cf['files']['input_file'])

    results_output_path=os.path.join(cf['files']['output_path'],'Results')
    paths_dict['results_output_path']=results_output_path
    if not os.path.isdir(results_output_path): os.makedirs(results_output_path)
    if cf['options']['output_plots']:
        plot_output_path=os.path.join(cf['files']['output_path'],'Plots')
        paths_dict['plot_output_path']=plot_output_path
        if not os.path.isdir(plot_output_path): os.makedirs(plot_output_path)

    # Prepare dictionary of user settings - drop strings or change to int / float
    for key in options_dict:
        if options_dict[key].isdigit():
            options_dict[key]=int(options_dict[key])
        else:
            try:
                options_dict[key]=float(options_dict[key])
            except ValueError:
                continue
    
    return paths_dict,variables_dict,options_dict

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_data(file_in,variables_dict,options_dict):
    
    # Get data from spreadsheets and put into generic dict
    d={}
    sheets=[x for x in variables_dict.keys() if x in ['profile','met']]

    # 1) Check that the sheet for the profile data is present
    if not 'profile' in sheets:
        print 'The """profile""" subitem is missing from the """variables""" item in the configuration file'
        return
        
    # Get the data        
    for sheet in sheets:
        sheet_name=variables_dict[sheet]['sheet_name']
        DateTime_name=variables_dict[sheet]['DateTime_name']
        d[sheet]=pd.read_excel(file_in,sheetname=sheet_name)
        d[sheet].index=pd.to_datetime(d[sheet][DateTime_name])
        d[sheet].drop(DateTime_name,axis=1,inplace=True)
        freq=str(variables_dict[sheet]['freq_mins'])+'T'
        new_index=pd.date_range(start=d[sheet].index[0],end=d[sheet].index[-1],freq=freq)
        d[sheet]=d[sheet].reindex(new_index)

    # 2) Check that the c_vars item is present
    if 'C_vars' in variables_dict['profile']:
        C_vars=variables_dict['profile']['C_vars'][1:-1].split(',')
    else:
        print 'The """C_vars""" subitem is missing from the """profile""" item in the configuration file'
        return
    
    # 3) Check there is something in the item
    if len(C_vars)==0:
        print 'No variable names specifying CO2 profile levels'
        return
 
    # 4) Check the named variables are in the dataframe
    if not all(i in d['profile'].columns for i in C_vars):
        print 'Variable names specifying CO2 levels not found in spreadsheet'
        return

    # Flag whether met is on separate sheet
    if 'met' in sheets:
        met_sheet='profile' if variables_dict['profile']['sheet_name']==variables_dict['met']['sheet_name'] else 'met'
    else:
        met_sheet='profile'
    
    # Check that specified variables are present, of appropriate number (same as profile or 1), and appear in spreadsheet
    met_vars=[]
    defaults_list=[]
    for met_var in ['T_vars','p_vars']:
        if not met_var in variables_dict[met_sheet]:
            defaults_list.append(met_var)
        else:
            temp_list=variables_dict[met_sheet][met_var][1:-1].split(',')
            if len(temp_list)!=len(C_vars) and len(temp_list)!=1:
                defaults_list.append(met_var)
            else:
                if not all(i in d[met_sheet].columns for i in temp_list):
                    defaults_list.append(met_var)
                else:
                    met_vars=met_vars+temp_list

    # If there is a separate met sheet but no valid data, can the sheet and switch the met sheet to 'profile'
    if met_sheet=='met' and not met_vars:
        d.pop('met')
        met_sheet='profile'
        sheets.remove('met')
        
    # Subset dataframe to remove all extraneous data before return
    for sheet in sheets:    
        if sheet=='profile':
            if met_sheet=='met':
                all_vars=C_vars
            else:
                all_vars=C_vars+met_vars
        else:
            if met_sheet=='met':            
                all_vars=met_vars
        d[sheet]=d[sheet][all_vars]
    
    # Insert defaults
    if defaults_list:
        if 'T_vars' in defaults_list:
            d[met_sheet]['T_default']=20
            variables_dict[met_sheet]['T_vars']='[T_default]'
        if 'p_vars' in defaults_list:
            d[met_sheet]['p_default']=101.3
            variables_dict[met_sheet]['p_vars']='[p_default]'
    
    # Drop the 'met' item from the variables dictionary if present
    if met_sheet=='profile':
        if 'met' in variables_dict.keys():
            variables_dict.pop('met')

    return d

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_layers(levels_dict):

    # Assemble dictionaries containing variable name:layer tuples
    levels=levels_dict['Cc'].keys()
    levels.sort()
    layers=[levels[0]]+[levels[i+1]-levels[i] for i in xrange(len(levels)-1)]
    layers_dict={}
    for var in levels_dict.keys():
        if len(set(levels_dict[var].values()))==len(levels_dict[var].values()):
            layer_names=[]
            for j in xrange(len(layers)):            
                if j==0:
                    layer_names.append(var+'_'+'0-'+str(levels[j])+'m')
                else:    
                    layer_names.append(var+'_'+str(levels[j-1])+'-'+str(levels[j])+'m')
            layers_dict[var]=dict(zip(layers,layer_names))
        else:
            layer_names=[levels_dict[var][levels[0]]]*len(levels)
            layers_dict[var]=dict(zip(layers,layer_names))    
    
    return layers_dict

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def get_levels(variables_dict,levels):

    # Assemble dictionaries containing variable name:level tuples
    C_vars=variables_dict['profile']['C_vars']
    C_vars=[i for i in C_vars[1:-1].split(',')]

    # Set variables containing meteorological data
    var='met' if 'met' in variables_dict.keys() else 'profile'

    # Temperature    
    T_vars=variables_dict[var]['T_vars']
    T_vars=[i for i in T_vars[1:-1].split(',')]
    if len(T_vars)==1:
        T_vars=T_vars*len(C_vars)

    # Pressure
    p_vars=variables_dict[var]['p_vars']
    p_vars=[i for i in p_vars[1:-1].split(',')]
    if len(p_vars)==1:
        p_vars=p_vars*len(C_vars)
      
    # Get heights
    if not levels==None:
        levels=[float(i) for i in levels[1:-1].split(',')]
        print 'Using levels supplied'
    else:
        levels=[float(re.findall(r'[-+]?\d*\.\d+|\d+',i)[0]) for i in C_vars]
        print 'Getting levels from variable names'
    
    # Create dictionaries with linked variable names and levels
    C_dict=dict(zip(levels,C_vars))
    temp_dict=dict(zip(levels,T_vars))
    p_dict=dict(zip(levels,p_vars))
    
    return {'Cc':C_dict,'Ta':temp_dict,'p':p_dict}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
def main():

    # Read configurations into dicts
    print 'Select control file containing run parameters'
    paths_dict,variables_dict,options_dict=get_configs()
      
    # Get data
    df_dict=get_data(paths_dict['file_in'],variables_dict,options_dict)
    if not df_dict:
        print 'Exiting'
        return

    # Build frequencies dictionary
    freq_dict={'profile':int(variables_dict['profile']['freq_mins']),
               'output':options_dict['output_freq_mins']}
    if 'met' in variables_dict.keys():
        freq_dict['met']=int(variables_dict['met']['freq_mins'])
    else:
        freq_dict['met']=int(variables_dict['profile']['freq_mins'])

    # Set profile heights
    levels=options_dict['levels'] if 'levels' in options_dict.keys() else None
    levels_dict=get_levels(variables_dict,levels)
    if not levels_dict:
        print 'Exiting'
        return
    
    # Set profile layers
    layers_dict=get_layers(levels_dict)

    # Truncate profile data to required frequency
    # Why do I need to return profile_df - why doesn't the operation in the function change
    # the df, since the passed argument is a reference to the dataframe?
    df_dict['profile']=truncate_data(df_dict['profile'],freq_dict,options_dict)
    if not isinstance(df_dict['profile'],pd.DataFrame):
        print 'Exiting'
        return
   
    # If separate meteorology, resample, align, fill as necessary then join
    if 'met' in df_dict.keys():
        df_dict['met']=resample_data(df_dict['met'],freq_dict)
        profile_df=merge_data(df_dict)
    else:
        profile_df=df_dict['profile']

    # Calculate storage
    storage_df=process_data(profile_df,levels_dict,layers_dict,freq_dict,paths_dict)

    # Return storage dataframe
    return storage_df

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def merge_data(df_dict):
    
    # Join dataframes - reindex meteorology to match CO2 profile then fill
    print 'Merging CO2 profile and meteorological data'
    profile_df=df_dict['profile']
    met_df=df_dict['met'].reindex(profile_df.index).fillna(method='bfill').fillna(method='ffill')
    profile_df=profile_df.join(met_df)
    print 'Done'
    
    return profile_df

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def process_data(df,levels_dict,layers_dict,freq_dict,paths_dict):

    # Set parameters
    r=8.3143 # Universal gas constant
    K_con=273.15 # Kelvin conversion
    
    # Get numeric levels and layers
    levels=levels_dict['Cc'].keys()
    levels.sort()
    layers=layers_dict['Cc'].keys()
    layers.sort()

    # For variables measured at all heights, create dictionaries and calculate averages for each layer
    # (assume value at lowest measurement height is average for that layer (i.e. 0-0.5m));
    # For single height or default variables, do nothing!
    for var in layers_dict.keys(): 
        if len(set(levels_dict[var].values()))==len(levels_dict[var].values()):
            for i in xrange(len(levels)):
                level_var=levels_dict[var][levels[i]]
                layer_var=layers_dict[var][layers[i]]
                if i==0:
                    df[layer_var]=df[level_var]
                else:
                    prev_level_var=levels_dict[var][levels[i-1]]
                    df[layer_var]=(df[level_var]+df[prev_level_var])/2
    
    # Create an intermediate df for calculations
#    CO2_dels_list=['del_'+str(i) for i in layers_dict['Cc'].values()]
    C_molar_dens_list=['C_md_'+str(i)+'m' for i in layers]
    C_layer_stor_list=['C_stor'+i[2:] for i in layers_dict['Cc'].values()]
    T_list=list(set(layers_dict['Ta'].values()))
    p_list=list(set(layers_dict['p'].values()))

    # Here if meteorological data recorded at slower interval, we use the mean from the two
    # periods neighbouring the CO2 measurement - for example, if the CO2 measurement is centred
    # on 1230, the mean of 1230 and 1300 temperature data (where 1230 represents the average for measurements
    # over the period 1200-1230) is the best estimate of the quasi-instantaneous 1230 temperature
    if freq_dict['met']>freq_dict['profile']:
        for i in xrange(len(T_list)):
            target_var=T_list[i]
            df[T_list[i]]=(df[target_var]+df[target_var].shift(-1))/2
            df[T_list[i]].fillna(method='ffill',inplace=True)
        for i in xrange(len(p_list)):
            target_var=p_list[i]
            df[p_list[i]]=(df[target_var]+df[target_var].shift(-1))/2
            df[p_list[i]].fillna(method='ffill',inplace=True)

    # Do the storage calculation for each layer, scaling to umol over layer thickness
    for i in xrange(len(layers)):
        source_var_Cc=layers_dict['Cc'][layers[i]]
        source_var_Ta=layers_dict['Ta'][layers[i]]
        source_var_p=layers_dict['p'][layers[i]]
        target_var=C_molar_dens_list[i]
        pdb.set_trace()
        df[C_molar_dens_list[i]]=(df[source_var_p]*10**3/(r*(K_con+df[source_var_Ta]))* # Gas molar density
                                  df[source_var_Cc]* # C molar density
                                  layers[i]) # Volume-integrated C molar density
        
    # Difference the time series and scale to 1s, then sum
    for i in xrange(len(layers)):
        df[C_layer_stor_list[i]]=(df[C_molar_dens_list[i]]-df[C_molar_dens_list[i]].shift())/(freq_dict['output']*60)
    df['C_stor_total']=df[C_layer_stor_list].sum(axis=1)
    vars_out_list=['C_stor_total']+C_layer_stor_list+T_list+p_list
 
    return df[vars_out_list]
        
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def resample_data(met_df,freq_dict):
    
    # Get frequencies from dictionary
    output_freq=freq_dict['output']
    met_freq=freq_dict['met']
    profile_freq=freq_dict['profile']
    
    # Resample if frequencies don't match
    if met_freq!=output_freq and met_freq!=profile_freq:
        print 'Resampling meteorological data to match CO2 profile data...'
        resample_freq=str(output_freq)+'T'
        met_df=met_df.resample(resample_freq)
        if met_freq>output_freq:
            met_df.interpolate(inplace=True)
        print 'Done'

    return met_df

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
def truncate_data(df,freq_dict,options_dict):
    
    output_freq=freq_dict['output']
    profile_freq=freq_dict['profile']
    lag=options_dict['lag']
    bin_size=options_dict['bin_size']
    
    if output_freq==profile_freq:
        print 'Output and profile frequencies match... skipping upsampling'
        return
    elif output_freq<profile_freq:
        print 'Profile data frequency is lower than requested output frequency! Exiting...'
        return
    
    print 'Resampling all data to '+str(output_freq)+'minute frequency:'

    # Find the time window to be used given the lag and bin size settings above
    temp_df=pd.DataFrame(index=df.index)
    temp_df['mins']=temp_df.index
    temp_df['mins']=temp_df['mins'].apply(lambda x: x.timetuple().tm_min)
    temp_df['mins']=temp_df['mins']%output_freq # Rewrite minutes to half hourly repetition if desired frequency is 30min
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
    step=(output_freq/2)
    for i in xrange(0,len(df),step):
        df['grp'].iloc[i:i+bin_size]=i

    # Create the output df and do the averaging
    df_index=pd.to_datetime(temp_df[temp_df['valid_tmstmp']].index)
    df=df.groupby('grp').mean()
    df.index=df_index
    
    print 'Done'
    
    return df
	
#------------------------------------------------------------------------------
	
