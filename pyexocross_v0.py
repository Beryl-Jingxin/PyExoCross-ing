
# Import all what we need.
# encoding: utf-8
import os
import bz2
import glob
import time
import requests
import argparse
import numpy as np
import pandas as pd
import numexpr as ne
from tqdm import tqdm
from io import StringIO
#import astropy.units as au
import matplotlib.pyplot as plt
#import threading, multiprocessing
from scipy.special import voigt_profile, wofz, erf

# The input file path
def parse_args():
    parse = argparse.ArgumentParser(description='PyExoCross Program')
    parse.add_argument('-p', '--path', type=str, metavar='', required=True, help='Input file path')
    args = parse.parse_args()
    inp_filepath = args.path
    return inp_filepath

# Report time
class Timer:    
    def start(self):
        self.start_CPU = time.process_time()
        self.start_sys = time.time()
        return self

    def end(self, *args):
        self.end_CPU = time.process_time()
        self.end_sys = time.time()
        self.interval_CPU = self.end_CPU - self.start_CPU
        self.interval_sys = self.end_sys - self.start_sys
        print('{:25s} : {}'.format('Running time on CPU', self.interval_CPU), 's')
        print('{:25s} : {}'.format('Running time on system', self.interval_sys), 's')

# Read information from input file

def inp_para(inp_filepath):
    
    # Find the maximum column for all the rows.
    with open(inp_filepath, 'r') as temp_f:
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1).
    column_names = [i for i in range(max(col_count)+1)] 
    inp_df = pd.read_csv(inp_filepath, sep='\\s+', header = None, names=column_names, usecols=column_names)
    col0 = inp_df[0]
    
    # Database
    database = inp_df[col0.isin(['Database'])][1].values[0]
    
    # Basic information
    molecule = inp_df[col0.isin(['Molecule'])][1].values[0]
    isotopologue = inp_df[col0.isin(['Isotopologue'])][1].values
    dataset = inp_df[col0.isin(['Dataset'])][1].values
    mol_iso_id = int(inp_df[col0.isin(['mol_iso_id'])][1])
    
    # File path
    read_path = inp_df[col0.isin(['ReadPath'])][1].values[0]
    save_path = inp_df[col0.isin(['SavePath'])][1].values[0]
    if os.path.exists(save_path):
        pass
    else:
        os.makedirs(save_path, exist_ok=True)
        
    # Functions 
    Conversion = int(inp_df[col0.isin(['Conversion'])][1])
    PartitionFunctions = int(inp_df[col0.isin(['PartitionFunctions'])][1])
    CoolingFunctions = int(inp_df[col0.isin(['CoolingFunctions'])][1])
    Lifetimes = int(inp_df[col0.isin(['Lifetimes'])][1])
    SpecificHeats = int(inp_df[col0.isin(['SpecificHeats'])][1])
    StickSpectra = int(inp_df[col0.isin(['StickSpectra'])][1])
    CrossSections = int(inp_df[col0.isin(['CrossSections'])][1])
    
    # Convert from one format to another
    if Conversion == 1:
        ConversionFormat = int(inp_df[col0.isin(['ConversionFormat'])][1])
        ConversionMinFreq = float(inp_df[col0.isin(['ConversionFrequncyRange'])][1])
        ConversionMaxFreq = float(inp_df[col0.isin(['ConversionFrequncyRange'])][2])
        ConversionUnc = float(inp_df[col0.isin(['ConversionUncertainty'])][1])
        GlobalQNLabel_list = list(inp_df[col0.isin(['GlobalQNLabel'])].iloc[0].dropna())[1:]
        GlobalQNFormat_list = list(inp_df[col0.isin(['GlobalQNFormat'])].iloc[0].dropna())[1:]
        LocalQNLabel_list = list(inp_df[col0.isin(['LocalQNLabel'])].iloc[0].dropna())[1:]
        LocalQNFormat_list = list(inp_df[col0.isin(['LocalQNFormat'])].iloc[0].dropna())[1:]
    else:
        ConversionFormat = 0
        ConversionMinFreq = 0
        ConversionMaxFreq = 0
        ConversionUnc = 100
        GlobalQNLabel_list = []
        GlobalQNFormat_list = []
        LocalQNLabel_list = []
        LocalQNFormat_list = []    
    
    # Calculate partition, cooling functions or specific heats 
    Ntemp = int(inp_df[col0.isin(['Ntemp'])][1])    # The number of temperature steps
    Tmax = int(inp_df[col0.isin(['Tmax'])][1])      # Maximal temperature in K (minimal T = 1 K )
     
    # Calculate lifetimes 
    # None 
    
    # Calculate stick spectra or cross-sections 
    broadeners = list(inp_df[col0.isin(['Broadeners'])].iloc[0])[1:]
    ratios = np.array(list(inp_df[col0.isin(['Ratios'])].iloc[0])[1:], dtype=float)
    
    T = int(inp_df[col0.isin(['Temperature'])][1])
    P = float(inp_df[col0.isin(['Pressure'])][1])
    min_wn = float(inp_df[col0.isin(['Range'])][1])
    max_wn = float(inp_df[col0.isin(['Range'])][2])
    NpointsORBinSize = inp_df[col0.isin(['Npoints/BinSize'])][1].values[0]
    if NpointsORBinSize in ['Npoints', 'Npoint', 'NPoints', 'NPoint', 'npoints', 'npoint']:
        N_point = int(inp_df[col0.isin(['Npoints/BinSize'])][2])
        bin_size = float((max_wn - min_wn)/(N_point-1))
    elif NpointsORBinSize in ['BinSize', 'Binsize', 'binsize', 'binSize', 'bin', 'size', 'Bin', 'Size']:
        bin_size = float(inp_df[col0.isin(['Npoints/BinSize'])][2])
        N_point = int((max_wn - min_wn)/bin_size+1)
    else:
        raise ImportError("Please type the correct grid choice 'Npoints' or 'BinSize' into the input file.")
        
    wn_grid = np.linspace(min_wn, max_wn, N_point)
    
    # cutoffYN = float(inp_df[col0.isin(['Cutoff'])]['c1'])
    cutoffYN = inp_df[col0.isin(['Cutoff(Y/N)'])][1].values[0]
    if cutoffYN in ['Y', 'Yes', 'yes', 'YES']:
        cutoff = float(inp_df[col0.isin(['Cutoff(Y/N)'])][2])
    elif cutoffYN in ['N', 'No', 'no', 'NO']:
        cutoff = 'None'
    else:
        raise ImportError("Please type the correct cutoff choice 'Y' or 'N' into the input file.")
    
    thresholdYN = inp_df[col0.isin(['Threshold(Y/N)'])][1].values[0]
    if thresholdYN in ['Y', 'Yes', 'yes', 'YES']:
        threshold = float(inp_df[col0.isin(['Threshold(Y/N)'])][2])
    elif thresholdYN in ['N', 'No', 'no', 'NO']:
        threshold = 'None'
    else:
        raise ImportError("Please type the correct threshold choice 'Y' or 'N' into the input file.")
    
    UncFilterYN = inp_df[col0.isin(['UncFilter(Y/N)'])][1].values[0]
    if UncFilterYN in ['Y', 'Yes', 'yes', 'YES']:
        UncFilter = float(inp_df[col0.isin(['UncFilter(Y/N)'])][2])
    elif UncFilterYN in ['N', 'No', 'no', 'NO']:
        UncFilter = 'None'
    else:
        raise ImportError("Please type the correct uncertainty filter choice 'Y' or 'N' into the input file.")
    
    QNsFilterYN = inp_df[col0.isin(['QNsFilter(Y/N)'])][1].values[0]
    QNslabel_list = list(inp_df[col0.isin(['QNslabel'])].iloc[0])[1:]
    QNsformat_list = list(inp_df[col0.isin(['QNsformat'])].iloc[0])[1:]
    if QNsFilterYN in ['Y', 'Yes', 'yes', 'YES']:
        QNsFilter = list(inp_df[col0.isin(['QNsFilter(Y/N)'])].iloc[0].dropna())[2:]
        QNs_label = []
        QNs_value = []
        for i in range(len(QNsFilter)):
            QNs_label.append(QNsFilter[i].split('[')[0])
            QNs_value.append(QNsFilter[i].split('[')[1].split(']')[0].split(','))
    elif QNsFilterYN in ['N', 'No', 'no', 'NO']:
        QNsFilter = []
        QNs_label = []
        QNs_value = []
    else:
        raise ImportError("Please type the correct quantum number filter choice 'Y' or 'N' into the input file.")

    abs_emi = inp_df[col0.isin(['Absorption/Emission'])][1].values[0]
    profile = inp_df[col0.isin(['Profile'])][1].values[0]
    wn_wl = inp_df[col0.isin(['Wavenumber(wn)/wavelength(wl)'])][1].values[0]   

    molecule_id = int(mol_iso_id/10)
    isotopologue_id = mol_iso_id - molecule_id * 10
        
    if database == 'ExoMol':
        # Read ExoMol definition file (.def) to get the mass.
        deffile_path = (read_path+'/'+molecule+'/'+isotopologue+'/'+dataset+'/'+isotopologue+'__'+dataset+'.def')
        def_df = pd.read_csv(deffile_path[0],sep='\\s+',usecols=[0,1,2,3,4],names=['0','1','2','3','4'],header=None)
        abundance = 1
        mass = float(def_df[def_df['4'].isin(['mass'])]['0'].values[0])     # ExoMol mass (Dalton)
        if def_df.to_string().find('Uncertainty') != -1:
            check_uncertainty = int(def_df[def_df['2'].isin(['Uncertainty'])]['0'].values[0])
        else:
            check_uncertainty = 0
        check_lifetime = int(def_df[def_df['2'].isin(['Lifetime'])]['0'].values[0])
        check_gfactor = int(def_df[def_df['3'].isin(['g-factor'])]['0'].values[0])
    elif database == 'HITRAN':
        isometa_url = 'https://hitran.org/docs/iso-meta/'
        iso_meta_table = pd.read_html(isometa_url)[molecule_id - 1]
        iso_meta_row = iso_meta_table[iso_meta_table['local ID'].isin([isotopologue_id])]
        abundance = float(iso_meta_row['Abundance'][0].replace('\xa0×\xa010','E'))
        mass = float(iso_meta_row['Molar Mass /g·mol-1'])                   # HITRAN molar mass (g/mol)
        check_uncertainty = 0
        check_lifetime = 0
        check_gfactor = 0
    else:
        raise ImportError("Please add the name of the database 'ExoMol' or 'HITRAN' into the input file.")
    
    return (database, molecule, isotopologue, dataset, read_path, save_path, 
            Conversion, PartitionFunctions, CoolingFunctions, Lifetimes, SpecificHeats, StickSpectra, CrossSections,
            ConversionFormat, ConversionMinFreq, ConversionMaxFreq, ConversionUnc, 
            GlobalQNLabel_list, GlobalQNFormat_list, LocalQNLabel_list, LocalQNFormat_list,
            Ntemp, Tmax, broadeners, ratios, T, P, min_wn, max_wn, N_point, bin_size, wn_grid, 
            cutoff, threshold, UncFilter, QNslabel_list, QNsformat_list, QNs_label, QNs_value, QNsFilter, 
            abs_emi, profile, wn_wl, molecule_id, isotopologue_id, abundance, mass,
            check_uncertainty, check_lifetime, check_gfactor)


# Constants and parameters
# Parameters for calculating
# Parameters for calculating
import astropy.constants as ac
#from astropy import constants, units as ac, au
Tref = 296.0                        # Reference temperature is 296 K
Pref = 1.0                          # Reference pressure is 1 bar
N_A = ac.N_A.value                  # Avogadro number (1/mol)
h = ac.h.to('erg s').value          # Planck's const (erg s)
c = ac.c.to('cm/s').value           # Velocity of light (cm/s)
kB = ac.k_B.to('erg/K').value       # Boltzmann's const (erg/K)
R = ac.R.to('J / (K mol)').value    # Molar gas constant (J/(K mol))
c2 = h * c / kB                     # Second radiation constant (cm K)

inp_filepath = parse_args()
(database, molecule, isotopologue, dataset, read_path, save_path, 
 Conversion, PartitionFunctions, CoolingFunctions, Lifetimes, SpecificHeats, StickSpectra, CrossSections,
 ConversionFormat, ConversionMinFreq, ConversionMaxFreq, ConversionUnc, 
 GlobalQNLabel_list, GlobalQNFormat_list, LocalQNLabel_list, LocalQNFormat_list,
 Ntemp, Tmax, broadeners, ratios, T, P, min_wn, max_wn, N_point, bin_size, wn_grid, 
 cutoff, threshold, UncFilter, QNslabel_list, QNsformat_list, QNs_label, QNs_value, QNsFilter, 
 abs_emi, profile, wn_wl, molecule_id, isotopologue_id, abundance, mass, 
 check_uncertainty, check_lifetime, check_gfactor) = inp_para(inp_filepath)


c2InvTref = c2 / Tref                 # c2 / T_ref (cm)
PI = np.pi
sinPI = np.sin(np.pi)
SqrtPI = np.sqrt(np.pi)
Sqrtln2 = np.sqrt(np.log(2))
OneminSqrtPIln2 = 1 - np.sqrt(np.pi * np.log(2))
Negln2 = -np.log(2)
Inv8Pic = 1 / (8 * np.pi * c)         # 8 * pi * c (s/cm)
Inv4Pi = 1 / (4 * np.pi)
InvSqrt2 = 1 / np.sqrt(2)
InvSqrt2Pi = 1 / np.sqrt(2 * np.pi)
InvSqrt2ln2 = 1 / np.sqrt(2 * np.log(2))
TwoSqrt2ln2 = 2 * np.sqrt(2 * np.log(2))
Sqrtln2InvPi = np.sqrt(np.log(2) / np.pi)
Sqrt2NAkBln2mInvc = np.sqrt(2 * N_A * kB * np.log(2) / mass) / c
binSize2 = bin_size * 2
binSizePI = bin_size * np.pi
binSizePI32 = bin_size * np.pi**1.5
binSizeHalf = bin_size / 2 
# Sqrt2NAkBln2mInvc = np.sqrt(2 * kB * np.log(2) / (mass*au.u.to(au.g))) / c


# Read input files
''' Read the parameters of the linelist in ExoMol or HITRAN format text file. 
    Return the dataframe of the data for the following calculations.'''
# Read ExoMol database files
# Read states file
def read_all_states(read_path):
    s_df = dict()
    states_df = pd.DataFrame()
    states_filenames = glob.glob((read_path + molecule + '/' + isotopologue + '/' + dataset 
                                  + '/' + isotopologue + '__' + dataset + '.states.bz2')[0])
    for states_filename in states_filenames:
        s_df[states_filename] = pd.read_csv(states_filename, compression='bz2', sep='\s+', header=None,
                                            chunksize=100000, iterator=True, low_memory=False, dtype=object)
        for chunk in s_df[states_filename]:
            states_df = pd.concat([states_df, chunk])                         
    return(states_df)

# Read transitions file
def get_transfiles(read_path):
    # Get all the transitions files from the folder including the older version files which are named by vn(version number).
    trans_filepaths_all = glob.glob((read_path + molecule + '/' + isotopologue + '/' + dataset + '/' + '*trans.bz2')[0])
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []    # The list of the lastest transitions files.
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0].split('_')    # Split the filenames.
        num = len(split_version)
        # There are four format filenames.
        # The lastest transitions files named in two formats:
        # 1. Filenames are named with the name of isotopologue and dataset. 
        #    End with .trans.bz2.
        #    e.g. 14N-16O__XABC.trans.bz2'
        # 2. Filenames are named with the name of isotopologue and dataset. 
        #    Also have the range of wavenumbers xxxxx-yyyyy.
        #    End with .trans.bz2.
        #    e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
        # 3. The older version transitions files are named with vn(version number) based on the first format of the lastest files.
        #    e.g. 14N-16O__XABC_v2.trans.bz2
        # 4. The older version transitions files are named with updated date (yyyymmdd).
        #    e.g. 1H3_p__MiZATeP__20170330.trans.bz2
        # After split the filenames:
        # The first format filenames only leave the dataset name, e.g. XABC.
        # The second format filenames only leave the range of the wavenumber, e.g. 00000-00100.
        # The third format filenames leave two parts(dataset name and version number), e.g. XABC and v2.
        # The fourth format filenames only leave the updated date, e.g. 20170330.
        # This program only process the lastest data, so extract the filenames named by the first two format.
        if num == 1:     
            if split_version[0] == dataset:        
                trans_filepaths.append(trans_filepaths_all[i])
            if len(split_version[0].split('-')) == 2:
                trans_filepaths.append(trans_filepaths_all[i])
    return(trans_filepaths)    

def read_all_trans(read_path):
    t_df = dict()
    trans_df = pd.DataFrame()
    trans_col_name = ['u', 'l', 'A', 'v']
    trans_filepaths = get_transfiles(read_path)
    for trans_filename in tqdm(trans_filepaths):
        t_df[trans_filename] = pd.read_csv(trans_filename, compression='bz2', sep='\s+', header=None,
                                           names=trans_col_name, chunksize=100000, iterator=True, low_memory=False)
        for chunk in t_df[trans_filename]:
            trans_df = pd.concat([trans_df,chunk])                   
    return(trans_df)

    
# Convert among the frequency, upper and lower state energy 
# Calculate frequency
def cal_v(Ep, Epp):
    v = ne.evaluate('abs(Ep - Epp)')
    return(v)

# Calculate upper state energy
def cal_Ep(Epp, v):
    Ep = ne.evaluate('abs(Epp + v)')
    return(Ep)

# Read partition function file from ExoMol database
# Read partition function with online webpage.
def read_exomolweb_pf(T):
    pf_url = ('http://www.exomol.com/db/' + molecule + '/' + isotopologue + '/' + dataset 
              + '/' + isotopologue + '__' + dataset + '.pf')[0]
    pf_content = requests.get(pf_url).text
    pf_col_name = ['T', 'Q']
    pf_df = pd.read_csv(StringIO(pf_content), sep='\\s+', names=pf_col_name, header=None)
    Q = pf_df['Q'][T-1]
    return(Q)

# Read partition function with local partition function file.
def read_exomol_pf(read_path, T):
    pf_filename = (read_path + molecule + '/' + isotopologue + '/' + dataset 
                   + '/' + isotopologue + '__' + dataset + '.pf')[0]
    pf_col_name = ['T', 'Q']
    pf_df = pd.read_csv(pf_filename, sep='\\s+', names=pf_col_name, header=None)
    Q = pf_df['Q'][T-1]
    return(Q)

# Read broadening file
def read_broad(read_path):
    broad_col_name = ['code', 'gamma_L', 'n_air', 'Jpp']
    default_broad_df = pd.DataFrame(columns=broad_col_name)
    air_broad_df = pd.DataFrame(columns=broad_col_name)
    self_broad_df = pd.DataFrame(columns=broad_col_name)
    H2_broad_df = pd.DataFrame(columns=broad_col_name)
    He_broad_df = pd.DataFrame(columns=broad_col_name)
    CO2_broad_df = pd.DataFrame(columns=broad_col_name)

    broad = []
    ratio = []
    for i in range(len(ratios)):
        if ratios[i] != 0.0:
            if broadeners[i] == 'Default':
                default_gamma_L = 0.07
                default_n_air = 0.5
                default_broad_df = pd.DataFrame([['code', default_gamma_L, default_n_air,'Jpp']],columns=broad_col_name)
            if broadeners[i] == 'Air':
                broad_path = (read_path + molecule + '/' + isotopologue + '/' + isotopologue + '__air.broad')[0]
                if not os.path.exists(broad_path):
                    raise ImportError('The boradening file ' + broad_path + ' does not exist.')    
                air_broad_df = pd.read_csv(broad_path, sep='\s+', names=broad_col_name, header=None, engine='python')
            if broadeners[i] == 'Self':
                broad_path = (read_path + molecule + '/' + isotopologue + '/' + isotopologue + '__self.broad')[0]
                if not os.path.exists(broad_path):
                    raise ImportError('The boradening file ' + broad_path + ' does not exist.')    
                self_broad_df = pd.read_csv(broad_path, sep='\s+', names=broad_col_name, header=None, engine='python')
            if broadeners[i] == 'H2':
                broad_path = (read_path + molecule + '/' + isotopologue + '/' + isotopologue + '__H2.broad')[0]
                if not os.path.exists(broad_path):
                    raise ImportError('The boradening file ' + broad_path + ' does not exist.')    
                H2_broad_df = pd.read_csv(broad_path, sep='\s+', names=broad_col_name, header=None, engine='python')
            if broadeners[i] == 'He':
                broad_path = (read_path + molecule + '/' + isotopologue + '/' + isotopologue + '__He.broad')[0]
                if not os.path.exists(broad_path):
                    raise ImportError('The boradening file ' + broad_path + ' does not exist.')    
                He_broad_df = pd.read_csv(broad_path, sep='\s+', names=broad_col_name, header=None, engine='python')
            if broadeners[i] == 'CO2':
                broad_path = (read_path + molecule + '/' + isotopologue + '/' + isotopologue + '__CO2.broad')[0]
                if not os.path.exists(broad_path):
                    raise ImportError('The boradening file ' + broad_path + ' does not exist.')    
                CO2_broad_df = pd.read_csv(broad_path, sep='\s+', names=broad_col_name, header=None, engine='python')
            
            broad.append(broadeners[i])
            ratio.append(ratios[i])
    broad = list(i for i in broad if i==i)
    ratio = list(i for i in ratio if i==i)
    print('Broadeners \t: ', broadeners[0],'\t',broadeners[1],'\t',broadeners[2],'\t',broadeners[3],'\t',broadeners[4],'\t',broadeners[5])
    print('Ratios \t\t: ', ratios[0],'\t\t',ratios[1],'\t',ratios[2],'\t',ratios[3],'\t',ratios[4],'\t',ratios[5])
    return(broad,ratio,default_broad_df,air_broad_df,self_broad_df,H2_broad_df,He_broad_df,CO2_broad_df)        
   

# Read HITRAN database files
# Process HITRAN linelist data
def read_hitran_parfile (read_path, parfile_df):
    '''
    Read the parameters of the molecular absorption features
    of HITRAN2020 format text file.
    
    Parameters
    ----------
    par_filepath : str
        Input file path for reading.
    Return
    ------
    hitran_df : DataFrame
        The DataFrame of HITRAN data for the molecule.
    '''    
    par_filename = read_path.split('/')[-1]
    if (len(str(parfile_df[0][0])) < 160):
        raise ImportError('The file ' + par_filename + ' is not a HITRAN2020 format data file.')
    #hitran_column_name = ['M','I','v','S','Acoeff','gamma_air','gamma_self',
    #                     'Epp','n_air','delta_air','Vp','Vpp','Qp','Qpp',
    #                     'Ierr','Iref','flag','gp','gpp']

    hitran_df = pd.DataFrame()
    hitran_df['M'] = pd.to_numeric(parfile_df[0].map(lambda x: x[0:2]), errors='coerce').astype('int32')                 # Molecule identification number
    hitran_df['I'] = pd.to_numeric(parfile_df[0].map(lambda x: x[2:3]), errors='coerce').astype('int32')                 # Isotopologue number
    hitran_df['v'] = pd.to_numeric(parfile_df[0].map(lambda x: x[3:15]), errors='coerce').astype('float64')              # Transition wavenumber (in cm^{-1})
    hitran_df['S'] = pd.to_numeric(parfile_df[0].map(lambda x: x[15:25]), errors='coerce').astype('float64')             # Intensity (cm^{-1} / (molecule cm^{-2}))
    hitran_df['A'] = pd.to_numeric(parfile_df[0].map(lambda x: x[25:35]), errors='coerce').astype('float64')             # The Einstein-A coefficient (s^{-1}) of a transition
    hitran_df['gamma_air'] = pd.to_numeric(parfile_df[0].map(lambda x: x[35:40]), errors='coerce').astype('float64')     # Air-broadened half-width at half maximum (HWHM) coefficient (cm^{-1} atm^{-1})
    hitran_df['gamma_self'] = pd.to_numeric(parfile_df[0].map(lambda x: x[40:45]), errors='coerce').astype('float64')    # Self-broadened half-width at half maximum (HWHM) coefficient (cm^{-1} atm^{-1})
    hitran_df['Epp'] = pd.to_numeric(parfile_df[0].map(lambda x: x[45:55]), errors='coerce').astype('float64')           # Lower state energy (cm^{-1})
    hitran_df['n_air'] = pd.to_numeric(parfile_df[0].map(lambda x: x[55:59]), errors='coerce').astype('float64')         # Temperature-dependent exponent for gamma_air
    hitran_df['delta_air'] = pd.to_numeric(parfile_df[0].map(lambda x: x[59:67]), errors='coerce').astype('float64')     # Air pressure_include line shift (cm^{-1} atm^{-1})
    #hitran_df['Vp'] = parfile_df[0].map(lambda x: x[67:82])                                                              # Upper-state "global" quanta
    #hitran_df['Vpp'] = parfile_df[0].map(lambda x: x[82:97])                                                             # Lower-state "global" quanta
    #hitran_df['Qp'] = parfile_df[0].map(lambda x: x[97:112])                                                             # Upper-state "local" quanta
    #hitran_df['Qpp'] = parfile_df[0].map(lambda x: x[112:127])                                                           # Lower-state "local" quanta
    hitran_df['gp'] = pd.to_numeric(parfile_df[0].map(lambda x: x[146:153]), errors='coerce').astype('float64')          # Statistical weight of the upper state
    
    hitran_df = hitran_df[hitran_df['M'].isin([molecule_id])]
    hitran_df = hitran_df[hitran_df['I'].isin([isotopologue_id])]
    hitran_df = hitran_df[hitran_df['v'].between(min_wn, max_wn)]
    if threshold != 'None':
        hitran_df = hitran_df[hitran_df['S'] >= threshold]
    return hitran_df


# Read HITRAN linelist file
def read_parfile(read_path):
    if not os.path.exists(read_path):
        raise ImportError('The input file ' + read_path + ' does not exist.')
    # Initialise the iterator object.
    read_par = pd.read_csv(read_path, chunksize=10000, iterator=True, header=None, encoding='utf-8')
    par_df = pd.DataFrame()
    for chunk in read_par:
        par_df = pd.concat([par_df, chunk])
    return(par_df)

# Read partition function file from HITRANOnline
def read_hitran_pf(T):
    isometa_url = 'https://hitran.org/docs/iso-meta/'
    iso_meta_table = pd.read_html(isometa_url)[molecule_id - 1]
    iso_meta_row = iso_meta_table[iso_meta_table['local ID'].isin([isotopologue_id])]
    Q_url = 'https://hitran.org/data/Q/' + iso_meta_row.loc[0][7]
    Q_content = requests.get(Q_url).text
    Q_col_name = ['T', 'Q']
    Q_df = pd.read_csv(StringIO(Q_content), sep='\\s+', names=Q_col_name, header=None)
    Q = Q_df['Q'][T - 1]   
    return(Q)

# Calculate parition function
def calculate_partition(En, gn, T):
    partition_func = ne.evaluate('sum(gn * exp(-c2 * En / T))') 
    return(partition_func)

# Partition function
def exomol_partition_func(states_df, Ntemp, Tmax):
    
    print('Calculate partition functions.')  
    t = Timer()
    t.start()

    En = states_df[1].astype('float').values
    gn = states_df[2].astype('int').values
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    
    partition_func = [calculate_partition(En, gn, T) for T in Ts]
    
    partition_func_df = pd.DataFrame()
    partition_func_df['T'] = Ts
    partition_func_df['partition function'] = partition_func
        
    pf_folder = save_path + '/partition/'
    if os.path.exists(pf_folder):
        pass
    else:
        os.makedirs(pf_folder, exist_ok=True)
    pf_path = (pf_folder + isotopologue + '__' + dataset + '.pf')[0]
    np.savetxt(pf_path, partition_func_df, fmt="%8.1f %15.4f")
    
    t.end()
    print('Partition functions has been saved!\n')  


# Specific heat
def calculate_specific_heats(En, gn, T):
    pf = ne.evaluate('sum(gn * exp(-c2 * En / T)) ')  
    pfp = ne.evaluate('sum(gn * np.exp(-c2 * En / T) * (c2 * En / T))')
    pfpp = ne.evaluate('sum(gn * np.exp(-c2 * En / T) * (c2 * En / T) ** 2)')
    specificheat_func = ne.evaluate('R * (pfpp / pf - (pfp / pf)**2) + 2.5 * R') 
    return(specificheat_func)

# Specific heat
def exomol_specificheat(states_df, Ntemp, Tmax):
    print('Calculate specific heats.')  
    t = Timer()
    t.start()
    
    En = states_df[1].astype('float').values
    gn = states_df[2].astype('int').values
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    
    specificheat_func = [calculate_specific_heats(En, gn, T) for T in Ts]
    
    specificheat_func_df = pd.DataFrame()
    specificheat_func_df['T'] = Ts
    specificheat_func_df['specific heat'] = specificheat_func
        
    cp_folder = save_path + '/specific_heat/'
    if os.path.exists(cp_folder):
        pass
    else:
        os.makedirs(cp_folder, exist_ok=True)  
    cp_path = (cp_folder + isotopologue + '__' + dataset + '.cp')[0]
    np.savetxt(cp_path, specificheat_func_df, fmt="%8.1f %15.4f")

    t.end()
    print('Specific heats has been saved!\n')  


# Lifetime
def exomol_lifetime(read_path, states_df, all_trans_df):
    
    print('Calculate lifetimes.')  
    t = Timer()
    t.start()
    
    sum_A = all_trans_df.groupby('u')['A'].sum()
    lifetime = ne.evaluate('1 / sum_A') 
    lt_df = lifetime.map('{: >12.4E}'.format).reset_index()
    lt_df.columns=['u','lt']
    uid = lt_df['u']
    add_u = pd.DataFrame()
    add_u['u'] = pd.concat([states_df[0].astype('int'), uid]).drop_duplicates(keep=False)
    add_u['lt'] = '         inf'
    lifetime_df = pd.concat([lt_df, add_u], ignore_index=True)
    lifetime_df.sort_values('u',inplace=True)
    
    states_filenames = glob.glob((read_path + molecule + '/' + isotopologue + '/' + dataset 
                                  + '/' + isotopologue + '__' + dataset + '.states.bz2')[0])
    s_df = pd.read_csv(states_filenames[0], compression='bz2', header=None, dtype=object)
    lifetime_list = list(lifetime_df['lt'])
    nrows = len(s_df)
    new_rows = []
    if check_uncertainty == 0:
        for i in range(nrows):
            new_rows.append(s_df[0][i][:41]+lifetime_list[i]+s_df[0][i][53:]+'\n')
    if check_uncertainty == 1:
        for i in range(nrows):
            new_rows.append(s_df[0][i][:53]+lifetime_list[i]+s_df[0][i][65:]+'\n')

    lf_folder = save_path + '/lifetime/'
    if os.path.exists(lf_folder):
        pass
    else:
        os.makedirs(lf_folder, exist_ok=True)  
        
    ####### bz2 ######
    lf_path = (lf_folder + isotopologue + '__' + dataset + '.states.bz2')[0]

    with bz2.open(lf_path, 'wt') as f:
        for i in range(nrows):
            f.write(new_rows[i])
        f.close
    ##################
    
    ##### states #####
    '''
    lf_path = (lf_folder + isotopologue + '__' + dataset + '.states')[0]

    with open(lf_path, 'wt') as f:
        for i in range(nrows):
            f.write(new_rows[i])
        f.close
    '''
    ##################

    t.end()
    print('Lifetimes has been saved!\n')   
    
    
# Calculate cooling function
def linelist_coolingfunc(states_df, all_trans_df):
    id_u = all_trans_df['u'].values
    id_l = all_trans_df['l'].values
    states_df[0] = pd.to_numeric(states_df[0])
    states_df.set_index([0], inplace=True, drop=False)
    id_s = states_df[0]
    all_trans_df.set_index(['u'], inplace=True, drop=False)
    id_us = list(set(id_u).intersection(set(id_s)))
    trans_us_df = all_trans_df.loc[id_us]
    id_l = trans_us_df['l'].values
    id_ls = list(set(id_l).intersection(set(id_s)))
    trans_us_df.set_index(['l'], inplace=True, drop=False)
    trans_s_df = trans_us_df.loc[id_ls]
    id_su = trans_s_df['u'].values
    id_sl = trans_s_df['l'].values
    states_u_df = states_df.loc[id_su]
    states_l_df = states_df.loc[id_sl]

    Ep = states_u_df[1].values.astype('float')
    gp = states_u_df[2].values.astype('int')
    A = trans_s_df['A'].values.astype('float')

    if pd.isna(all_trans_df['v']).iloc[0] == False:
        v = trans_s_df['v'].values.astype('float')
    else:
        Epp = states_l_df[1].astype('float') # Upper state energy
        v = cal_v(Ep, Epp) 
    return (A, v, Ep, gp)

def calculate_cooling(A, v, Ep, gp, T, Q):
    _sum = ne.evaluate('sum(A * h * c * v * gp * exp(-c2 * Ep / T))')  
    cooling_func = ne.evaluate('_sum / (4 * PI * Q)')
    return(cooling_func)

# Cooling function
def exomol_cooling_func(read_path, states_df, all_trans_df, Ntemp, Tmax):
    
    print('Calculate cooling functions.')  
    t = Timer()
    t.start()
    
    A, v, Ep, gp = linelist_coolingfunc(states_df, all_trans_df)
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    Qs = [read_exomol_pf(read_path, T) for T in Ts]
    
    cooling_func = [calculate_cooling(A, v, Ep, gp, T, Q) for T,Q in zip(Ts,Qs)]
    
    cooling_func_df = pd.DataFrame()
    cooling_func_df['T'] = Ts
    cooling_func_df['cooling function'] = cooling_func

    cf_folder = save_path + '/cooling/'
    if os.path.exists(cf_folder):
        pass
    else:
        os.makedirs(cf_folder, exist_ok=True)  
    cf_path = (cf_folder + isotopologue + '__' + dataset + '.cf')[0]    
    np.savetxt(cf_path, cooling_func_df, fmt="%8.1f %20.8E")
    
    t.end()
    print('Cooling functions has been saved!\n')     


# Process data
def read_part_states(states_df):
    if UncFilter != 'None':
        states_part_df = states_df[states_df[4].astype(float) <= UncFilter]
        states_part_df[0] = pd.to_numeric(states_part_df[0])
        states_part_df.set_index([0], inplace=True, drop=False)
    else:
        states_part_df = states_df
        states_part_df[0] = pd.to_numeric(states_part_df[0])
        states_part_df.set_index([0], inplace=True, drop=False)
    if check_uncertainty == 1:
        col_unc = ['unc']
    else:
        col_unc = []
    if check_lifetime == 1:
        col_lifetime = ['tau']
    else:
        col_lifetime = []
    if check_gfactor == 1:
        col_gfac = ['gfac']
    else:
        col_gfac = []
    colname = ['id','E','g','J'] + col_unc + col_lifetime + col_gfac + QNslabel_list
    states_part_df.drop(states_part_df.columns[len(colname):], axis=1, inplace=True)
    states_part_df.columns = colname
    if QNsFilter !=[]:    
        for i in range(len(QNs_label)):
            states_part_df = states_part_df[states_part_df[QNs_label[i]].isin(QNs_value[i])]
    return(states_part_df)

def get_part_transfiles(read_path):
    # Get all the transitions files from the folder including the older version files which are named by vn(version number).
    trans_filepaths_all = glob.glob((read_path + molecule + '/' + isotopologue + '/' + dataset + '/' + '*trans.bz2')[0])
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []    # The list of the lastest transitions files.
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0].split('_')    # Split the filenames.
        num = len(split_version)
        # There are four format filenames.
        # The lastest transitions files named in two formats:
        # 1. Filenames are named with the name of isotopologue and dataset. 
        #    End with .trans.bz2.
        #    e.g. 14N-16O__XABC.trans.bz2'
        # 2. Filenames are named with the name of isotopologue and dataset. 
        #    Also have the range of wavenumbers xxxxx-yyyyy.
        #    End with .trans.bz2.
        #    e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
        # 3. The older version transitions files are named with vn(version number) based on the first format of the lastest files.
        #    e.g. 14N-16O__XABC_v2.trans.bz2
        # 4. The older version transitions files are named with updated date (yyyymmdd).
        #    e.g. 1H3_p__MiZATeP__20170330.trans.bz2
        # After split the filenames:
        # The first format filenames only leave the dataset name, e.g. XABC.
        # The second format filenames only leave the range of the wavenumber, e.g. 00000-00100.
        # The third format filenames leave two parts(dataset name and version number), e.g. XABC and v2.
        # The fourth format filenames only leave the updated date, e.g. 20170330.
        # This program only process the lastest data, so extract the filenames named by the first two format.
        if num == 1:     
            if split_version[0] == dataset:        
                trans_filepaths.append(trans_filepaths_all[i])
            elif len(split_version[0].split('-')) == 2:
                trans_filepaths.append(trans_filepaths_all[i])
        
    if len(trans_filepaths) == 1:
        filenames = trans_filepaths
    else:
        filenames = []
        for trans_filename in tqdm(trans_filepaths):
            lower = int(trans_filename.split('__')[2].split('.')[0].split('-')[0])
            upper = int(trans_filename.split('__')[2].split('.')[0].split('-')[1]) 
            if (lower <= int(min_wn) <= upper):
                filenames.append(trans_filename)
            if (lower >= int(min_wn) and upper <= int(max_wn)):
                filenames.append(trans_filename)
            if (lower <= int(max_wn) <= upper):
                filenames.append(trans_filename)   
    return(filenames)     

def read_part_trans(read_path):
    trans_filenames = get_part_transfiles(read_path)
    t_df = dict()
    trans_part_df = pd.DataFrame()
    # Initialise the iterator object.
    trans_col_name = ['u', 'l', 'A', 'v']
    for trans_filename in tqdm(trans_filenames):
        t_df[trans_filename] = pd.read_csv(trans_filename, compression='bz2', sep='\s+', header=None, 
                                           names=trans_col_name, chunksize=10000, iterator=True, encoding='utf-8')
        for chunk in t_df[trans_filename]:
            trans_part_df = pd.concat([trans_part_df, chunk])
    return(trans_part_df)


# Intensity
def cal_abscoefs(v, gp, A, Epp, Q, abundance):
    abscoef = ne.evaluate('gp * A * exp(- c2 * Epp / T) * (1 - exp(- c2 * v / T)) * Inv8Pic / (v ** 2 * Q) * abundance')  
    return abscoef

def cal_emicoefs(v, gp, A, Ep, Q, abundance):
    emicoef = ne.evaluate('gp * A * v * exp(- c2 * Ep / T) * Inv4Pi / Q * abundance')
    return emicoef


# Uncertainty
def cal_uncertainty(unc_u, unc_l):
    unc = ne.evaluate('sqrt(unc_u ** 2 + unc_l ** 2)')
    return unc


# Stick spectra
def linelist(cutoff,states_part_df,trans_part_df):
    if pd.isna(trans_part_df['v']).values[0] == False:
        if cutoff == 'None':
            trans_part_df = trans_part_df[trans_part_df['v'].between(min_wn, max_wn)] 
        else:
            trans_part_df = trans_part_df[trans_part_df['v'].between(min_wn - cutoff, max_wn + cutoff)] 
        id_u = trans_part_df['u'].values
        id_s = states_part_df['id'].values

        trans_part_df.set_index(['u'], inplace=True, drop=False)
        id_us = list(set(id_u).intersection(set(id_s)))
        trans_us_df = trans_part_df.loc[id_us]

        id_l = trans_us_df['l'].values
        id_ls = list(set(id_l).intersection(set(id_s)))
        trans_us_df.set_index(['l'], inplace=True, drop=False)
        trans_s_df = trans_us_df.loc[id_ls]
        trans_s_df.sort_values(by=['v'], inplace=True)

        id_su = trans_s_df['u'].values
        id_sl = trans_s_df['l'].values
        states_u_df = states_part_df.loc[id_su]
        states_l_df = states_part_df.loc[id_sl]

        Ep = states_u_df['E'].values.astype('float')
        Epp = states_l_df['E'].values.astype('float')
        gp = states_u_df['g'].values.astype('int')
        Jp = states_u_df['J'].values.astype('float')
        Jpp = states_l_df['J'].values.astype('float')
        A = trans_s_df['A'].values.astype('float')
        v = trans_s_df['v'].values.astype('float')
        QNp = pd.DataFrame()
        QNpp = pd.DataFrame()
        for i in range(len(QNslabel_list)):
            QNp[QNslabel_list[i]+"'"] = states_u_df[QNslabel_list[i]].values
            QNpp[QNslabel_list[i]+'"'] = states_l_df[QNslabel_list[i]].values
        stick_qn_df = pd.concat([QNp,QNpp],axis='columns')
    else:
        id_u = trans_part_df['u'].values
        id_s = states_part_df['id'].values

        trans_part_df.set_index(['u'], inplace=True, drop=False)
        id_us = list(set(id_u).intersection(set(id_s)))
        trans_us_df = trans_part_df.loc[id_us]

        id_l = trans_us_df['l'].values
        id_ls = list(set(id_l).intersection(set(id_s)))
        trans_us_df.set_index(['l'], inplace=True, drop=False)
        trans_s_df = trans_us_df.loc[id_ls]

        id_su = trans_s_df['u'].values
        id_sl = trans_s_df['l'].values
        states_u_df = states_part_df.loc[id_su]
        states_l_df = states_part_df.loc[id_sl]

        trans_s_df['Ep'] = states_u_df['E'].values.astype('float')
        trans_s_df['Epp'] = states_l_df['E'].values.astype('float')
        trans_s_df['gp'] = states_u_df['g'].values.astype('int')
        trans_s_df['v'] = cal_v(trans_s_df['Ep'].values, trans_s_df['Epp'].values)
        if cutoff == 'None':
            trans_s_df = trans_s_df[trans_s_df['v'].between(min_wn, max_wn)]
        else:
            trans_s_df = trans_s_df[trans_s_df['v'].between(min_wn - cutoff, max_wn + cutoff)]
        trans_s_df.sort_values(by=['v'], inplace=True)
        
        Epp = trans_s_df['Epp'].values
        gp = trans_s_df['gp'].values
        A = trans_s_df['A'].values
        v = trans_s_df['v'].values
        
        id_sl = trans_s_df['l'].values
        states_l_df = states_part_df.loc[id_sl]
    return (A, v, Ep, Epp, gp, Jp, Jpp, stick_qn_df)

# Stick spectra
def exomol_stick_spectra(read_path, states_part_df, trans_part_df, T):
    
    print('Calculate stick spectra.')  
    t = Timer()
    t.start()

    A, v, Ep, Epp, gp, Jp, Jpp, stick_qn_df = linelist(cutoff,states_part_df,trans_part_df)
    Q = read_exomol_pf(read_path, T)
    I = cal_abscoefs(v, gp, A, Epp, Q, abundance)
    stick_st_dic = {'v':v, 'I':I, "J'":Jp, "E'":Ep, 'J"':Jpp, 'E"':Epp}
    stick_st_df = pd.DataFrame(stick_st_dic)
    stick_spectra_df = pd.concat([stick_st_df, stick_qn_df], axis='columns')
    
    if isinstance(stick_spectra_df['J"'][0],int) == 'True':
        J_format = '%7s'
    else:
        J_format = '%7.1f'
    QNs_format = str(QNsformat_list).replace("', '"," ").replace("['","").replace("']","").replace('d','s').replace('.1f','s')
    
    ss_folder = save_path + '/stick_spectra/stick/'
    if os.path.exists(ss_folder):
        pass
    else:
        os.makedirs(ss_folder, exist_ok=True)
    ss_path = (ss_folder + isotopologue + '__' + dataset + '.stick')[0]
    ss_colname = stick_spectra_df.columns
    fmt = '%12.8E  %12.8E '+J_format+' %12.4f '+J_format+' %12.4f '+QNs_format+' '+QNs_format
    np.savetxt(ss_path, stick_spectra_df, fmt=fmt, header='')
    
    # Plot cross sections and save it as .png.
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    parameters = {'axes.labelsize': 14, 
                  'legend.fontsize': 14,
                  'xtick.labelsize': 12,
                  'ytick.labelsize': 12}
    plt.rcParams.update(parameters)
    ss_plot_folder = save_path + '/stick_spectra/plots/'
    if os.path.exists(ss_plot_folder):
        pass
    else:
        os.makedirs(ss_plot_folder, exist_ok=True)
    plt.figure(figsize=(8, 6))
    plt.ylim([1e-30, 10*max(I)])
    plt.plot(v, I, label='T = '+str(T)+' K', linewidth=0.4)
    plt.semilogy()
    #plt.title(database+' '+molecule+' intensity') 
    plt.xlabel('Wavenumber, cm$^{-1}$')
    plt.ylabel('Intensity, cm/molecule')
    plt.legend()
    leg = plt.legend()                  # Get the legend object.
    for line in leg.get_lines():
        line.set_linewidth(1.0)         # Change the line width for the legend.
    plt.savefig(ss_plot_folder+molecule+'__T'+str(T)+'__'+str(min_wn)
                +'-'+str(max_wn)+'__'+database+'.png', dpi=500)
    plt.show()
    print('Cross sections plot saved.')
    
    t.end()
    print('Stick spectra has been saved!\n')  



# Cross section
def extract_broad(broad_df, states_l_df):
    J_df = pd.DataFrame()
    max_broad_J = max(broad_df['Jpp'])
    J_df['Jpp'] = states_l_df['J'].values.astype('float')
    J_df['Jpp'][J_df.Jpp > max_broad_J] = max_broad_J
    if J_df['Jpp'][0] % 1 == 0:
        id_broad = J_df['Jpp']
    else:
        id_broad = J_df['Jpp'] - 0.5
    gamma_L = broad_df['gamma_L'][id_broad].values
    n_air = broad_df['n_air'][id_broad].values
    
    return(gamma_L, n_air)

def linelist_exomol_abs(broad,ratio,cutoff,states_part_df,trans_part_df,default_broad_df,
                        air_broad_df,self_broad_df,H2_broad_df,He_broad_df, CO2_broad_df):
    
    if pd.isna(trans_part_df['v'])[0] == False:
        if cutoff == 'None':
            trans_part_df = trans_part_df[trans_part_df['v'].between(min_wn, max_wn)] 
        else:
            trans_part_df = trans_part_df[trans_part_df['v'].between(min_wn - cutoff, max_wn + cutoff)] 

        id_u = trans_part_df['u'].values
        id_s = states_part_df['id'].values

        trans_part_df.set_index(['u'], inplace=True, drop=False)
        id_us = list(set(id_u).intersection(set(id_s)))
        trans_us_df = trans_part_df.loc[id_us]

        id_l = trans_us_df['l'].values
        id_ls = list(set(id_l).intersection(set(id_s)))
        trans_us_df.set_index(['l'], inplace=True, drop=False)
        trans_s_df = trans_us_df.loc[id_ls]
        trans_s_df.sort_values(by=['v'], inplace=True)

        id_su = trans_s_df['u'].values
        id_sl = trans_s_df['l'].values
        states_u_df = states_part_df.loc[id_su]
        states_l_df = states_part_df.loc[id_sl]

        Epp = states_l_df['E'].values.astype('float')
        gp = states_u_df['g'].values.astype('int')
        A = trans_s_df['A'].values.astype('float')
        v = trans_s_df['v'].values.astype('float')
    else:
        id_u = trans_part_df['u'].values
        id_s = states_part_df['id'].values

        trans_part_df.set_index(['u'], inplace=True, drop=False)
        id_us = list(set(id_u).intersection(set(id_s)))
        trans_us_df = trans_part_df.loc[id_us]

        id_l = trans_us_df['l'].values
        id_ls = list(set(id_l).intersection(set(id_s)))
        trans_us_df.set_index(['l'], inplace=True, drop=False)
        trans_s_df = trans_us_df.loc[id_ls]

        id_su = trans_s_df['u'].values
        id_sl = trans_s_df['l'].values
        states_u_df = states_part_df.loc[id_su]
        states_l_df = states_part_df.loc[id_sl]

        trans_s_df['Ep'] = states_u_df['E'].values.astype('float')
        trans_s_df['Epp'] = states_l_df['E'].values.astype('float')
        trans_s_df['gp'] = states_u_df['g'].values.astype('int')
        trans_s_df['v'] = cal_v(trans_s_df['Ep'].values, trans_s_df['Epp'].values)
        if cutoff == 'None':
            trans_s_df = trans_s_df[trans_s_df['v'].between(min_wn, max_wn)]
        else:
            trans_s_df = trans_s_df[trans_s_df['v'].between(min_wn - cutoff, max_wn + cutoff)]
        trans_s_df.sort_values(by=['v'], inplace=True)
        
        Epp = trans_s_df['Epp'].values
        gp = trans_s_df['gp'].values
        A = trans_s_df['A'].values
        v = trans_s_df['v'].values
        
        id_sl = trans_s_df['l'].values
        states_l_df = states_part_df.loc[id_sl]
     
    gamma_L = pd.DataFrame()
    n_air = pd.DataFrame()
    rows = len(id_sl)
    for i in range(len(broad)):
        if broad[i] == 'Default':
            gamma_L[i] = np.full((1,rows),default_broad_df['gamma_L'][0])[0] * ratio[i]
            n_air[i] = np.full((1,rows),default_broad_df['n_air'][0])[0] * ratio[i]
        if broad[i] == 'Air':
            gamma_L[i] = extract_broad(air_broad_df,states_l_df)[0] * ratio[i]
            n_air[i] = extract_broad(air_broad_df,states_l_df)[1] * ratio[i]
        if broad[i] == 'Self':
            gamma_L[i] = extract_broad(self_broad_df,states_l_df)[0] * ratio[i]
            n_air[i] = extract_broad(self_broad_df,states_l_df)[1] * ratio[i]
        if broad[i] == 'H2':
            gamma_L[i] = extract_broad(H2_broad_df,states_l_df)[0] * ratio[i]
            n_air[i] = extract_broad(H2_broad_df,states_l_df)[1] * ratio[i]
        if broad[i] == 'He':
            gamma_L[i] = extract_broad(He_broad_df,states_l_df)[0]  * ratio[i] 
            n_air[i] = extract_broad(He_broad_df,states_l_df)[1]  * ratio[i] 
        if broad[i] == 'CO2':
            gamma_L[i] = extract_broad(CO2_broad_df,states_l_df)[0]  * ratio[i] 
            n_air[i] = extract_broad(CO2_broad_df,states_l_df)[1]  * ratio[i]
    
    return (A, v, Epp, gp, gamma_L, n_air)

def linelist_exomol_emi(broad,ratio,cutoff,states_part_df,trans_part_df,default_broad_df,
                        air_broad_df,self_broad_df,H2_broad_df,He_broad_df,CO2_broad_df):
    
    if pd.isna(trans_part_df['v'])[0] == False:
        if cutoff == 'None':
            trans_part_df = trans_part_df[trans_part_df['v'].between(min_wn, max_wn)] 
        else:
            trans_part_df = trans_part_df[trans_part_df['v'].between(min_wn - cutoff, max_wn + cutoff)] 

        id_u = trans_part_df['u'].values
        id_s = states_part_df['id'].values

        trans_part_df.set_index(['u'], inplace=True, drop=False)
        id_us = list(set(id_u).intersection(set(id_s)))
        trans_us_df = trans_part_df.loc[id_us]

        id_l = trans_us_df['l'].values
        id_ls = list(set(id_l).intersection(set(id_s)))
        trans_us_df.set_index(['l'], inplace=True, drop=False)
        trans_s_df = trans_us_df.loc[id_ls]
        trans_s_df.sort_values(by=['v'], inplace=True)

        id_su = trans_s_df['u'].values
        id_sl = trans_s_df['l'].values
        states_u_df = states_part_df.loc[id_su]
        states_l_df = states_part_df.loc[id_sl]

        Ep = states_u_df['E'].values.astype('float')
        gp = states_u_df['g'].values.astype('int')
        A = trans_s_df['A'].values.astype('float')
        v = trans_s_df['v'].values.astype('float')
    else:
        id_u = trans_part_df['u'].values
        id_s = states_part_df['id'].values

        trans_part_df.set_index(['u'], inplace=True, drop=False)
        id_us = list(set(id_u).intersection(set(id_s)))
        trans_us_df = trans_part_df.loc[id_us]

        id_l = trans_us_df['l'].values
        id_ls = list(set(id_l).intersection(set(id_s)))
        trans_us_df.set_index(['l'], inplace=True, drop=False)
        trans_s_df = trans_us_df.loc[id_ls]

        id_su = trans_s_df['u'].values
        id_sl = trans_s_df['l'].values
        states_u_df = states_part_df.loc[id_su]
        states_l_df = states_part_df.loc[id_sl]

        trans_s_df['Ep'] = states_u_df['E'].values.astype('float')
        trans_s_df['Epp'] = states_l_df['E'].values.astype('float')
        trans_s_df['gp'] = states_u_df['g'].values.astype('int')
        trans_s_df['v'] = cal_v(trans_s_df['Ep'].values, trans_s_df['Epp'].values)
        if cutoff == 'None':
            trans_s_df = trans_s_df[trans_s_df['v'].between(min_wn, max_wn)]
        else:
            trans_s_df = trans_s_df[trans_s_df['v'].between(min_wn - cutoff, max_wn + cutoff)]
        trans_s_df.sort_values(by=['v'], inplace=True)
        
        Ep = trans_s_df['Ep'].values
        gp = trans_s_df['gp'].values
        A = trans_s_df['A'].values
        v = trans_s_df['v'].values
        
        id_sl = trans_s_df['l'].values
        states_l_df = states_part_df.loc[id_sl]
     
    gamma_L = pd.DataFrame()
    n_air = pd.DataFrame()
    rows = len(id_sl)
    for i in range(len(broad)):
        if broad[i] == 'Default':
            gamma_L[i] = np.full((1,rows),default_broad_df['gamma_L'][0])[0] * ratio[i]
            n_air[i] = np.full((1,rows),default_broad_df['n_air'][0])[0] * ratio[i]
        if broad[i] == 'Air':
            gamma_L[i] = extract_broad(air_broad_df,states_l_df)[0] * ratio[i]
            n_air[i] = extract_broad(air_broad_df,states_l_df)[1] * ratio[i]
        if broad[i] == 'Self':
            gamma_L[i] = extract_broad(self_broad_df,states_l_df)[0] * ratio[i]
            n_air[i] = extract_broad(self_broad_df,states_l_df)[1] * ratio[i]
        if broad[i] == 'H2':
            gamma_L[i] = extract_broad(H2_broad_df,states_l_df)[0] * ratio[i]
            n_air[i] = extract_broad(H2_broad_df,states_l_df)[1] * ratio[i]
        if broad[i] == 'He':
            gamma_L[i] = extract_broad(He_broad_df,states_l_df)[0]  * ratio[i] 
            n_air[i] = extract_broad(He_broad_df,states_l_df)[1]  * ratio[i] 
        if broad[i] == 'CO2':
            gamma_L[i] = extract_broad(CO2_broad_df,states_l_df)[0]  * ratio[i] 
            n_air[i] = extract_broad(CO2_broad_df,states_l_df)[1]  * ratio[i] 
    
    return (A, v, Ep, gp, gamma_L, n_air)

def linelist_hitran_abs(hitran_df):
    '''
    Read HITRAN .par file as the input file.
    Return the data for calculating wavennumbers and cross sections with line profiles.
    
    '''
    A = hitran_df['A'].values
    Epp = hitran_df['Epp'].values
    n_air = hitran_df['n_air'].values
    gamma_air = hitran_df['gamma_air'].values
    gamma_self = hitran_df['gamma_self'].values
    delta_air = hitran_df['delta_air'].values
    gp = hitran_df['gp'].values
    v = hitran_df['v'].values
    #if broad == 'Air':
    #    v = hitran_df['v'].values + delta_air * (P - P_ref) / P
    #else:
    #    v = hitran_df['v'].values

    return (A, v, Epp, gp, n_air, gamma_air, gamma_self, delta_air)


# Line profile
def Doppler_HWHM(v,T):
    '''Return the Doppler half-width at half-maximum (HWHM) -- alpha.'''
    # alpha = np.sqrt(2 * N_A * kB * T * np.log(2) / mass) * v / c
    alpha = ne.evaluate('Sqrt2NAkBln2mInvc * sqrt(T) * v')
    return alpha

def Gaussian_standard_deviation(alpha):
    '''Return the Gaussian standard deviation -- sigma.'''
    # sigma = alpha / np.sqrt(2 * np.log(2))
    sigma = ne.evaluate('alpha * InvSqrt2ln2')
    return sigma

def Lorentzian_HWHM(gamma_L, n_air,T,P):
    '''Return the Lorentzian half-width at half-maximum (HWHM) -- gamma.'''
    gamma = ne.evaluate('gamma_L * (Tref / T)**n_air * (P / Pref)')
    return gamma

def Doppler_profile(dv, alpha):
    '''Return Doppler line profile at dv with HWHM alpha.'''
    # Doppler_profile = np.sqrt(np.log(2) / np.pi) / alpha * np.exp(-np.log(2) * (dv / alpha)**2) 
    DopplerProfile = ne.evaluate('Sqrtln2InvPi / alpha * exp(Negln2 * (dv / alpha)**2)')
    return DopplerProfile

def Gaussian_profile(dv, sigma):
    '''Return Gaussian line profile at dv with HWHM alpha.'''
    # Gaussian_profile = np.exp(- (dv / sigma)**2 / 2)  / (np.sqrt(2 * np.pi) * sigma) 
    GaussianProfile = ne.evaluate('InvSqrt2Pi * exp(-(dv / sigma)**2 * 0.5) / sigma')
    return GaussianProfile

def Lorentzian_profile(dv, gamma):
    '''Return Lorentzian line profile at dv with HWHM gamma.'''
    LorentzianProfile = ne.evaluate('gamma / PI / (dv**2 + gamma**2)')
    return LorentzianProfile

from scipy.special import voigt_profile
def SciPyVoigt_profile(dv, sigma, gamma):
    '''Return the Voigt line profile with Lorentzian component HWHM gamma and Gaussian component HWHM alpha.'''
    SciPyVoigtProfile = voigt_profile(dv, sigma, gamma)
    return SciPyVoigtProfile

from scipy.special import wofz
def SciPyWofzVoigt_profile(dv, sigma, gamma):
    '''Return the Voigt line profile with Lorentzian component HWHM gamma and Gaussian component HWHM alpha.'''
    # scipy_wofz_Voigt_profile = np.real(wofz((dv + 1j*gamma)/sigma/np.sqrt(2))) / sigma / np.sqrt(2*np.pi)
    z = ne.evaluate('(dv + 1j*gamma)/sigma*InvSqrt2')
    wz = wofz(z)
    SciPyWofzVoigtProfile = ne.evaluate('real(wz) / sigma * InvSqrt2Pi')
    return SciPyWofzVoigtProfile

def PseudoVoigt_profile(dv, alpha, gamma, eta, hV):
    '''
    fv is the total FWHM parameter and the Voigt full-width at half-maximum (FWHM), 
    which can be found from the widths of the associated Gaussian and Lorentzian widths.
    eta is a function of Lorentz (fL), Gaussian (fG) and total (fV) full width at half maximum (FWHM) parameters.
    '''
    GaussianProfile = Doppler_profile(dv, hV)
    LorentzianProfile = Lorentzian_profile(dv, hV)
    PseudoVoigtProfile = ne.evaluate('eta * LorentzianProfile + (1 - eta) * GaussianProfile')
    return PseudoVoigtProfile

def PseudoVoigt(alpha, gamma):
    '''
    fv is the total FWHM parameter and the Voigt full-width at half-maximum (FWHM), 
    which can be found from the widths of the associated Gaussian and Lorentzian widths.
    eta is a function of Lorentz (fL), Gaussian (fG) and total (fV) full width at half maximum (FWHM) parameters.
    '''  
    hV = ne.evaluate('(alpha**5+2.69269*alpha**4*gamma+2.42843*alpha**3*gamma**2+4.47163*alpha**2*gamma**3+0.07842*alpha*gamma**4+gamma**5)**0.2')
    eta = ne.evaluate('1.36603*(gamma/hV) - 0.47719*(gamma/hV)**2 + 0.11116*(gamma/hV)**3')
    return (eta, hV)

def PseudoKielkopfVoigt(alpha, gamma):
    '''
    fv is the total FWHM parameter and the Voigt full-width at half-maximum (FWHM), 
    which can be found from the widths of the associated Gaussian and Lorentzian widths.
    eta is a function of Lorentz (fL), Gaussian (fG) and total (fV) full width at half maximum (FWHM) parameters.
    '''
    hV = ne.evaluate('0.5346 * gamma + sqrt(0.2166 * gamma**2 + alpha**2)')
    eta = ne.evaluate('1.36603*(gamma/hV) - 0.47719*(gamma/hV)**2 + 0.11116*(gamma/hV)**3')
    return (eta, hV)

def PseudoOliveroVoigt(alpha, gamma):
    '''
    fv is the total FWHM parameter and the Voigt full-width at half-maximum (FWHM), 
    which can be found from the widths of the associated Gaussian and Lorentzian widths.
    eta is a function of Lorentz (fL), Gaussian (fG) and total (fV) full width at half maximum (FWHM) parameters.
    '''
    d = ne.evaluate('(gamma-alpha)/(gamma+alpha)')
    hV = ne.evaluate('(1-0.18121*(1-d**2)-(0.023665*exp(0.6*d)+0.00418*exp(-1.9*d))*sinPI*d)*(alpha+gamma)')
    eta = ne.evaluate('1.36603*(gamma/hV) - 0.47719*(gamma/hV)**2 + 0.11116*(gamma/hV)**3')
    return (eta, hV)

def PseudoLiuLinVoigt(alpha, gamma):
    '''
    fv is the total FWHM parameter and the Voigt full-width at half-maximum (FWHM), 
    which can be found from the widths of the associated Gaussian and Lorentzian widths.
    eta is a function of Lorentz (fL), Gaussian (fG) and total (fV) full width at half maximum (FWHM) parameters.
    '''
    d = ne.evaluate('(gamma-alpha)/(gamma+alpha)')
    hV = ne.evaluate('(1-0.18121*(1-d**2)-(0.023665*exp(0.6*d)+0.00418*exp(-1.9*d))*sinPI*d)*(alpha+gamma)')
    eta = ne.evaluate('0.68188+0.61293*d-0.18384*d**2-0.11568*d**3')
    return (eta, hV)

def PseudoRoccoVoigt(alpha, gamma):
    '''
    hv is the half width parameter and the Voigt half-width at half-maximum (HWHM), 
    which can be found from the widths of the associated Gaussian and Lorentzian widths.
    eta is a function of Lorentz (fL), Gaussian (fG) and (hV) half width at half maximum (HWHM) parameters.
    '''  
    y = gamma*Sqrtln2/alpha
    erfy = erf(y)
    bhalfy = ne.evaluate('y+Sqrtln2*exp(-0.6055*y+0.0718*y**2-0.0049*y**3+0.000136*y**4)')
    Vy = ne.evaluate('bhalfy*exp(y**2)*(1-erfy)')
    eta = ne.evaluate('(Vy-Sqrtln2)/(Vy*OneminSqrtPIln2)')
    hV = ne.evaluate('alpha / Sqrtln2 * bhalfy')
    return (eta, hV)

def BinnedGaussian_profile(dv, alpha):
    '''Return binned Gaussian line profile at dv with HWHM alpha.'''
    erfxpos = erf(ne.evaluate('Sqrtln2*(dv+binSizeHalf)/alpha'))
    erfxneg = erf(ne.evaluate('Sqrtln2*(dv-binSizeHalf)/alpha'))
    BinnedGaussianProfile = ne.evaluate('(erfxpos-erfxneg)/binSize2')
    return BinnedGaussianProfile

def BinnedLorentzian_profile(dv, gamma):
    '''Return binned Lorentzian line profile at dv with HWHM gamma.'''
    BinnedLorentzianProfile = ne.evaluate('(arctan((binSizeHalf-dv)/gamma)-arctan((-binSizeHalf-dv)/gamma))/binSizePI')
    return BinnedLorentzianProfile

def BinnedVoigt_profile(dv, gamma):
    '''Return binned Voigt line profile at dv with HWHM gamma.'''
    BinnedVoigtProfile = ne.evaluate('(arctan((binSizeHalf-dv)/gamma)-arctan((-binSizeHalf-dv)/gamma))/binSizePI32')
    return BinnedVoigtProfile


# Calculate cross sections
def cross_section_Doppler(wn_grid, v, alpha, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with Doppler profile.
    
    '''  
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):   
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            Doppler = Doppler_profile(dv, alpha)
            _xsec[idx] = ne.evaluate('sum(coef * Doppler)')
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _alpha = alpha[filter]
            _coef = coef[filter]
            Doppler = Doppler_profile(_dv, _alpha)
            _xsec[idx] = ne.evaluate('sum(_coef * Doppler)')
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _alpha = alpha[filter]
            _coef = coef[filter]
            Doppler = Doppler_profile(_dv, _alpha)
            _xsec[idx] = ne.evaluate('sum(_coef * Doppler)')
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _alpha = alpha[filter]
            _coef = coef[filter]
            Doppler = Doppler_profile(_dv, _alpha)
            _xsec[idx] = ne.evaluate('sum(_coef * Doppler)')      
            
    xsec[start:end] += _xsec
    return (xsec)

def cross_section_Gaussian(wn_grid, v, sigma, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with Gaussian profile.
    
    '''  
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):   
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            Gaussian = Gaussian_profile(dv, sigma)
            _xsec[idx] = ne.evaluate('sum(coef * Gaussian)')
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _sigma = sigma[filter]
            _coef = coef[filter]
            Gaussian = Gaussian_profile(_dv, _sigma)
            _xsec[idx] = ne.evaluate('sum(_coef * Gaussian)')
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _sigma = sigma[filter]
            _coef = coef[filter]
            Gaussian = Gaussian_profile(_dv, _sigma)
            _xsec[idx] = ne.evaluate('sum(_coef * Gaussian)')
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _sigma = sigma[filter]
            _coef = coef[filter]
            Gaussian = Gaussian_profile(_dv, _sigma)
            _xsec[idx] = ne.evaluate('sum(_coef * Gaussian)')      
            
    xsec[start:end] += _xsec
    return (xsec)

def cross_section_Lorentzian(wn_grid, v, gamma, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with Lorentzian profile.
    
    '''
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            Lorentzian = Lorentzian_profile(dv, gamma)
            _xsec[idx] = ne.evaluate('sum(coef * Lorentzian)')       
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            Lorentzian = Lorentzian_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * Lorentzian)')       
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            Lorentzian = Lorentzian_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * Lorentzian)')        
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            Lorentzian = Lorentzian_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * Lorentzian)')        
            
    xsec[start:end] += _xsec
    return (xsec)

def cross_section_SciPyVoigt(wn_grid, v, sigma, gamma, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with SciPy Voigt profile.
    
    '''
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            SciPyVoigt = SciPyVoigt_profile(dv, sigma, gamma)
            _xsec[idx] = ne.evaluate('sum(coef * SciPyVoigt)') 
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _sigma = sigma[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            SciPyVoigt = SciPyVoigt_profile(_dv, _sigma, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * SciPyVoigt)') 
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _sigma = sigma[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            SciPyVoigt = SciPyVoigt_profile(_dv, _sigma, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * SciPyVoigt)') 
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _sigma = sigma[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            SciPyVoigt = SciPyVoigt_profile(_dv, _sigma, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * SciPyVoigt)') 
            
    xsec[start:end] += _xsec
    return (xsec)

def cross_section_SciPyWofzVoigt(wn_grid, v, sigma, gamma, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with SciPy Wofz Voigt profile.

    '''
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            SciPyWofzVoigt = SciPyWofzVoigt_profile(dv, sigma, gamma)
            _xsec[idx] = ne.evaluate('sum(coef * SciPyWofzVoigt)')
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _sigma = sigma[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            SciPyWofzVoigt = SciPyWofzVoigt_profile(_dv, _sigma, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * SciPyWofzVoigt)') 
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _sigma = sigma[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            SciPyWofzVoigt = SciPyWofzVoigt_profile(_dv, _sigma, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * SciPyWofzVoigt)')
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _sigma = sigma[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            SciPyWofzVoigt = SciPyWofzVoigt_profile(_dv, _sigma, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * SciPyWofzVoigt)')
            
    xsec[start:end] += _xsec
    return (xsec)    

def cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, hV, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with Pseudo Voigt profile.

    '''
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            PseudoVoigt = PseudoVoigt_profile(dv, alpha, gamma, eta, hV)
            _xsec[idx] = ne.evaluate('sum(coef * PseudoVoigt)')
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _alpha = alpha[filter]
            _gamma = gamma[filter]
            _eta = eta[filter]
            _hV = hV[filter]
            _coef = coef[filter]
            PseudoVoigt = PseudoVoigt_profile(_dv, _alpha, _gamma, _eta, _hV)
            _xsec[idx] = ne.evaluate('sum(_coef * PseudoVoigt)')
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _alpha = alpha[filter]
            _gamma = gamma[filter]
            _eta = eta[filter]
            _hV = hV[filter]
            _coef = coef[filter]
            PseudoVoigt = PseudoVoigt_profile(_dv, _alpha, _gamma, _eta, _hV)
            _xsec[idx] = ne.evaluate('sum(_coef * PseudoVoigt)')
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _alpha = alpha[filter]
            _gamma = gamma[filter]
            _eta = eta[filter]
            _hV = hV[filter]
            _coef = coef[filter]
            PseudoVoigt = PseudoVoigt_profile(_dv, _alpha, _gamma, _eta, _hV)
            _xsec[idx] = ne.evaluate('sum(_coef * PseudoVoigt)')
            
    xsec[start:end] += _xsec
    return (xsec)       

def cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with binned Gaussian profile.
    
    '''  
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):   
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            BinnedGaussian = BinnedGaussian_profile(dv, alpha)
            _xsec[idx] = ne.evaluate('sum(coef * BinnedGaussian)')
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _alpha = alpha[filter]
            _coef = coef[filter]
            BinnedGaussian = BinnedGaussian_profile(_dv, _alpha)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedGaussian)')
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _alpha = alpha[filter]
            _coef = coef[filter]
            BinnedGaussian = BinnedGaussian_profile(_dv, _alpha)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedGaussian)')
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _alpha = alpha[filter]
            _coef = coef[filter]
            BinnedGaussian = BinnedGaussian_profile(_dv, _alpha)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedGaussian)')     
            
    xsec[start:end] += _xsec
    return (xsec)

def cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with binned Lorentzian profile.
    
    '''
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            BinnedLorentzian = BinnedLorentzian_profile(dv, gamma)
            _xsec[idx] = ne.evaluate('sum(coef * BinnedLorentzian)')        
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter] 
            BinnedLorentzian = BinnedLorentzian_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedLorentzian)')    
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            BinnedLorentzian = BinnedLorentzian_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedLorentzian)')        
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            BinnedLorentzian = BinnedLorentzian_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedLorentzian)')          
            
    xsec[start:end] += _xsec
    return (xsec)

def cross_section_BinnedVoigt(wn_grid, v, gamma, coef, cutoff, threshold):
    '''
    Read ExoMol .states, .trans, .pf and .broad files as the input files.
    Return the wavennumbers and cross sections with binned Voigt profile.
    
    '''
    xsec = np.zeros_like(wn_grid)
    if (cutoff == 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            BinnedVoigt = BinnedVoigt_profile(dv, gamma)
            _xsec[idx] = ne.evaluate('sum(coef * BinnedVoigt)')        
    elif (cutoff == 'None') & (threshold != 'None'):
        start = max(0,wn_grid.searchsorted(v.min())-1)
        end = min(wn_grid.searchsorted(v.max()),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = coef >= threshold
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            BinnedVoigt = BinnedVoigt_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedVoigt)')            
    elif (cutoff != 'None') & (threshold == 'None'):
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            BinnedVoigt = BinnedVoigt_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedVoigt)')       
    else: 
        start = max(0,wn_grid.searchsorted(v.min()-cutoff)-1)
        end = min(wn_grid.searchsorted(v.max()+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter_cutoff = np.abs(dv) <= cutoff
            filter_threshold = coef >= threshold
            filter = filter_cutoff & filter_threshold
            _dv = dv[filter]
            _gamma = gamma[filter]
            _coef = coef[filter]
            BinnedVoigt = BinnedVoigt_profile(_dv, _gamma)
            _xsec[idx] = ne.evaluate('sum(_coef * BinnedVoigt)')            
            
    xsec[start:end] += _xsec
    return (xsec)


# Plot and Save Results
def plot_xsec(wn, xsec, database, profile):
    
    plots_foldername = save_path+'/xsecs/plots/'+molecule+'/'+database+'/'
    if os.path.exists(plots_foldername):
        pass
    else:
        os.makedirs(plots_foldername, exist_ok=True)       
    
    xsecs_foldername = save_path+'/xsecs/files/'+molecule+'/'+database+'/'
    if os.path.exists(xsecs_foldername):
        pass
    else:
        os.makedirs(xsecs_foldername, exist_ok=True)

    print('{:25s} : {:<6}'.format('Temperature selected', T), 'K')
    print('{:25s} : {:<6}'.format('Pressure selected', P), 'bar')
    
    '''
    wn_xsec = pd.DataFrame()
    wn_xsec['v'] = wn
    wn_xsec['xsec'] = xsec
    wn_xsec = wn_xsec[wn_xsec.replace([np.inf, -np.inf], np.nan).notnull().all(axis=1)]
    wn = wn_xsec['v'].values
    xsec = wn_xsec['xsec'].values
    '''
    
    #plt.legend(fancybox=True, framealpha=0.0)
    parameters = {'axes.labelsize': 14,
                'legend.fontsize': 14,
                'xtick.labelsize': 12,
                'ytick.labelsize': 12}
    plt.rcParams.update(parameters)
    
    if (wn_wl == 'wn'):
        print('{:25s} : {:<6}'.format('Cutoff is', cutoff), u'cm\u207B\u00B9')
        print('{:25s} : {:<6}'.format('Threshold is', threshold), u'cm\u207B\u00B9/(molecule cm\u207B\u00B2)')
        print('{:25s} : {} {} {} {}'.format('Wavenumber range selected', min_wn, u'cm\u207B\u00B9 -', max_wn, 'cm\u207B\u00B9'))
        
        # Plot cross sections and save it as .png.
        plt.figure(figsize=(12, 6))
        plt.ylim([1e-30, 10*max(xsec)])
        plt.plot(wn, xsec, label='T = '+str(T)+' K, '+profile, linewidth=0.4)
        plt.semilogy()
        #plt.title(database+' '+molecule+' '+abs_emi+' Cross-Section with '+ profile) 
        plt.xlabel('Wavenumber, cm$^{-1}$')
        plt.ylabel('Cross-section, cm$^{2}$/molecule')
        plt.legend()
        leg = plt.legend()                  # Get the legend object.
        for line in leg.get_lines():
            line.set_linewidth(1.0)         # Change the line width for the legend.
        plt.savefig(plots_foldername+molecule+'__T'+str(T)+'__'+wn_wl+str(min_wn)+'-'+str(max_wn)+'__'
                   +database+'__'+abs_emi+'__'+profile+'.png', dpi=500)
        plt.show()
        print('Cross sections plot saved.')

        # Save cross sections into .xsec file.
        xsec_df = pd.DataFrame()
        xsec_df['wavenumber'] = wn
        xsec_df['cross-section'] = xsec
        xsec_filename = (xsecs_foldername+molecule+'__T'+str(T)+'__'+wn_wl+str(min_wn)+'-'+str(max_wn)+'__'
                         +database+'__'+abs_emi+'__'+profile+'.xsec')
        np.savetxt(xsec_filename, xsec_df, fmt="%.06f %20.8e")
        print('Cross sections file saved.')
        
    elif wn_wl == 'wl':
        wl = 10000 / wn
        min_wl = '%.02f' % (10000 / max_wn)
        max_wl = '%.02f' % (10000 / min_wn)
        print('{:25s} : {:<6}'.format('Cutoff is', 10000/cutoff),u'\xb5m')
        print('{:25s} : {:<6}'.format('Threshold is',10000/threshold),u'\xb5m/(moleculeu \xb5m\u00B2)')
        print('{:25s} : {} {} {} {}'.format('Wavelength range selected',min_wl,u'\xb5m -',max_wl,u'\xb5m'))

        # Plot cross sections and save it as .png.
        plt.figure(figsize=(12, 6))
        plt.ylim([1e-30, 10*max(xsec)])
        plt.plot(wl, xsec, label='T = '+str(T)+' K, '+profile, linewidth=0.4)
        plt.semilogy()
        #plt.title(database+' '+molecule+' '+abs_emi+' Cross-Section with '+ profile) 
        plt.xlabel(u'Wavelength, \xb5m')
        plt.ylabel(u'Cross-section, \xb5m\u207B\u00B2/molecule')
        plt.legend()
        leg = plt.legend()                  # Get the legend object.
        for line in leg.get_lines():
            line.set_linewidth(1.0)         # Change the line width for the legend.
        plt.savefig(plots_foldername+molecule+'__T'+str(T)+'__'+wn_wl+str(min_wl)+'-'+str(max_wl)+'__'
                    +database+'__'+abs_emi+'__'+profile+'.png', dpi=500)
        plt.show()
        print('Cross sections plot saved.')
        
        # Save cross sections into .xsec file.
        xsec_df = pd.DataFrame()
        xsec_df['wavelength'] = wl
        xsec_df['cross-section'] = xsec
        xsec_filename = (xsecs_foldername+molecule+'__T'+str(T)+'__'+wn_wl+str(min_wl)+'-'+str(max_wl)+'__'
                         +database+'__'+abs_emi+'__'+profile+'.xsec')
        np.savetxt(xsec_filename, xsec_df, fmt="%.06f %20.8e")
        print('Cross sections file saved.')

    else:
        print('Please type in correct format: wn or wl.')

def get_crosssection(read_path, states_part_df, trans_part_df, hitran_df):
    
    print('Calculate cross-sections.')
    t = Timer()
    t.start()
    
    # ExoMol or HITRAN
    if database == 'ExoMol':
        broad,ratio,default_broad_df,air_broad_df,self_broad_df,H2_broad_df,He_broad_df,CO2_broad_df = read_broad(read_path)
        
        # Read partition function from the ExoMol website
        Q = read_exomolweb_pf(T)
        # Read partition function from local partition function file
        # Q = read_exomol_pf(read_path, T)
        
        # Absorption or emission cross section
        if abs_emi == 'Absorption':
            print('Absorption cross section')
            A, v, Epp, gp, gamma_L, n_air = linelist_exomol_abs(broad,ratio,cutoff,states_part_df,trans_part_df,default_broad_df,
                                                                air_broad_df,self_broad_df,H2_broad_df,He_broad_df,CO2_broad_df)
            coef = cal_abscoefs(v, gp, A, Epp, Q, abundance)
        elif abs_emi == 'Emission': 
            print('Emission cross section')
            A, v, Ep, gp, gamma_L, n_air = linelist_exomol_emi(broad,ratio,cutoff,states_part_df,trans_part_df,default_broad_df,
                                                               air_broad_df,self_broad_df,H2_broad_df,He_broad_df,CO2_broad_df)
            coef = cal_emicoefs(v, gp, A, Ep, Q, abundance)
        else:
            raise ImportError("Please choose one from: 'Absoption' or 'Emission'.")    
        
        # Line profile: Gaussion, Lorentzian or Voigt
        if profile in ['Doppler', 'BinnedDoppler', 'BinnedGaussian']:
            alpha = Doppler_HWHM(v,T)
        elif profile == 'Gaussian':
            alpha = Doppler_HWHM(v,T)
            sigma = Gaussian_standard_deviation(alpha)
        elif profile in ['Lorentzian', 'BinnedLorentzian', 'BinnedVoigt']:
            gamma = pd.DataFrame()
            for i in range(len(broad)):
                gamma[i] = Lorentzian_HWHM(gamma_L[i].values, n_air[i].values,T,P)
            gamma = gamma.sum(axis=1).values
        elif (profile == 'PseudoRoccoVoigt'):
            alpha = Doppler_HWHM(v,T)
            alpha = np.full(len(alpha), 3)
            gamma = pd.DataFrame()
            for i in range(len(broad)):
                gamma[i] = Lorentzian_HWHM (gamma_L[i].values, n_air[i].values,T,P)
            gamma = gamma.sum(axis=1).values  
        elif (profile == 'SciPyVoigt' or profile == 'SciPyWofzVoigt'):
            alpha = Doppler_HWHM(v,T)
            sigma = Gaussian_standard_deviation(alpha)
            gamma = pd.DataFrame()
            for i in range(len(broad)):
                gamma[i] = Lorentzian_HWHM (gamma_L[i].values, n_air[i].values,T,P)
            gamma = gamma.sum(axis=1).values
        else:
            alpha = Doppler_HWHM(v,T)
            gamma = pd.DataFrame()
            for i in range(len(broad)):
                gamma[i] = Lorentzian_HWHM (gamma_L[i].values, n_air[i].values,T,P)
            gamma = gamma.sum(axis=1).values          

    elif database == 'HITRAN':
        print('HITRAN database')
        A, v, Epp, gp, n_air, gamma_air, gamma_self, delta_air = linelist_hitran_abs(hitran_df)
        Q = read_hitran_pf(T)
        # Absorption or emission cross section
        if abs_emi == 'Absorption': 
            print('Absorption cross section') 
            coef = cal_abscoefs(v, gp, A, Epp, Q, abundance)
        elif abs_emi == 'Emission': 
            print('Emission cross section')
            Ep = cal_Ep(Epp, v)
            coef = cal_emicoefs(v, gp, A, Ep, Q, abundance)
        else:
            raise ImportError("Please choose one from: 'Absoption' or 'Emission'.") 
         
        # Line profile: Gaussion, Lorentzian or Voigt
        if profile in ['Doppler', 'BinnedDoppler', 'BinnedGaussian']:
            alpha = Doppler_HWHM(v,T)
        elif profile == 'Gaussian':
            alpha = Doppler_HWHM(v,T)
            sigma = Gaussian_standard_deviation(alpha)
        elif profile in ['Lorentzian', 'BinnedLorentzian', 'BinnedVoigt']:
            gamma_L = gamma_air*ratios[1]+ gamma_self*ratios[2]
            gamma = Lorentzian_HWHM(gamma_L, n_air,T,P)
        elif (profile == 'PseudoRoccoVoigt'):
            alpha = Doppler_HWHM(v,T)
            alpha = np.full(len(alpha), 3)
            gamma = pd.DataFrame()
            for i in range(len(broad)):
                gamma[i] = Lorentzian_HWHM (gamma_L[i].values, n_air[i].values,T,P)
            gamma = gamma.sum(axis=1).values  
        elif (profile == 'SciPyVoigt' or profile == 'SciPyWofzVoigt'):
            alpha = Doppler_HWHM(v,T)
            sigma = Gaussian_standard_deviation(alpha)
            gamma_L = gamma_air*ratios[1]+ gamma_self*ratios[2]
            gamma = Lorentzian_HWHM(gamma_L, n_air,T,P)
        else:
            alpha = Doppler_HWHM(v,T)
            gamma_L = gamma_air*ratios[1]+ gamma_self*ratios[2]
            gamma = Lorentzian_HWHM(gamma_L, n_air,T,P)         

    else:
        raise ImportError("Please add the name of the database 'ExoMol' or 'HITRAN' into the input file.")
    
     
    # Line profiles
    if profile =='Doppler':
        print('Doppler profile')
        xsec = cross_section_Doppler(wn_grid, v, alpha, coef, cutoff, threshold)
    if profile == 'Gaussian':
        print('Gaussion profile')
        xsec = cross_section_Gaussian(wn_grid, v, sigma, coef, cutoff, threshold)
    elif profile == 'Lorentzian':
        print('Lorentzian profile')
        xsec = cross_section_Lorentzian(wn_grid, v, gamma, coef, cutoff, threshold)
    elif profile == 'SciPyVoigt':
        print('SciPy Voigt profile')
        xsec = cross_section_SciPyVoigt(wn_grid, v, sigma, gamma, coef, cutoff, threshold)
    elif profile == 'SciPyWofzVoigt':
        print('SciPy wofz Voigt profile')
        xsec = cross_section_SciPyWofzVoigt(wn_grid, v, sigma, gamma, coef, cutoff, threshold)
    elif profile == 'PseudoVoigt':
        print('Pseudo Voigt profile')
        eta, hV = PseudoVoigt(alpha, gamma)
        xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, hV, coef, cutoff, threshold)       
    elif profile == 'PseudoKielkopfVoigt':
        print('Kielkopf Pseudo Voigt profile')
        eta, hV = PseudoKielkopfVoigt(alpha, gamma)
        xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, hV, coef, cutoff, threshold)       
    elif profile == 'PseudoOliveroVoigt':
        print('Olivero Pseudo Voigt profile')
        eta, hV = PseudoOliveroVoigt(alpha, gamma)
        xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, hV, coef, cutoff, threshold)       
    elif profile == 'PseudoLiuLinVoigt':
        print('Liu-Lin Pseudo Voigt profile')
        eta, hV = PseudoLiuLinVoigt(alpha, gamma)
        xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, hV, coef, cutoff, threshold)       
    elif profile == 'PseudoRoccoVoigt':
        print('Rocco Pseudo Voigt profile')
        eta, hV = PseudoRoccoVoigt(alpha, gamma)
        xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, hV, coef, cutoff, threshold) 
    elif profile == 'BinnedDoppler':
        print('Binned Doppler profile')
        xsec = cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff, threshold)        
    elif profile == 'BinnedGaussian':
        print('Binned Gaussion profile')
        xsec = cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff, threshold)
    elif profile == 'BinnedLorentzian':
        print('Binned Lorentzian profile')
        xsec = cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff, threshold)
    elif profile == 'BinnedVoigt':
        print('Binned Voigt profile')
        xsec = cross_section_BinnedVoigt(wn_grid, v, gamma, coef, cutoff, threshold)           
    else:
        raise ImportError('Please choose line profile from the list.')

    plot_xsec(wn_grid, xsec, database, profile)
    
    t.end()
    
    pass


# Get Results
def get_results(read_path):
    
    t_tot = Timer()
    t_tot.start()
    
    states_part_df = pd.DataFrame()
    trans_part_df = pd.DataFrame() 
    hitran_df = pd.DataFrame()
    
    # ExoMol or HITRAN
    if database == 'ExoMol':
        print('ExoMol database')
        # All functions need whole states.
        states_df = read_all_states(read_path)
        # Only calculating lifetimes and cooling functions need whole transitions.
        NeedAllTrans = Lifetimes + CoolingFunctions
        if NeedAllTrans != 0: 
            all_trans_df = read_all_trans(read_path)
        # Only calculating stick spectra and cross sections need part of the states and transitions.
        NeedPartTrans = StickSpectra + CrossSections
        if NeedPartTrans != 0:
            states_part_df = read_part_states(states_df)
            trans_part_df = read_part_trans(read_path)

            
        # Functions
        Nfunctions = (PartitionFunctions + SpecificHeats + Lifetimes + CoolingFunctions 
                      + StickSpectra + Conversion + CrossSections)
        if Nfunctions > 0:
            if PartitionFunctions == 1:
                exomol_partition_func(states_df, Ntemp, Tmax)
            if SpecificHeats == 1:
                exomol_specificheat(states_df, Ntemp, Tmax)
            if Lifetimes == 1:
                exomol_lifetime(read_path, states_df, all_trans_df)
            if CoolingFunctions == 1:
                exomol_cooling_func(read_path, states_df, all_trans_df, Ntemp, Tmax)
            if StickSpectra == 1:
                exomol_stick_spectra(read_path, states_part_df, trans_part_df, T)
            if CrossSections ==1:
                get_crosssection(read_path, states_part_df, trans_part_df, hitran_df)
        else:   
            raise ImportError("Please choose functions which you want to calculate.")


    elif database == 'HITRAN':
        print('HITRAN database')
        parfile_df = read_parfile(read_path)
        hitran_df = read_hitran_parfile (read_path, parfile_df)
        if CrossSections ==1:
            get_crosssection(read_path, states_part_df, trans_part_df, hitran_df)
 
    else:
        raise ImportError("Please add the name of the database 'ExoMol' or 'HITRAN' into the input file.")
        
    print('\nThe program total running time:')    
    t_tot.end()
    print('\nFinished!')
    
    pass

get_results(read_path)


