#    ColdBrew
#    Copyright (C) 2024 Justin Seffernick and Marcus Fischer, St. Jude Children's Research Hospital

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import os
import subprocess
import pandas as pd
from biopandas.pdb import PandasPdb
from functions.validation import *

def read_in_RSCC(pdb_id__, outdir__, rscc_raw_datafile_):
    """
    Reads in RSCC values from a raw data file and assigns them to water molecules in a PDB file.

    Args:
        pdb_id__ (str): Identifier for the PDB structure.
        outdir__ (str): Directory where the processed PDB files are stored.
        rscc_raw_datafile_ (str): Path to the raw RSCC data file.

    Outputs:
        Saves a new PDB file with updated RSCC values assigned to water molecules.
    Returns:
        None
    """

    # load renumbered water pdb
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']

    # open the file and get RSCC values for all water molecules
    f = open( rscc_raw_datafile_, 'r' ).readlines()
    RSCC_values = []
    for i, line in enumerate(f):
        if 'HOH' in line[:15]:
            spline = line.split()

            RSCC = float(spline[-3])
            RSCC_values.append(RSCC)

    # make sure the number of RSCC values is correct
    assert len(df_wat['b_factor'])==len(RSCC_values), print('RSCC error', len(df_wat['b_factor']), len(RSCC_values))

    # save results
    df_wat['b_factor'] = RSCC_values
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_RSCC_original.pdb', records=['HETATM'])


#def calc_avg_stdev(df_, col_code):
#    return [df_[col_code].mean(), df_[col_code].std()]

def calc_avg_stdev_bfactor(df_):
    """
    Calculates the average and standard deviation of the B-factor values from a DataFrame.

    Args:
        df_ (pandas.DataFrame): DataFrame containing the 'b_factor' column.

    Outputs:
        list: A list containing the mean and standard deviation of the 'b_factor' values.
    Returns:
        list: [mean, standard deviation]
    """
    return [df_['b_factor'].mean(), df_['b_factor'].std()]

def read_in_B_norm(pdb_file__, pdb_id__, outdir__):
    """
    Reads in a PDB file, normalizes the B-factor values of water molecules,
    and saves the updated PDB file.

    Args:
        pdb_file__ (str): Path to the protein PDB file.
        pdb_id__ (str): Identifier for the PDB structure (used for file naming).
        outdir__ (str): Directory where the processed PDB files are stored.

    Outputs:
        Saves a new PDB file with normalized B-factor values (Bnorm) for water molecules.
    Returns:
        None
    """
    
    # load renumbered water pdb
    ppdb = PandasPdb().read_pdb( pdb_file__ )
    df_prot = ppdb.df['ATOM']

    # calculate average and standard deviation of protein b factor
    B_avg_stdev_protein = calc_avg_stdev_bfactor(df_prot)

    #calculate normalized b factor and save
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']
    df_wat['b_factor'] = ( df_wat['b_factor'] - B_avg_stdev_protein[0] ) / B_avg_stdev_protein[1]
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_Bnorm.pdb', records=['HETATM'])

def read_in_SASA(pdb_id__, outdir__, SASA_raw_datafile_pdb_):
    """
    Reads in a PDB file with SASA (solvent-accessible surface area) data,
    copies the occupancy values to the B-factor column, and saves the updated PDB file.

    Args:
        pdb_id__ (str): Identifier for the PDB structure (used for file naming).
        outdir__ (str): Directory where the processed PDB files are stored.
        SASA_raw_datafile_pdb_ (str): Path to the PDB file containing SASA data.

    Outputs:
        Saves a new PDB file with SASA values copied to the B-factor column for water molecules.
    Returns:
        None
    """

    # load pdb from naccess, reformat, and save
    ppdb = PandasPdb().read_pdb( SASA_raw_datafile_pdb_ )
    df_wat = ppdb.df['HETATM']
    df_wat['b_factor'] = df_wat['occupancy'] #naccess writes the value as occupancy, but for consistency, copy that to B factor
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_SASA.pdb', records=['HETATM'])

def read_in_EDIA(pdb_id__, outdir__, edia_raw_datafile_):
    """
    Reads in EDIA values from a CSV file, assigns them to water molecules in a PDB file,
    and saves the updated PDB file.

    Args:
        pdb_id__ (str): Identifier for the PDB structure (used for file naming).
        outdir__ (str): Directory where the processed PDB files are stored.
        edia_raw_datafile_ (str): Path to the raw EDIA data file in CSV format.

    Outputs:
        Saves a new PDB file with EDIA values assigned to the B-factor column of water molecules.
    Returns:
        None
    """

    # load renumbered water pdb
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']

    # read in EDIA data
    df_temp = pd.read_csv(edia_raw_datafile_)
    df_edia_wat = df_temp.loc[ (df_temp['Structure specifier']=='w') & (df_temp['Substructure name']=='HOH') ]

    # make sure the number of EDIA values is correct
    assert len(df_wat['b_factor']) == len(df_edia_wat['EDIA'])

    # save results
    df_wat['b_factor'] = list(df_edia_wat['EDIA'])
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_EDIA.pdb', records=['HETATM'])

def read_in_HB(pdb_id__, outdir__, HB_raw_datafile_):
    """
    Reads in hydrogen bond (HB) information from a raw data file, assigns the number of hydrogen bonds 
    to water molecules in a PDB file based on their interactions with proteins (main chain, side chain) 
    and other water molecules, and saves the updated PDB files.

    Args:
        pdb_id__ (str): Identifier for the PDB structure (used for file naming).
        outdir__ (str): Directory where the processed PDB files are stored.
        HB_raw_datafile_ (str): Path to the raw hydrogen bond data file.

    Outputs:
        Saves new PDB files with updated B-factor values corresponding to the number of hydrogen bonds 
        for water molecules interacting with protein.
    Returns:
        None
    """
    
    # load renumbered water pdb
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']
    n_wats = len(df_wat.index)

    # read in parse HB data 
    dict_wat_ID_to_n_HB_M = {} #HB to main chain protein
    dict_wat_ID_to_n_HB_S = {} #HB to side chain protein
    for i in range(1, n_wats+1):
        dict_wat_ID_to_n_HB_M[i] = 0
        dict_wat_ID_to_n_HB_S[i] = 0

    water_partner_pairs = []

    f = open(HB_raw_datafile_, 'r')
    for line in f:
        if 'HOH' in line:
            spline = line.split()    
            if len(spline)==13:
                HB_type = spline[5].upper()
                if 'HOH' in spline[0]: 
                    wat_num = int(spline[0][1:5])
                    ID_wat = spline[0] 
                    ID_part = spline[2] 
                    if (ID_wat, ID_part) not in water_partner_pairs:
                        if HB_type in ['HM','MH']:
                            dict_wat_ID_to_n_HB_M[wat_num] = dict_wat_ID_to_n_HB_M[wat_num] + 1
                        elif HB_type in ['HS','SH']:
                            dict_wat_ID_to_n_HB_S[wat_num] = dict_wat_ID_to_n_HB_S[wat_num] + 1
                        water_partner_pairs.append((ID_wat, ID_part))
                if 'HOH' in spline[2]: 
                    wat_num = int(spline[2][1:5])
                    ID_wat = spline[2] 
                    ID_part = spline[0] 
                    if (ID_wat, ID_part) not in water_partner_pairs:
                        if HB_type in ['HM','MH']:
                            dict_wat_ID_to_n_HB_M[wat_num] = dict_wat_ID_to_n_HB_M[wat_num] + 1
                        elif HB_type in ['HS','SH']:
                            dict_wat_ID_to_n_HB_S[wat_num] = dict_wat_ID_to_n_HB_S[wat_num] + 1
                        water_partner_pairs.append((ID_wat, ID_part))
    # save results
    # main chain (M)
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']
    df_wat['b_factor'] = df_wat['residue_number'].map(dict_wat_ID_to_n_HB_M)
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_HB_pymolH_M.pdb', records=['HETATM'])

    # side chain (S)
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']
    df_wat['b_factor'] = df_wat['residue_number'].map(dict_wat_ID_to_n_HB_S)
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_HB_pymolH_S.pdb', records=['HETATM'])


def parse_raw_datafiles(pdb_file_, pdb_id_, outdir_):
    """
    Parses raw data files for a given PDB ID, validates their existence, and processes them into parsed data files.

    Args:
        pdb_file_ (str): Path to the PDB file.
        pdb_id_ (str): Identifier for the PDB structure (used for file naming).
        outdir_ (str): Directory where the parsed data files will be saved.

    Outputs:
        Saves parsed data files in the specified output directory.
    Returns:
        None
    """
    
    dict_file_suffixes = {
	'EDIA': '_renumberatomscores.csv',
	'RSCC': '_original_edited.txt',
	'HB': '_renumber_pymolH.hb2',
	'SASA': '_renumber_no_header.asa'
    }

    # check raw datafiles and copy SASA result to pdb file
    check_raw_datafiles(pdb_id_, outdir_, dict_file_suffixes)
    SASA_raw_datafile = outdir_ + '/raw_data_files/' + pdb_id_ + dict_file_suffixes['SASA']
    SASA_raw_datafile_pdb = outdir_ + '/raw_data_files/' + pdb_id_ + '_renumber_asa.pdb'
    subprocess.run('cp ' + SASA_raw_datafile + ' ' + SASA_raw_datafile_pdb, shell=True, check=True)

    print('parsing raw datafiles...')

    # parse the datafiles
    subprocess.run('mkdir -p ' + outdir_ + '/parsed_data_files', shell=True, check=True)
    read_in_RSCC(pdb_id_, outdir_, outdir_ + '/raw_data_files/' + pdb_id_ + dict_file_suffixes['RSCC'])
    read_in_B_norm(pdb_file_, pdb_id_, outdir_)
    read_in_SASA(pdb_id_, outdir_, SASA_raw_datafile_pdb)
    read_in_EDIA(pdb_id_, outdir_, outdir_ + '/raw_data_files/' + pdb_id_ + dict_file_suffixes['EDIA'])
    read_in_HB(pdb_id_, outdir_, outdir_ + '/raw_data_files/' + pdb_id_ + dict_file_suffixes['HB'])
