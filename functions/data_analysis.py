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
import pandas as pd
from biopandas.pdb import PandasPdb
import joblib

def read_in_parsed_data(pdb_id__, outdir__):
    """
    Reads parsed data from renumbered water molecules and extracts relevant metrics.

    This function:
    1. Loads the renumbered water molecule PDB file.
    2. Extracts water residue numbers.
    3. Iterates over a set of predefined metrics, reading corresponding parsed PDB files.
    4. Constructs a DataFrame containing the extracted data.

    Args:
        pdb_id__ (str): Identifier for the PDB file.
        outdir__ (str): Directory where parsed data files are stored.

    Returns:
        pd.DataFrame: A DataFrame containing water residue IDs, chains, and various computed metrics.
    """

    # load renumbered PDB file
    ppdb = PandasPdb().read_pdb( outdir__ + '/wats_' + pdb_id__ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']
    n_wats = len(df_wat.index)
    list_wat_ID = list(df_wat['residue_number'])

    # read in parsed data for each metric
    metrics = ['RSCC_original', 'Bnorm', 'SASA', 'EDIA', 'HB_pymolH_M', 'HB_pymolH_S', 'HB_pymolH_wat']
    dict_metric_to_list = {}
    for metric in metrics:
        parsed_metric_file = outdir__ + '/parsed_data_files/wats_' + pdb_id__ + '_renumber_' + metric + '.pdb'
        if os.path.isfile(parsed_metric_file):
            ppdb = PandasPdb().read_pdb( parsed_metric_file )
            df_wat = ppdb.df['HETATM']
            dict_metric_to_list[metric] = list(df_wat['b_factor'])
        else:
            dict_metric_to_list[metric] = ['M']*n_wats

    # create df with results
    df_out_cur = pd.DataFrame( {'pdb':[pdb_id__]*n_wats, 'wat_ID':list_wat_ID, 'chain':list(df_wat['chain_id']), 'RSCC':dict_metric_to_list['RSCC_original'], 'B_norm':dict_metric_to_list['Bnorm'], 'SASA':dict_metric_to_list['SASA'], 'EDIA':dict_metric_to_list['EDIA'], 'HB_M': dict_metric_to_list['HB_pymolH_M'], 'HB_S':dict_metric_to_list['HB_pymolH_S'] }  )
    df_out_cur['HB'] = df_out_cur['HB_M'] + df_out_cur['HB_S']

    return df_out_cur

def calculate_CB_prob(pdb_file_, pdb_id_, outdir_, df_out_cur):
    """
    Calculates ColdBrew probabilities for water molecules using a pre-trained model.

    This function:
    1. Loads the ColdBrew model from a pre-specified directory.
    2. Applies the model to compute probabilities based on selected features.
    3. Assigns probabilities to water molecules and modifies output data accordingly.
    4. Saves updated probability values into new PDB files.
    5. Saves results for all metrics to csv file.

    Args:
        pdb_file_ (str): Path to the original PDB file.
        pdb_id_ (str): Identifier for the PDB structure.
        outdir_ (str): Directory where processed files are stored.
        df_out_cur (pd.DataFrame): DataFrame containing extracted metrics for water molecules.

    Returns:
        None
    """
    
    print('calculating ColdBrew probabilities...')
    model = joblib.load( os.path.join( (os.path.dirname(os.path.abspath(__file__))) , '..', 'model', 'model.joblib') )

    #apply the model on the df
    feature_cols = ['RSCC', 'B_norm', 'SASA', 'HB', 'EDIA']
    X_ = df_out_cur[feature_cols] # Features
    y_pred_proba = model.predict_proba(X_)[::,1]
    
    df_out_cur['ColdBrew_probability'] = y_pred_proba
            
    #if EDIA = -1, then set CB prob also to -1
    df_out_cur.loc[df_out_cur['EDIA'] == -1, 'ColdBrew_probability'] = -1
            
    #save the results to pdb
    print('saving results')
    #wat pdb
    ppdb = PandasPdb().read_pdb( outdir_ + '/wats_' + pdb_id_ + '_renumber.pdb' )
    df_wat = ppdb.df['HETATM']
    assert len(df_wat.index) == len(df_out_cur.index)
    df_wat = df_wat.copy()
    df_wat['b_factor'] = df_out_cur['ColdBrew_probability'].values
    ppdb.df['HETATM'] = df_wat
    ppdb.to_pdb(path=outdir_ + '/wats_' + pdb_id_ + '_renumber_ColdBrew_probability.pdb')

    #raw pdb
    ppdb = PandasPdb().read_pdb( pdb_file_ )
    df_het = ppdb.df['HETATM']
    df_wat = df_het.loc[ (df_het['residue_name']=='HOH') & (df_het['atom_name']=='O')   ]
    df_het_other = df_het.loc[df_het['residue_name']!='HOH']
    
    assert len(df_wat.index) == len(df_out_cur.index), pdb_ID + ' ' + str(len(df_wat.index)) + ' ' + str(len(df_out_cur.index))
    df_wat = df_wat.copy()
    df_wat['b_factor'] = df_out_cur['ColdBrew_probability'].values
    ppdb.df['HETATM'] = pd.concat([df_wat,df_het_other])
    ppdb.to_pdb(path=outdir_ + '/' + pdb_id_ + '_ColdBrew_probability.pdb', records=['ATOM', 'HETATM'])

    #save results to csv file
    df_out_cur.rename(columns={'wat_ID':'wat_ID_renumbered'}, inplace=True)
    df_out_cur['wat_ID'] = list(df_wat['residue_number'])
    df_out_cur = df_out_cur.reindex(columns=['wat_ID', 'wat_ID_renumbered', 'chain', 'RSCC', 'B_norm', 'SASA', 'EDIA', 'HB', 'ColdBrew_probability'])
    df_out_cur.to_csv(outdir_ + '/' + pdb_id_ + '_ColdBrew_results.csv')
