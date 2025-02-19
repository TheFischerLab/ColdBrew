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
import sys
import subprocess

def edit_RSCC(pdb_id__, outdir__):
    """
    Edits the RSCC raw data file by removing lines corresponding to hydrogens.
                                                                                                                                                                                                                   
    Args:                                                                                                                                                                                                          
        pdb_id__ (str): Identifier for the PDB structure.                                                                                                                                                           
        outdir__ (str): Directory where the raw data files are stored and the edited file will be saved.                                                                                                                                                                                                                   
    Outputs:                                                                                                                                                                                                        
        Saves an edited version of the RSCC data file where certain lines are removed.                                                                                                                              
    Returns:                                                                                                                                                                                                        
        None                                                                                                                                                                                                        
    """

    # read in raw data and remove hydrogen lines
    f = open(outdir__ + '/raw_data_files/' + pdb_id__ + '_original.txt', 'r').readlines()
    fout = open(outdir__ + '/raw_data_files/' + pdb_id__ + '_original_edited.txt','w')
    for line in f:
        if 'HOH' in line and ('H1' in line or 'H2' in line):
            print('Delete:', line.replace('\n', ''))
        else:
            fout.write(line)
    fout.close()

def run_calculations(pdb_file_, pdb_id_, mtz_file_, ccp4_file_, outdir_):
    """                                                                                                                                                                                                            
    Runs various calculations (RSCC, SASA, EDIA, HB) on the given PDB and related files, and saves the raw output data in a specified directory.

    Args:                                                                                                                                                                                                         
        pdb_file_ (str): Path to the PDB file for the structure.                                                                                                                                                   
        pdb_id_ (str): Identifier for the PDB structure.                                                                                                                                                          
        mtz_file_ (str): Path to the MTZ file containing reflection data.                                                                                                                                           
        ccp4_file_ (str): Path to the CCP4 file containing the electron density map.                                                                                                                               
        outdir_ (str): Directory to store the raw data files generated during calculations.                                                                                                                         
                                                                                                                                                                                                                  
    Outputs:                                                                                                                                                                                                        
        Saves various raw data files in the 'raw_data_files' directory.                                                                                                                                           
    Returns:                                                                                                                                                                                                        
        None                                                                                                                                                                                                        
    """
    
    print('running calculations...')
    subprocess.run('mkdir -p ' + outdir_ + '/raw_data_files', shell=True, check=True)

    # RSCC
    print('calculating RSCC...')
    subprocess.run('$PHENIX_BIN/phenix.real_space_correlation ' + pdb_file_ + ' ' + mtz_file_ + ' > ' + outdir_ + '/raw_data_files/' + pdb_id_ + '_original.txt', shell=True, check=True)
    edit_RSCC(pdb_id_, outdir_)

    # SASA
    print('calculating SASA...')
    subprocess.run('$NACCESS_EXE ' + outdir_ + '/' + pdb_id_ + '_renumber_no_header.pdb -w > ' + outdir_ + '/raw_data_files/naccess.log', shell=True, check=True)
    subprocess.run('mv ' + pdb_id_ + '_renumber_no_header.asa ' + outdir_ + '/raw_data_files', shell=True, check=True)
    subprocess.run('mv ' + pdb_id_ + '_renumber_no_header.rsa ' + outdir_ + '/raw_data_files', shell=True, check=True)
    subprocess.run('mv ' + pdb_id_ + '_renumber_no_header.log ' + outdir_ + '/raw_data_files', shell=True, check=True)

    # EDIA
    print('calculating EDIA...')
    subprocess.run('$EDIASCORER_EXE --license $EDIASCORER_LICENSE --target ' + outdir_ + '/' + pdb_id_ + '_renumber.pdb --outputfolder ' + outdir_ + '/raw_data_files --densitymap ' + ccp4_file_ + ' > ' + outdir_ + '/raw_data_files/EDIAscorer.log 2>&1', shell=True, check=True)

    # HB
    print('calculating HB...')
    subprocess.run('$HBPLUS_EXE ' + outdir_ + '/' + pdb_id_ + '_renumber_pymolH.pdb ' + pdb_file_ + ' > ' + outdir_ + '/raw_data_files/hbplus.log', shell=True, check=True)
    subprocess.run('mv ' + pdb_id_ + '_renumber_pymolH.hb2 ' + outdir_ + '/raw_data_files', shell=True, check=True)
    subprocess.run('mv hbdebug.dat ' + outdir_ + '/raw_data_files', shell=True, check=True)
