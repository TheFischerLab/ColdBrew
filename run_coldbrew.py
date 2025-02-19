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
import argparse
import joblib
import pandas as pd
from biopandas.pdb import PandasPdb
import warnings
warnings.filterwarnings("ignore", module="biopandas.pdb")
from functions.validation import *
from functions.configuration import *
from functions.execution import *
from functions.data_parsing import *
from functions.data_analysis import *


def cmd_lineparser():
    """
    Parses command-line arguments for the script.

    Arguments:
    -pdb   : Path to the PDB file (default: 'NA').
    -ccp4  : Path to the CCP4 map file (default: 'NA').
    -mtz   : Path to the MTZ file (default: 'NA').
    -o     : Output directory (default: current directory).

    Returns:
        argparse.Namespace: Parsed arguments with attributes 
        `pdb_file`, `ccp4_file`, `mtz_file`, and `outdir`.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-pdb', dest='pdb_file', type=str, action='store', help='Path to the input PDB file (cryo crystal structure containing water molecules).')
    parser.add_argument('-ccp4', dest='ccp4_file', type=str, action='store', help='Path to the CCP4 map file.')
    parser.add_argument('-mtz', dest='mtz_file', type=str, action='store', help='Path to the MTZ file containing structure factor data.')
    parser.add_argument('-o', dest='outdir',type=str, action='store', help='Path to the output directory where results will be saved. (full path preferred)')

    return parser.parse_args()

def main():
    """
    Main execution function for the pipeline.
    - Parses command-line arguments.
    - Performs initial checks on the environment and input files.
    - Sets up the files needed for calculations.
    - Runs calculations on the input data.
    - Parses and processes raw output files.
    - Computes final results and saves them.
    """

    # parse command-line arguments
    args = cmd_lineparser()

    # validate environment variables, check files, and assign arguments as global variables
    check_env_variables()
    check_argument_files(args)
    for arg_name, value in vars(args).items():
        globals()[arg_name] = value

    # get ID to use for output and setup files for calculation
    pdb_id = pdb_file[pdb_file.rfind('/') + 1: -4]
    print('using ' + pdb_id + ' as the ID...')
    do_setup(pdb_file, pdb_id, outdir)

    # run calculations
    run_calculations(pdb_file, pdb_id, mtz_file, ccp4_file, outdir)

    # parse raw datafiles
    parse_raw_datafiles(pdb_file, pdb_id, outdir)

    # read in parsed data    
    df_out = read_in_parsed_data(pdb_id, outdir)

    # calculate CB prob and save results
    calculate_CB_prob(pdb_file, pdb_id, outdir, df_out)

if __name__ == "__main__":
    main()
