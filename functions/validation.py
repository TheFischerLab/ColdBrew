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

def check_env_variables():
    """
    Checks that all required environment variables are set.
    If any environment variable is missing, it prints an error message and exits the program.

    Args:
        None

    Outputs:
        Prints error messages if environment variables are missing.
    
    Returns:
        None
    """
    
    print('checking environment variables...')
    required_vars = ['PYMOL_EXE', 'PHENIX_BIN', 'NACCESS_EXE', 'HBPLUS_EXE', 'EDIASCORER_EXE', 'EDIASCORER_LICENSE']
    for var in required_vars:
        if not os.getenv(var):
            print(f'Error: Required environment variable {var} is not set. See github page for instructions.')
            sys.exit(1)
    print('environment variables found...')

def check_argument_files(args_):
    """
    Checks the validity of the file paths and extensions provided in the arguments.
    Ensures that each file exists, is of the correct type, and the directory is valid.

    Args:
        args_ (Namespace): An object containing the command-line arguments, which include:
            - pdb_file (str): Path to the PDB file.
            - ccp4_file (str): Path to the CCP4 file.
            - mtz_file (str): Path to the MTZ file.
            - outdir (str): Path to the output directory.

    Outputs:
        Raises appropriate errors if any file or directory is invalid or if file extensions are incorrect.
    
    Returns:
        None
    """
    
    files_and_dirs = {
        args_.pdb_file: FileNotFoundError,
        args_.ccp4_file: FileNotFoundError,
        args_.mtz_file: FileNotFoundError,
        args_.outdir: NotADirectoryError
    }

    for path, error in files_and_dirs.items():
        if not os.path.isfile(path) if error == FileNotFoundError else not os.path.isdir(path):
            raise error(path)

    if not args_.pdb_file.endswith('.pdb'):
        raise ValueError('Invalid file extension for ' + args_.pdb_file + '. Expected a ".pdb" file.')
    if not args_.ccp4_file.endswith('.ccp4'):
        raise ValueError('Invalid file extension for ' + args_.ccp4_file + '. Expected a ".ccp4" file.')
    if not args_.mtz_file.endswith('.mtz'):
        raise ValueError('Invalid file extension for ' + args_.mtz_file + '. Expected a ".mtz" file.')
    args_.outdir = os.path.abspath(args_.outdir)

def check_file_exists(raw_datafile_, metric_name_):
    """
    Checks if a file exists at the specified path and raises an error if it is not found.

    Args:
        raw_datafile_ (str): The path to the raw data file that needs to be checked.
        metric_name_ (str): The name of the metric (e.g., RSCC, SASA) associated with the file.

    Outputs:
        Raises a FileNotFoundError if the file does not exist.

    Returns:
        None
    """
    
    if not os.path.isfile(raw_datafile_):
        raise FileNotFoundError('The file ' + raw_datafile_ + ' was not found. There was an error in ' + metric_name_ + ' calculation.')

def check_raw_datafiles(pdb_id__, outdir__, dict_file_suffixes_):
    """
    Checks the existence of raw data files for various metrics (e.g., RSCC, SASA) based on the given suffixes.

    Args:
        pdb_id__ (str): The identifier for the PDB structure.
        outdir__ (str): The directory where the raw data files are located.
        dict_file_suffixes_ (dict): A dictionary mapping metric names to file suffixes.

    Outputs:
        Calls the check_file_exists function for each file. Raises an error if any file is missing.

    Returns:
        None
    """

    for metric_name, suffix in dict_file_suffixes_.items():
        raw_datafile = outdir__ + '/raw_data_files/' + pdb_id__ + suffix
        check_file_exists(raw_datafile, metric_name)
