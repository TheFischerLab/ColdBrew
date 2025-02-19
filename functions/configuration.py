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

def do_setup(pdb_file_, pdb_id_, outdir_): 
    """
    Prepares the input PDB file for further analysis by performing the following steps:
    
    1. Renumbers water molecules in the HETATM records to ensure sequential numbering.
    2. Saves renumbered PDB files, including one without the header.
    3. Generates a separate PDB file containing only water molecules.
    4. Runs an external PyMOL script to add hydrogen atoms to the structure.

    Args:
        pdb_file_ (str): Path to the input PDB file.
        pdb_id_ (str): Identifier for the PDB file, used for output naming.
        outdir_ (str): Directory where output files will be saved.

    Returns:
        None
    """
    
    # renumber waters
    ppdb = PandasPdb().read_pdb(pdb_file_ )
    df_prot = ppdb.df['ATOM']
    df_het = ppdb.df['HETATM']
    chains_cryo = df_prot['chain_id'].unique()
    hets_cryo = df_het['residue_name'].unique()
    df_het_non_HOH = df_het.loc[df_het['residue_name']!='HOH']
    df_het_HOH = df_het.loc[ (df_het['residue_name']=='HOH') & (df_het['element_symbol']=='O') ]
    df_het_HOH = df_het_HOH.copy()
    df_het_HOH['residue_number'] = range(1, len(df_het_HOH) + 1)
    df_het_HOH.loc[:, 'alt_loc'] = ''
    ppdb.df['HETATM'] = pd.concat([df_het_non_HOH, df_het_HOH])
    ppdb.to_pdb(path=outdir_ + '/' + pdb_id_ + '_renumber.pdb', records=['ATOM', 'HETATM', 'OTHERS'] )
    ppdb.to_pdb(path=outdir_ + '/' + pdb_id_ + '_renumber_no_header.pdb', records=['ATOM', 'HETATM'] )

    # save pdb of only water
    df_het_HOH.loc[:, 'atom_number'] = range(1, len(df_het_HOH) + 1)
    ppdb.df['HETATM'] = df_het_HOH
    ppdb.to_pdb(path=outdir_ + '/wats_' + pdb_id_ + '_renumber.pdb', records=['HETATM'] )

    # add hydrogens
    subprocess.run('$PYMOL_EXE -cq ' + os.path.join( os.path.dirname(os.path.abspath(__file__)), '..', 'scripts', 'add_hydrogens.py') + ' -- -id ' + pdb_id_ + ' -o ' + outdir_ + ' > ' + outdir_ + '/pymol_HB_add.log', shell=True, check=True)
