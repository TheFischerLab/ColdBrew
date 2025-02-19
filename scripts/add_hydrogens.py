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


from pymol import cmd
import argparse

#arg parse functions

def cmd_lineparser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-id', dest='pdb_id', default='NA', type=str,
                        action='store', help='')
    parser.add_argument('-o', dest='outdir', default='.', type=str,
			action='store', help='')

    return parser.parse_args()

def add_hydrogens(pdb_id, outdir):

    cmd.load( outdir + '/' + pdb_id + '_renumber.pdb' )
    cmd.do('h_add')
    cmd.save( outdir + '/' + pdb_id + '_renumber_pymolH.pdb' )

def main():

    args = cmd_lineparser()

    add_hydrogens(args.pdb_id, args.outdir)
    
if __name__ == "pymol":   
    main()
