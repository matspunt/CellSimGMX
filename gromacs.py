import os
import datetime
import subprocess as sp
from distutils.spawn import find_executable
import shutil

class Gromacs_IO:
    """
    Methods:
    --------
    convert_xyz_to_gro(xyz_file, gro_file):
        Converts an XYZ file to a GRO file         

    build_GMX_top(gro_path, bond_types, lj_pairs):
        Builds a GROMACS topology file (.itp) based on a given coordinate file (in .gro format). 
    """
    @staticmethod
    def convert_xyz_to_gro(xyz_file, gro_file):
        """Converts an XYZ file to a GRO file (hardcodes resname as 'CELL' but is otherwise non-specific)

        Args:
            xyz_file (str): The path to the input XYZ file.
            gro_file (str): The path to the output GRO file.
        """
        try:
            with open(xyz_file, 'r') as file:
                # skip first two lines
                file.readline()
                file.readline()

                atoms = []
                for line in file:
                    name = line[0:4].strip()
                    x, y, z = map(float, line[13:].split())
                    atoms.append((name, x, y, z))
        except FileNotFoundError:
            raise Exception(f"ERROR: input file '{xyz_file}' not found. Please check that PACKMOL built a suitable coordinate file and try again.  ")

        with open(gro_file, 'w') as gro:
            now = datetime.datetime.now()
            header = "GRO file built from {} at {}\n".format(xyz_file, now.strftime("%H:%M:%S"))
            gro.write(header)
            gro.write(f'  {len(atoms)}\n')

            for i, atom in enumerate(atoms):
                name, x, y, z = atom
                x /= 10.0  # .xyz file format works with Angstroms instead of nanometers
                y /= 10.0
                z /= 10.0
                gro.write(f'{i + 1:>5d}CELL{name:>5s}{i + 1:>5d}{x:>8.3f}{y:>8.3f}{z:>8.3f}\n') #.GRO only saves 3 decimals (single precision)

            gro.write(f'{0:>10.5f}{0:>10.5f}{0:>10.5f}\n') #add an empty box at the end of the .GRO

    def build_GMX_top(gro_path, bead_mass, bond_types, lj_pairs):
        """
            Builds itp based on a given coordinate file (in .gro format). Currently assumes the Nucleus bead is at the end
            (in how the bonds are drawn) but this can be easily made general if needed. 
            
            Input:
            str "gro_path", int "bead mass", dict of "bond_types" and dict of "lj_pairs"
        """
        with open(gro_path, "r") as gro:
            abspath = os.path.abspath(gro_path)
            gro_list = [line.split() for line in gro.readlines()[2:-1]]  # Read and split lines, excluding the first and last lines

            with open("CELL.top", "w") as top:
                #Write the topology header
                now = datetime.datetime.now()
                header = "; Topology file generated from {} at {}\n".format(abspath, now.strftime("%H:%M:%S"))
                top.write(header)
                
                #Write the [defaults] directive. nbfunc = 1 (LJ), comb-rule = 2 (sigma/eps notation). Default LJ settings
                top.write("\n[ defaults ]\n; nbfunc     comb-rule    gen-pairs  fudgeLJ fudgeQQ\n    1           2          yes          1.0     1.0\n")
                bead_names = set(line[1] for line in gro_list if len(line) >= 2)  # Extract only the unique bead identities
                
                #Write the [atomtypes] directive based on the unique beads in the .GRO
                top.write("\n[ atomtypes ]\n; name mass charge   ptype   sigma   epsilon\n")
                top.writelines("{:<4s}{}       0.0    A      0.0       0.0\n".format(name, bead_mass) for name in bead_names)

                # Write the [nonbond_params] directive
                top.write("\n[ nonbond_params ]\n;  i   j  func sigma epsilon\n")
                for pair in lj_pairs:
                    params = lj_pairs[pair]["params"]
                    top.write("   {:<3s} {:<3s} {:<3s} {:<6s} {:<6s}\n".format(pair.split("-")[0], pair.split("-")[1], str(lj_pairs[pair]["func"]), str(params[0]), str(params[1])))

                # Write the [moleculetype] directive and the [atoms] directive based on the .GRO
                top.write("\n[ moleculetype ]\n; Name        nrexcl\n  CELL           1\n\n[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n")
                for lines in gro_list:
                    atom_name, atom_nr = str(lines[1]), str(lines[2])
                    atom_type = atom_name[0]
                    top.write("  {:<3s}   {:<3s}    1    CELL    {:<3s}   {:<3s} 0.0000   {:<3s} \n".format(atom_nr, atom_type, atom_name, atom_nr, str(bead_mass)))
                    
                # Write the [bonds] directive based on the .GRO --> this is currently quite ugly and needs to be rewritten when the force field definitions are finished 
                top.write("\n[ bonds ]\n; i j func  r0 fk\n")
                for line in gro_list:
                    atom_nr = str(line[2])
                    atom_type = str(line[1][0])
                    n_atom = str(len(gro_list))
                    bond_type = None
                    if atom_type == "A":
                        bond_type = "type A-N"
                    elif atom_type == "B":
                        bond_type = "type B-N"
                    if bond_type is not None:
                        params = bond_types[bond_type]["params"]
                        top.write("  {:<2s} {:<3s} {:<3s} {:<3s} {:<3s} \n".format(atom_nr, str(n_atom), str(bond_types[bond_type]["func"]), str(params[0]), str(params[1])))
        
                # Write system information for simulation at the end
                top.write("\n[ system ]\nCELL MODEL\n\n[ molecules ]\nCELL       1")
            top.close()

bead_mass = 72
#predefine the bonded types
bond_types = {"type A-N": {"func": 1, "params": [1.57, 5000]},
              "type B-N": {"func": 1, "params": [1.57, 6000]}}

#LJ is in sigma epsilon form/order.
lj_pairs = {"A-A": {"func": 1, "params": [0.47, 2.0]},
              "B-B": {"func": 1, "params": [0.47, 2.0]},
              "A-B": {"func": 1, "params": [0.47, 2.0]}}

Gromacs_IO().convert_xyz_to_gro('CELL.xyz', 'CELL.gro')
Gromacs_IO.build_GMX_top(gro_path='CELL.gro', bead_mass=72, bond_types=bond_types, lj_pairs=lj_pairs)