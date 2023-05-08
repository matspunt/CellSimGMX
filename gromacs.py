import os
import datetime
import subprocess as sp
from distutils.spawn import find_executable
import shutil

class GromacsIO:
    def convert_xyz_to_gro(self, xyz_file, gro_file):
        """Converts an XYZ file to a GRO file (hardcodes resname as 'CELL' but is otherwise non-specific)

        Args:
            xyz_file (str): The path to the input XYZ file.
            gro_file (str): The path to the output GRO file.
        """
        with open(xyz_file, 'r') as file:
            num_atoms = int(file.readline())
            file.readline()

            atoms = []
            for line in file:
                name = line[0:4].strip()
                x, y, z = map(float, line[14:].split())
                atoms.append((name, x, y, z))

        with open(gro_file, 'w') as gro:
            now = datetime.datetime.now()
            header = "GRO file build from {} at {}\n".format(xyz_file, now.strftime("%H:%M:%S"))
            gro.write(header)
            gro.write(f'  {len(atoms)}\n')

            for i, atom in enumerate(atoms):
                name, x, y, z = atom
                x /= 10.0 #.xyz file format works with Angstroms instead of nanometers
                y /= 10.0
                z /= 10.0
                gro.write(f'{i+1:>5d}CELL{name:>5s}{i+1:>5d}{x:>8.3f}{y:>8.3f}{z:>8.3f}\n')
            
            gro.write(f'{0:>10.5f}{0:>10.5f}{0:>10.5f}\n')




####Todo, add other functions to Gromacs_IO class. Compress the .itp code a bit. 



    
def build_itp(gro_file_path, bead_mass, bond_types, lj_pairs):
    """
        Builds itp based on a given coordinate file (in .gro format). This function assumes the nucleus
        bead is at the END of the .gro file (final bead)
        
        Input:
        str "gro_file_path", int "bead mass", dict of "bond_types" and dict of "lj_pairs"
    """
    with open(gro_file_path, "r") as f:
        abspath = os.path.abspath(gro_file_path)
        f.readline()
        n_atoms = int(f.readline())
    
        # Create the .itp file
        itp_file = open("CELL.itp", "w")
        itp_file.write("; Topology file generated from " + abspath + "\n")
    
        # Write the [defaults] section (hardcoded)
        itp_file.write("\n[ defaults ]\n")
        itp_file.write("; nbfunc     comb-rule    gen-pairs  fudgeLJ fudgeQQ\n")
        itp_file.write("    1           2          yes          1.0     1.0\n")
    
        # Write the [atomtypes] section
        itp_file.write("\n[ atomtypes ]\n")
        itp_file.write("; name mass charge   ptype   sigma   epsilon\n")
        itp_file.write("   {:<3s}".format("N") + str(bead_mass) +  "       0.0    A      0.0       0.0\n")
        itp_file.write("   {:<3s}".format("A") + str(bead_mass) +  "       0.0    A      0.0       0.0\n")
        itp_file.write("   {:<3s}".format("B") + str(bead_mass) +  "       0.0    A      0.0       0.0\n")
    
        # Write the [nonbond_params] section
        itp_file.write("\n[ nonbond_params ]\n")
        itp_file.write(";  i   j  func sigma epsilon\n")
        for pair in lj_pairs:
            params = lj_pairs[pair]["params"]
            itp_file.write("   {:<3s} {:<3s} {:<3s} {:<6s} {:<6s}\n".format(pair.split("-")[0], pair.split("-")[1], str(lj_pairs[pair]["func"]), str(params[0]), str(params[1])))
    
        # Write the moleculetype section
        itp_file.write("\n[ moleculetype ]\n")
        itp_file.write("; Name        nrexcl\n")
        itp_file.write("  CELL           1\n\n")
        itp_file.write("[ atoms ]\n")
        itp_file.write("; nr type resnr residue atom cgnr charge mass\n")

        for i in range(n_atoms):
            lines = f.readline()
            #extract the atom names and numbers from the .gro to populate the .itp
            split_lines = lines.split()
            atom_name = str(split_lines[1])
            atom_nr = str(split_lines[2])
            atom_type = atom_name[0]
            itp_file.write("{:<3s}   {:<3s}    1    CELL    {:<3s}   {:<3s} 0.0000   {:<3s} \n".format(atom_nr, atom_type, atom_name, atom_nr, str(bead_mass)))
                
        # Write the [bonds] directive
        itp_file.write("\n[ bonds ]\n")
        itp_file.write("; i j func  r0 fk\n")
    
        f.seek(0) #go to beginning of .gro file to figure out the bonds
        f.readline()
        n_atoms = int(f.readline())
    
        for i in range(n_atoms):
            lines = f.readline()
            split_lines = lines.split()
            atom_nr = str(split_lines[2])
            atom_name = str(split_lines[1])
            atom_type = atom_name[0]
            #Determine the bond type based on the atom names
            bond_type = None
            if atom_type == "A":
                bond_type = "type A-N"
            elif atom_type == "B":
                bond_type = "type B-N"
            # Write the bond parameters if the bond type is defined
            if bond_type is not None:
                params = bond_types[bond_type]["params"]
                itp_file.write("  {:<2s} {:<3s} {:<3s} {:<3s} {:<3s} \n".format(atom_nr, str(n_atoms), str(bond_types[bond_type]["func"]), str(params[0]), str(params[1])))
    itp_file.close()
    

        
def build_top(top_name='system.top'):
    """ 
        Constructs a topology holder based on the .itp
    """
    top_file = open(top_name, "w")
    top_file.write("#include \"toppar/CELL.itp\"\n\n")
    top_file.write("[ system ]\nCell model\n\n[ molecules ]\nCELL       1")
    top_file.close()




#####UNIT TESTING
io = GromacsIO()
io.convert_xyz_to_gro('CELL.xyz', 'CELL.gro')