
##############################################################################
#    This file ("cell_layers.py") contains modules that are useful for the CELL project, but most likely won't be be included in the final Python package.
#    The modules are therefore not well documented, but I think using them is quite straightforward. 
###############################################################################

import os
import shutil
from distutils.spawn import find_executable
import numpy as np
import pandas as pd
from itertools import product #for Cartesian product
import subprocess as sp

##############################################################################
# A. MDP settings tools
###############################################################################

class Gromacs_MDP:
    @staticmethod
    def write_min_mdp_file(filename, em_algorithm, nsteps, emtol, emstep):
            #we don't need to change any other minimization settings than the ones exposed here. 
        min_mdp = {
            'integrator': em_algorithm,
            'nsteps': nsteps,
            'emtol': emtol,
            'emstep': emstep,
            'nstcomm': '100',
            'nstxout': '0',
            'nstvout': '0',
            'nstfout': '0',
            'nstlog': '1000',
            'nstenergy': '100',
            'nstxout-compressed': '1000',
            'compressed-x-precision': '100',
            'compressed-x-grps': 'System',
            'energygrps': 'System',
            'cutoff-scheme': 'Verlet',
            'nstlist': '20',
            'ns_type': 'grid',
            'pbc': 'xyz',
            'verlet-buffer-tolerance': '0.005',
            'coulombtype': 'cutoff',
            'coulomb-modifier': 'Potential-shift-verlet',
            'rcoulomb': '1.1',
            'epsilon_r': '15',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'tcoupl': 'no',
            'pcoupl': 'no',
            'constraints': 'none',
            }
        
        with open(filename, "w") as f:
            for key, value in min_mdp.items():
                f.write(f"{key} = {value}\n")
            
        if not os.path.exists('toppar'):
            os.makedirs('toppar')
            shutil.move(filename, 'toppar/'+filename)
        else:
            shutil.move(filename, 'toppar/'+filename)       
              
    @staticmethod
    def write_run_mdp_file(filename, nsteps, params_dict):

        run_mdp = {
            'integrator': 'md',
            'dt': '0.02',
            'nsteps': nsteps,
            'nstcomm': '100',
            'nstxout': '0',
            'nstvout': '0',
            'nstfout': '0',
            'nstlog': '1000',
            'nstenergy': '100',
            'nstxout-compressed': '1000',
            'compressed-x-precision': '100',
            'compressed-x-grps': 'System',
            'energygrps': 'System',
            'cutoff-scheme': 'Verlet',
            'nstlist': '20',
            'ns_type': 'grid',
            'pbc': 'xyz',
            'verlet-buffer-tolerance': '0.005',
            'coulombtype': 'cutoff',
            'coulomb-modifier': 'Potential-shift-verlet',
            'rcoulomb': '1.1',
            'epsilon_r': '15',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'tcoupl': 'Berendsen',
            'tc-grps': 'System',
            'tau_t': '1.0',
            'ref_t': '298',
            'Pcoupl': 'no',
            'Pcoupltype': 'isotropic',
            'tau_p': '12.0',
            'compressibility': '3e-4',
            'ref_p': '1',
            'constraints': 'none',
            'constraint_algorithm': 'Lincs',
            'gen_vel': 'yes',
            'gen_seed': '-1'
        }
        
        #where params_dict is a dict with the desired change in parameters which overrules the existing .mdp
        run_mdp.update(params_dict)
        
        with open(filename, "w") as f:
            for key, value in run_mdp.items():
                f.write(f"{key} = {value}\n")
        
        #move everything to a folder "toppar" to keep .itp and .mdp nicely contained. 
        if not os.path.exists('toppar'):
            os.makedirs('toppar')
            shutil.move(filename, 'toppar/'+filename)
        else:
            shutil.move(filename, 'toppar/'+filename)
    

##############################################################################
# B. LOGIC TO WORK WITH TOPOLOGIES
###############################################################################

class Gromacs_TOP:
    """
        Class to construct a force field.itp, including LJ definitions with either a single bead type or two bead types
        Since this is a string, we can store the string entries we would like to replace in a dictionary. 
    """
    @staticmethod
    def parse_forcefield_two_bead_types(forcefield_dict):
    
        forcefield= '''[ defaults ]\n1     2    yes    1.0     1.0\n\n[ atomtypes ]\n; name mass   charge   ptype   sigma   epsilon\n    M1   M1_mass       0.0      A      0.0       0.0 \
         \n    M2   M2_mass       0.0      A      0.0       0.0\n    C    C_mass        0.0      A      0.0       0.0\n[ nonbond_params ]\n;  i   j  func sigma epsilon \
            \n    M1    M1   1   M1_M1_sigma   M1_M1_epsilon\n    M2    M2   1   M2_M2_sigma   M2_M2_epsilon\n    M1    M2   1   M1_M2_sigma   M1_M2_epsilon'''

        for key, value in forcefield_dict.items():
            forcefield = forcefield.replace(key, value)

        with open("forcefield.itp", "w") as file:
            file.write(forcefield)

        if not os.path.exists('toppar'):
            os.makedirs('toppar')
            shutil.move('forcefield.itp', 'toppar/forcefield.itp')
        else:
            shutil.move('forcefield.itp', 'toppar/forcefield.itp')

    @staticmethod
    def parse_forcefield_one_bead_type(forcefield_dict):
    
        forcefield= '''[ defaults ]\n1     2    yes    1.0     1.0\n\n[ atomtypes ]\n; name mass   charge   ptype   sigma   epsilon\n   M1   M1_mass       0.0      A      0.0       0.0 \
            \n   C    C_mass        0.0      A      0.0       0.0\n\n[ nonbond_params ]\n;  i   j  func sigma epsilon\n   M1    M1   1   M1_M1_sigma   M1_M1_epsilon\n'''

        for key, value in forcefield_dict.items():
            forcefield = forcefield.replace(key, value)

        with open("forcefield.itp", "w") as file:
            file.write(forcefield)

        if not os.path.exists('toppar'):
            os.makedirs('toppar')
            shutil.move('forcefield.itp', 'toppar/forcefield.itp')
        else:
            shutil.move('forcefield.itp', 'toppar/forcefield.itp')
    
    @staticmethod
    def build_GMX_itp_from_gro(gro_name, bond_dict):
        
        with open(gro_name, "r") as gro:
            #read the .gro and extract the atom information.
            abspath = os.path.abspath(gro_name)
            gro_list = [line.split() for line in gro.readlines()[2:-1]]  # Read and split lines, excluding the first and last lines

            with open("toppar/CELL.itp", "w") as top:
                top.write("[ moleculetype ]\n; Name        nrexcl\n  CELL        1\n\n[ atoms ]\n; nr type resnr residue atom cgnr  charge\n")
                for lines in gro_list:
                    atom_type, atom_nr = str(lines[1]), str(lines[2]) #save the atom name and index from the .GRO
                    atom_name = atom_type[0]
                    top.write("  {:<3s}  {:<3s}  1    CELL    {:<3s}  {:<3s}  0.0000 \n".format(atom_nr, atom_type, atom_name, atom_nr))
                    # we don't need to add the mass since GROMACS pulls it from the force field
                    
                top.write("\n[ bonds ]\n; i j func  r0 fk\n")
                #we extract the atom names (except Center bead) from the bondedtypes dictionary we give the function. 
                atom_names = set()
                for key in bond_dict.keys():
                    atom_names.add(key[1])
                    
                #now we look for the atom in the .gro and parse the right bonded type accordingly
                for lines in gro_list:
                    if lines[1] in atom_names:
                        atom_type, atom_nr = str(lines[1]), str(lines[2]) #save the atom name and index from the .GRO
                        n_atom = str(len(gro_list))  # if the last atom is always the C-bead, we can easily find it like this
                        values = bond_dict[('C', atom_type)] 
                        top.write("  {:<2s} {:<3s} {:<3s} {:<3s} {:<3s} \n".format(atom_nr, n_atom, str(values['func']), str(values['r0']), str(values['fk'])))
            top.close()           

    @staticmethod
    def parse_topol(nr_of_cells='1'):
        """"
            This currently assumes each cell uses the same .itp, which might not be the case in further tests. But then it's easily rewritten. 
        """

        with open("topol.top", "w") as top:
            include_ff =  f"#include \"toppar/forcefield.itp\"\n"
            include_itp = f"#include \"toppar/CELL.itp\"\n"
            system = f"\n[ system ]\nCELL MODEL\n\n[ molecules ]\nCELL    {nr_of_cells}"
            top.write(include_ff)
            top.write(include_itp)
            top.write(system)
            
##############################################################################
# C. PARAMETER SEARCH AND SIMULATION SETUP FOR GROMACS
###############################################################################

class GromacsRun:
    @staticmethod
    def run_GMX_basic(filename='CELL.gro', timeout=2400, nr_of_threads=12, top_name='system.top'):
        """
            Basic logic to start and run GROMACS in a fixed folder structure where "toppar"
            contains an .itp and the respective .mdps. Timeout = 2400 s, is appropriate for 
            about a microsecond simulation to deal with crashes/divisions by zero.
            Assumes single minimization and single run step. 
        """
        GMX = find_executable("gmx")
        if not GMX:
            raise RuntimeError("Cannot find GROMACS executable 'gmx' in PATH")
        try:
            grompp_min_cmd = ['gmx', 'grompp', '-p', top_name, '-f', 'toppar/min.mdp', '-c', filename, '-o', '1-min', '-maxwarn', '1']
            sp.run(grompp_min_cmd, check=True, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)
            mdrun_min_cmd = ['gmx', 'mdrun', '-nt', str(nr_of_threads), '-pin', 'on', '-deffnm', '1-min', '-v']
            sp.run(mdrun_min_cmd, check=True, timeout=timeout, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)

            grompp_run_cmd = ['gmx', 'grompp', '-p', top_name, '-f', 'toppar/run.mdp', '-c', '1-min.gro', '-o', '2-run', '-maxwarn', '1']
            sp.run(grompp_run_cmd, check=True, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)

            mdrun_run_cmd = ['gmx', 'mdrun', '-nt', str(nr_of_threads), '-pin', 'on', '-deffnm', '2-run', '-v']
            sp.run(mdrun_run_cmd, check=True, timeout=timeout, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)

        except sp.TimeoutExpired:
            print('The simulation timed out and was terminated. Log saved as "gmx_run.log"')
            return
        except sp.CalledProcessError as e:
            print(f'The simulation failed with exit code {e.returncode}. Log saved as "gmx_run.log"')
            return
        else:
            print(f'Simulation complete, no problems encountered. Log saved as "gmx_run.log"')

