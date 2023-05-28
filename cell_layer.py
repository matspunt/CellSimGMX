
##############################################################################
#    This file ("cell_layers.py") contains modules that are useful for the CELL project, but most likely won't be be included in the final Python package. The code may therefore
#    suck even harder than is normal:-)
###############################################################################

import shutil
import os
import subprocess as sp
from distutils.spawn import find_executable

##############################################################################
# A. MDP settings tools
###############################################################################

class Gromacs_basic_tools:
    def write_mdp_file(mdp_dict, filename):
        min_mdp = {
            'integrator': 'cg',
            'nsteps': '5000',
            'emtol': '5',
            'emstep': '0.01',
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
               
        run_mdp = {
            'integrator': 'md',
            'dt': '0.02',
            'nsteps': '350000',
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
            'constraint_algorithm': 'Lincs'
            'gen_vel: yes'
            'gen_seed: -1'
        }
        
        ####separate nsteps out of this and create classes for Langevin, thermostat settings and Barostat settings
        
        with open(filename, "w") as f:
            for key, value in mdp_dict.items():
                f.write(f"{key} = {value}\n")
            
        #move mdp to toppar dir
        if not os.path.exists('toppar'):
            os.makedirs('toppar')
            shutil.move(filename, 'toppar/'+filename)
        else:
            shutil.move(filename, 'toppar/'+filename)
    
  
    def run_GMX_basic(filename='CELL.gro', timeout=40, nr_of_threads=12, top_name='system.top'):
        """
            Basic logic to start and run GROMACS in a fixed folder structure where "toppar"
            contains an .itp and the respective .mdps. 
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




