import gmxapi as gmx
import os
import logging

# Get the root gmxapi logger.
gmx_logger = logging.getLogger('gmxapi')
# Set a low default logging level
gmx_logger.setLevel(logging.WARNING)
# Make some tools very verbose
#  by descending the hierarchy
gmx_logger.getChild('commandline').setLevel(logging.DEBUG)
#  or by direct reference
logging.getLogger('gmxapi.mdrun').setLevel(logging.DEBUG)

work_dir = '/wrk/matspunt/coding/CELL_MODEL/testing/gmxapi' 

#Creating a .tpr
min_grompp = gmx.commandline_operation('gmx',
                                    arguments=['grompp', '-maxwarn', '-1'],
                                    input_files={'-c': os.path.join(work_dir, 'CELL.gro'),
                                                 '-p': os.path.join(work_dir, 'system.top'),
                                                 '-f': os.path.join(work_dir, 'min.mdp')},
                                    output_files={'-o': os.path.join(work_dir, '1-min.tpr'),
                                                  '-po': os.path.join(work_dir, '1-min.mdp'),
                                                  })

#'commandline_operation' returns 0 if succesful, 1 if failed
if min_grompp.output.returncode.result() != 0:
    #if it failed, show the user shell output to find error
    print(min_grompp.output.stderr.result())

#The .tpr can be direclty run, or .mdp settings can still be changed:
min_tpr = min_grompp.output.file['-o'].result()
min_input = gmx.read_tpr(min_tpr)
modified_tpr = gmx.modify_input(input=min_input, parameters={'nsteps': 10000})
md = gmx.mdrun(input=modified_tpr)
md.run()


