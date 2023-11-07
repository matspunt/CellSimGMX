.. role:: console(code)
  :language: console
  :class: highlight

Installation and basic usage
============

Introduction
**********

``CellSimGMX`` is a discrete element simulator for GROMACS, aimed especially at epithelial tissues. It can pack particles in predefined shapes on different surfaces, assemble these 'cells' into tissues and stack them on top of a matrix element if needed. In esence, it can be considered a supra coarse-grained method, modelling a single epithelial cell (tens of microns) as a couple hundred particles.

While setup of the necessary input files (coordinates and topologies) is happening in GROMACS formats (.itp/.top/.gro), this would allow conversion to other MD engines (e.g. LAMMPS, CHARMM). Especially as (currently) only simple harmonic bonds and Lennard-Jones potentials are used to parametrize the model, which are widely supported. 

Cloning and installing through :console:`pip` or :console:`conda`
**********

1. Clone the repository in a local directory: :console:`git clone https://github.com/matspunt/CellSimGMX.git`
2. By default, you should be on the :console:`master` branch, with the latest stable code. If you want to make sure you are on this branch, run :console:`git checkout master`. If you want to install or work on the latest code, run :console:`git checkout development`. 
3. Use a virtual environment: :console:`python3 -m venv cellsimgmx` or :console:`conda create -n cellsimgmx`
4. Activate the venv: :console:`source cellsimgmx/bin/activate` or :console:`conda activate cellsimgmx`
5. (pip installation in a pip venv): simply run :console:`pip install .` in the source directory (where :console:`setup.py` is located) 
6. **Recommended method** (pip installation in a conda venv): Make sure you are using the *local* pip version. Run :console:`conda info --envs` and explicitly link to the pip from that environment to prevent problems: :console:`/path/to/miniconda3/envs/cellsimgmx/bin/pip pip install .`. If your console scripts are not loaded (i.e. :console:`cellsimgmx` is not recognized as command) add them to your :console:`~/.bashrc` by running :console:`export PATH="~/.local/bin:$PATH"` and :console:`source ~/.bashrc`
7. (conda installation in a conda venv): Install conda-build if you don't have it: :console:`conda install conda-build`. Build the conda package: :console:`conda build meta.yaml`. The archive will be saved somewhere in your environment files, outputted by the terminal. Install this file: :console:`conda install /home/matspunt/miniconda3/envs/cellsimgmx/conda-bld/linux-64/cellsimgmx-0.0.1-py311_0.tar.bz2 --use-local --offline`. 
8. The package can now be called with :console:`cellsimgmx` or :console:`python3 -m cellsimgmx`. 

.. warning::
   Support for installations through conda is experimental, note that conda works with explicit versioning. To prepare a conda compatible installation, a commit need to be versioned (git tag -a <version>). During development phase, using pip is thus recommended

Running the programme
**********

Running :console:`cellsimgmx -h` will output:

.. code-block:: console

    CellSimGMX: A 3D discrete element framework simulator for epithelial tissues using GROMACS
    
    -h, --help            show this help message and exit
    --ff-dir FF_DIR, -ff FF_DIR
                        Path to the directory where the forcefield is stored (accepts arbitrary .itp file names)
    --input-dir INPUT_DIR, -in INPUT_DIR
                        Path to the directory where the simulation input is stored (accepts arbitrary .json file name)
    --output-dir OUTPUT_DIR, -out OUTPUT_DIR
                        Path where files should be generated and simulation should be run
    --no-sim, -nosim      Does not run a simulation, but does generate all input files like a dryrun.
    --verbose, -v         Optional argument. When enabled, prints detailed logging. Useful for debugging output.

Recommended is to use absolute paths. The basic programme can be run using e.g.:

:console:`cellsimgmx --ff-dir /home/matspunt/CELL/CellSimGMX/cellsimgmx/ff --input-dir /home/matspunt/CELL/CellSimGMX/cellsimgmx/ff --output-dir /home/matspunt/CELL/CellSimGMX/output_dir`

Inputs
**********
The programme uses two basic inputs to construct systems and simulate them. 

1. :console:`input.json`. This is a fixed format configuration file (similar to gromacs .mdp) in which details about the simulation, such as the number of cells, but also Gromacs simulation settings can be set. All options have a description and internally, parameter names are identical to the json inputs. 
2. :console:`forcefield.itp`. This is a GROMACS compatible forcefield.itp file in which the different bead parameters and interactions are defined. These are then parsed, based on input.json settings, to construct the topologies of each cell and the matrix. 