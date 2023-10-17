from setuptools import setup, find_packages

setup(
    name='cellsimgmx',
    version='0.0.1',    
    description='A 3D discrete element framework simulator for epithelial tissues using GROMACS',
    url='https://github.com/matspunt/CellSimGMX',
    author='Mats Punt',
    author_email='mats.punt@helsinki.fi',
    license='GPL V3',
    packages=find_packages(),
    install_requires=['numpy',                     
                      'matplotlib', 
                      'scipy'],
    entry_points = {
        "console_scripts": [
            "cellsimgmx = cellsimgmx.__main__:main",
        ]
    }
)

