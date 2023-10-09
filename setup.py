from setuptools import setup, find_packages

setup(
    name='CellSimGMX',
    packages='cellsimgmx',
    entry_points = {
        "console_scripts": [
            "cellsimgmx = script_proj.__main__:main",
        ]
    }
    version='0.0.1',    
    description='A 3D discrete element framework simulator for epithelial tissues using',
    url='https://github.com/matspunt/CellSimGMX',
    author='Mats Punt',
    author_email='mats.punt@helsinki.fi',
    license='GPL V3',
    packages=['cellsimgmx'],
    install_requires=['numpy',                     
                      ],
    entry_points = {
    'console_scripts': ['cellsimgx.__main__:main']
}
)

