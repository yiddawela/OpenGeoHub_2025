# Processing PACE data with Python: get a polarized hyperspectral view of the Earth (EO Summer School 2025)

This repository contains the SPEXone data tutorial prepared for the EO Summer School 2025.

## Data access

PACE data can be downloaded from NASA. For this purpose you need a [EARTHDATA login](https://urs.earthdata.nasa.gov/). In the tutorial we will access the data with [earthaccess](https://earthaccess.readthedocs.io/en/latest/).

## Software Requirements

This tutorial was prepared using [Python 3.13.2](https://www.python.org/downloads/release/python-3122/).

The required python packages are installed with [pip](https://pip.pypa.io/en/stable/).

The virtual environment is created with [venv](https://docs.python.org/3/library/venv.html).

Version control is done with [git](https://git-scm.com/).

## Installation

### Clone the repository
To install from source clone the repository. Use either an ssh key (consult [gitlab documentation](https://docs.gitlab.com) for more information) or provide a token:

```bash
    git clone https://<TOKEN_NAME>:<TOKEN_PASSWD>@gitlab.sron.nl/isg_optx/tutorials/eo_summer_school_2025_spexone.git
```

### Create the virtual environment
First navigate to the eo_summer_school_2025_spexone folder in your terminal.

To create the virtual environment use the following command:
```bash
    python -m venv .venv
```

If you have several python executables you can select the relevant one by providing the full path, e.g.
```bash
    /path_to_all_my_python_installations/Python-3.13.2/python -m venv .venv
```

.venv ensures the virtual environment is created in the project folder. Of course the user can provide a different location if desired.

### Dependencies

To install the dependencies first activate the virtual environment.

Using bash (linux/mac os):
```bash
    source .venv/bin/activate
```

Using PowerShell (windows):
```bash
    .venv/Scripts/Activate.ps1
```

Using cmd.exe (windows):
```bash
    .venv/Scripts/Activate.bat
```

Then install the dependencies with pip:
```bash
    python -m pip install -r requirements.txt
```


The repository should have the following structure at the end of the installation process:

```
georegistration
|  README.md                            This file
|  requirements.txt                     File containing the required python packages
|  tutorial.ipynb                       The tutorial.
|  .gitignore                           File letting git know which files should not be added to version control
|--.git                                 Folder
|--.venv                                Folder containing the virtual environment (the python installation)
```

## Run jupyter notebook
First make sure the virtual environment is activated and you navigated to the project folder. Then run:
```bash
    jupyter lab
```
