# Processing PACE data with Python: get a polarized hyperspectral view of the Earth (EO Summer School 2025)

This README describes the preparatory activities to setup the correct Python environment for the EO Summer School 2025.

## Data access

PACE data can be downloaded from NASA. For this purpose you need a [EARTHDATA login](https://urs.earthdata.nasa.gov/). In the tutorial we will access the data with [earthaccess](https://earthaccess.readthedocs.io/en/latest/).

## Software Requirements

This tutorial was prepared using [Python 3.13.2](https://www.python.org/downloads/release/python-3122/).

The required python packages are installed with [pip](https://pip.pypa.io/en/stable/).

The virtual environment is created with [venv](https://docs.python.org/3/library/venv.html).

Version control is done with [git](https://git-scm.com/).

## Installation

### Create the virtual environment
Create a "eo_summer_school_2025_spexone"-folder on your machine at a suitable location.
Navigate to this folder in your terminal.

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


## Run jupyter notebook
First make sure the virtual environment is activated and you navigated to the project folder. Then run:
```bash
    jupyter lab
```
