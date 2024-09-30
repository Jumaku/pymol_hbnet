# pymol_hbnet

This repository contains a collection of Python functions and tools for analyzing hydrogen bond networks using PyMOL. The core functionality is centered around identifying hydrogen bonds in PDB structures and visualizing them in PyMOL. The analysis includes hydrogen bond search (`hbsearch`), hydrogen bond network generation (`hbnet`) and cluster visualisation (`shownet`).

## Prerequisites

- **Python 3.7+** 
- **PyMOL**:
- **Conda**

## Installation

1. Clone the repo
   ```sh
   git clone https://github.com/Jumaku/pymol_hbnet.git
   ```
2. Create conda environment
   ```sh
   conda env create --name {NAME} --file=requirements.yml
   ```
3. Activate environment
   ```shell
   conda activate {NAME}
   ```

## Usage

Since performance issues occured, it is strongly advised to run this script with a pymol session that is launched outside the conda environment. This requires pymol to be installed independently of the conda environment. Run pymol in the root directory of the script after activating the conda environment:

```python
pymol
```

Then run following in the pymol command line: 

```python
import pymol_hbnet # do not put the ending ".py" here
```

In case you want to run the script with pymol outside the root directory you need to move the `pymol_hbnet.py` file to the destination where you need it. 

**Change the following line then in `pymol_hbnet.py`**

```python
ROOT_PATH = None
```

Change this variable to the actual root path after you moved it, so the script knows where to find the subprograms. 

### Usage with jupyter notebook 

In case you do not want to follow the approach above, there is also a jupyter notebook file (`pymol_hbnet.ipynb`) provided in this repository. Run:

```shell
pymol -R
```

This runs pymol in remote mode with the default port of 9123.  Following lines in the notebook connect with the pymol api:

**CAUTION: Make sure pymol uses the port 9123**

```python
import xmlrpc.client as xmlrpclib
cmd = xmlrpclib.ServerProxy('http://localhost:9123')
```


## Available Commands
---

**CAUTION: The commands may take some time. Do not loose patience. Progress can be monitored in console if pymol was launched accordingly.**

**CAUTION: When using custom paths make sure not to use relative paths and not to use the '~' for the home directory.**

---
### hbsearch

**Description**: Initiates a hydrogen bond search for a given PDB structure.

**Usage**:

```python
hbsearch(pdb_str, pdb_save_dir=None, hb_file='path/to/hb-define.txt', pse_file='path/to/period-table-info.txt', solvent_key='NONE', connections='0', remove_pdb=0)

```


- `pdb_str` (str): PDB file path or PDB ID.
- `pdb_save_dir` (str, optional): Directory to save the PDB file. Defaults to `pdb_files/` in the script directory.
- `hb_file` (str, optional): Path to the hydrogen bond definition file. Defaults to `hb-define.txt` in the script directory.
- `pse_file` (str, optional): Path to the periodic table info file. Defaults to `period-table-info.txt` in the script directory.
- `solvent_key` (str, optional): Solvent key for hydrogen bond identification. Default: `None`
- `connections` (str, optional): Connection settings for HB search. Default: `0`
- `remove_pdb` (bool, optional): Deletes pdb_file after run if set to true. Default: `0`

**Example**:

```python
# Within pymol
hbsearch 4akr, hb_save_dir=/path/you/want

# Within jupyter notebook
hbsearch('4akr',hb_save_dir='/path/you/want') 
```
---
### hbnet

**Description**: Analyzes and visualizes the hydrogen bond network of a molecule.

**Usage**:

```python
# Within pymol
hbnet 4akr, clus_size=5

# Within jupyter notebook
hbnet('4akr', clus_size=5)
```



**Parameters**:

- `mol` (str): Molecule identifier
- `clus_size` (int, optional): Minimum number of unique atoms in a cluster to visualize.
---
### shownet

**Description**: Visualizes the hydrogen bond network related to a specific atom selection.

**Usage**:

```python
# Within pymol
shownet sel_str

# Within jupyter notebook
shownet(sel_str)
```

**Parameters**:

- `sel_str` (str): PyMOL selection string (e.g., `/4akr//A/6/OE1`).

---
## Directory Structure

The expected directory structure for the repository is as follows:

```bash
pymol_hbnet
├── hb-define.txt
├── hb_results/
├── pdb_files/
├── period-table-info.txt
├── pymol_hbnet.py
├── pymol_hbnet.ipynb
├── requirements.yml
├── README.md
└── [OS-specific directories]
    ├── Darwin/
    │   ├── get-water-cluster
    │   ├── hb-network
    │   └── hb-search
    ├── Linux/
    │   ├── get-water-cluster
    │   ├── hb-network
    │   └── hb-search
    └── Windows/
        ├── get-water-cluster.exe
        ├── hb-network.exe
        └── hb-search.exe
```
