#!/usr/bin/env python
# coding: utf-8

# In[218]:


ROOT_PATH = None


# # Init

# In[219]:


import contextlib
import os
import platform
import subprocess
from datetime import datetime

import pandas as pd
from pymol import cmd

# # Info Handling

# In[221]:


class DictionaryManager:
    """This class serves the purpose to share information across implemented functions"""

    def __init__(self):
        self.my_dict = {}
        if not ROOT_PATH:
            self.my_dict["root_path"] = os.path.abspath("")
        else:
            self.my_dict["root_path"] = ROOT_PATH

    def update_dict(self, dct):
        """
        Update the existing dictionary with new key-value pairs.
        New keys are added, and existing keys are updated.
        """
        self.my_dict.update(dct)

    def get_dict(self):
        """
        Return the contents of the dictionary.
        """
        return self.my_dict


# In[222]:


info_class = DictionaryManager()


# ## Error Classes

# In[223]:


class SubprocessError(Exception):
    """Raised when subprocess execution fails."""

    def __init__(self, message, output=None):
        self.message = message
        self.output = output
        super().__init__(self.message)


class InvalidSelectionError(Exception):
    """Raised when selection is not appropiate to pymol syntax."""

    def __init__(self, message, output=None):
        self.message = message
        self.output = output
        super().__init__(self.message)


class DirectoryStructureError(Exception):
    """Raised when the directory structure does not match the expected structure."""

    def __init__(self, missing_items):
        self.missing_items = missing_items
        joined_missing = "\n ".join(missing_items)
        self.message = f"The following items are missing: {joined_missing}"
        super().__init__(self.message)


class PipelineError(Exception):
    """Raised when subsequent steps in the program are skipped. I.E hbnet is executed before hbsearch."""

    def __init__(self, message, output=None):
        self.message = message
        self.output = output
        super().__init__(self.message)


class PDBError(Exception):
    """Base class for PDB related errors."""

    pass


class PDBFileNotFoundError(PDBError):
    """Raised when the PDB file does not exist."""

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.message = f"ERROR: {pdb_file} does not exist!"
        super().__init__(self.message)


class InvalidPDBIDError(PDBError):
    """Raised when an invalid PDB ID is provided."""

    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        self.message = f"ERROR: {pdb_id} is not a valid PDB-ID!"
        super().__init__(self.message)


class PDBFetchError(PDBError):
    """Raised when fetching a PDB ID fails."""

    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        self.message = f'ERROR-fetch: unable to fetch "{pdb_id}"'
        super().__init__(self.message)


# ## Hbsearch related

# In[224]:


def get_date():
    """Returns the current date and time"""
    return datetime.now().strftime("%Y-%m-%d-%H:%M:%S")


# In[225]:


def form4cmd(pdbID, line):
    """
    Part of lambda function to format dataframe to Pymol compatible form:
    from "A:183:LEU:O" to "/2akr//A/LEU`183/O"
    :param pdbID: pdb code of handled structure
    :param line: line within dataframe
    """
    chain, resID, amiaci, atom = line.split(":")
    ret = f"/{pdbID}//{chain}/{resID}/{atom}"
    return ret


# In[226]:


def formID(pdbID, x):
    """
    Part of lambda function to form distance ID for pymol object naming:
    from "/2akr//A/LEU`183/O" to "A183O"
    :param pdbstr: pdb code of handled structure
    :param x: entry within dataframe
    """
    pdbstr, chain, resID, atom = filter(None, x.split("/"))

    ret = f"{chain}{resID}{atom}"

    return ret


# In[227]:


def check_directory(script_path):
    """
    Function to validate the directory and file structure for a project based on the provided `script_path`.

    :param script_path: The path to the root directory of the project where the structure will be checked.
    """

    # Expected directory and file structure
    # Dictionary where keys represent different operating systems (Darwin for macOS, Linux, Windows)
    # and the value is a list of expected files or executables for that platform.
    # The empty string "" represents the base directory where general files are expected.
    expected_structure = {
        "Darwin": ["get-water-cluster", "hb-network", "hb-search"],
        "Linux": ["get-water-cluster", "hb-network", "hb-search"],
        "Windows": ["get-water-cluster.exe", "hb-network.exe", "hb-search.exe"],
        "": [
            "hb-define.txt",
            "hb_results",
            "pdb_files",
            "pymol_hbnet.py",
            "pymol_hbnet.ipynb",
            "period-table-info.txt",
            "README.md",
        ],
    }

    # List to store missing files or directories
    missing_items = []

    # Loop through the expected directory structure
    for directory, expected_files in expected_structure.items():
        # Construct the full path to the directory
        dir_path = os.path.join(script_path, directory)

        # If the directory is not the base directory (represented by ""),
        # check whether the directory exists.
        if directory and not os.path.isdir(dir_path):
            # If the directory does not exist, record it as missing
            missing_items.append(f"Missing directory: {dir_path}")
            # Skip checking files inside the directory if the directory itself is missing
            continue

        # Check if all expected files exist in the current directory
        for file in expected_files:
            # Construct the full path to the file
            file_path = os.path.join(dir_path, file)
            # If the file does not exist, record it as missing
            if not os.path.exists(file_path):
                missing_items.append(f"Missing file: {file_path}")

    # If any missing items were found, raise a custom error with the list of missing items
    if missing_items:
        raise DirectoryStructureError(missing_items)
    else:
        # If all directories and files are present, print a success message
        print("Directory structure matches the expected structure.")


# In[228]:


def startHBsearch(pdb_path, hb_file, pse_file, solvent_key, connections):
    """
    Function to initiate an HB search using a specific PDB file, hydrogen bond definition, solvent key,
    and connection settings for the hb-search command. This function dynamically sets the environment
    and selects the appropriate executable based on the operating system.

    :param pdb_path: Path to the PDB file to be processed.
    :param hb_file: Hydrogen bond definition file.
    :param pse_file: PSE file (periodic table information).
    :param solvent_key: Solvent key to define solvent atoms.
    :param connections: Connection file for HB search.
    :return: Output of the hb-search subprocess.
    """

    root_path = info_class.get_dict()["root_path"]

    # Setting environment variable
    os.environ["PSE_FILE"] = pse_file

    # Determine operating system
    osys = platform.system()

    try:
        # Executing hb_search
        hbs = subprocess.run(
            os.path.normpath(
                f"{root_path}/{osys}/hb-search -hb {hb_file} -solv {solvent_key} -con {connections} {pdb_path}"
            ),
            capture_output=True,
            shell=True,
            check=True,
            text=True,
        ).stdout

    except subprocess.CalledProcessError as e:
        # Raise custom error for hb-search failure
        raise SubprocessError(f"HB-search execution failed with error: {e}", e.output)

    return hbs


# In[229]:


def readInHBS(pdb_id, hbs):
    """
    Function to process the output of an HB search, extracting relevant information for hydrogen bonds.
    The function splits the HB search output, filters for rows indicating hydrogen bonds (HBOND),
    and applies transformations to format and identify donor/acceptor atoms.

    :param pdb_id: The PDB ID of the structure being analyzed.
    :param hbs: The HB search output containing hydrogen bond information.
    :return: A pandas DataFrame containing formatted acceptor and donor data along with their IDs.
    """
    try:
        # Split the file into rows and split each row into columns
        hbs_data = [line.split() for line in hbs.strip().split("\n") if line.strip()]

        # Define the column headers
        headers = ["IDENT", "ACC", "sep1", "DONO", ":", "x", "y", "z", "sep2", "a", "b"]

        # Create a DataFrame, filter rows where IDENT is "HBOND", and select relevant columns
        df = pd.DataFrame(hbs_data, columns=headers)
        df = df[df["IDENT"] == "HBOND"][["ACC", "DONO"]]

        # Apply transformations using form4cmd and formID
        df["ACC"] = df["ACC"].apply(lambda x: form4cmd(pdb_id, x))
        df["DONO"] = df["DONO"].apply(lambda x: form4cmd(pdb_id, x))

        df["ACC_ID"] = df["ACC"].apply(lambda x: formID(pdb_id, x))
        df["DONO_ID"] = df["DONO"].apply(lambda x: formID(pdb_id, x))

        return df

    except Exception as e:
        raise ValueError(f"Error processing HBS file: {e}")


# In[245]:


def check_pdb(pdb_str: str, save_loc):
    """
    Function to check if the given PDB string is a file path or a PDB ID and handle accordingly.
    If it's a file path, the function loads the file into PyMOL and returns the PDB ID.
    If it's a 4-character PDB ID, it fetches the structure from an external source.

    :param pdb_str: The PDB file path or PDB ID to be checked.
    :param save_loc: The directory path where the PDB file will be saved if fetched.
    :return: A tuple containing the PDB ID and the file path.
    """

    # Check if input is a file path
    # Check if the input is a file path (by checking if it's a valid file path with a .pdb or .cif extension)
    if os.path.isfile(pdb_str) and pdb_str.lower().endswith((".pdb", ".cif")):
        if not os.path.exists(pdb_str):
            raise PDBFileNotFoundError(pdb_str)
        # Load the PDB file and extract PDB ID
        cmd.load(pdb_str)
        pdb_id = os.path.basename(os.path.normpath(pdb_str)).split(".")[0]
        return pdb_id, pdb_str

    # Check if PDB object is already loaded in PyMOL
    if pdb_str in cmd.get_object_list("all"):
        save_loc = os.path.join(save_loc, pdb_str + ".pdb")
        cmd.save(save_loc)
        return pdb_str, save_loc

    # Handle 4-character numeric PDB IDs
    elif len(pdb_str) == 4 and pdb_str[0].isnumeric():
        # Fetch the PDB ID from an external source
        with contextlib.chdir(save_loc):
            if cmd.fetch(pdb_str) == -1:
                raise PDBFetchError(pdb_str)
        save_loc = os.path.join(save_loc, pdb_str + ".pdb")
        cmd.save(save_loc)
        return pdb_str, save_loc

    # Raise an error for invalid PDB IDs
    else:
        raise InvalidPDBIDError(pdb_str)


# In[231]:


def generate_bond_names(df, pdbID, incr):
    """
    Function to generate unique bond names and distance commands for each hydrogen bond in the dataframe.

    :param df: DataFrame containing hydrogen bond data.
    :param pdbID: PDB ID of the structure being processed.
    :param incr: Incremental number to ensure unique naming.
    :return: DataFrame with bond names and distance commands.
    """
    df["bond_name"] = (
        df["ACC_ID"] + "-" + df["DONO_ID"] + "_" + str(pdbID) + "_" + str(incr)
    )
    df["distance_cmd"] = (
        "distance " + df["bond_name"] + ", " + df["ACC"] + ", " + df["DONO"]
    )
    return df


def execute_distance_commands(df):
    """
    Function to execute distance commands in PyMOL for visualizing hydrogen bonds.

    :param df: DataFrame containing the distance commands to be executed.
    """
    cmd.do(";\n".join(df["distance_cmd"].tolist()))


def create_selections(df):
    """
    Function to create PyMOL selections for donor and acceptor atoms involved in hydrogen bonds.

    :param df: DataFrame containing hydrogen bond data with ACC and DONO columns.
    :return: DataFrame with new columns for donor and acceptor selections.
    """
    df["dono_sel"] = df["DONO"].str.split("/").str[:5].str.join("/")
    df["acc_sel"] = df["ACC"].str.split("/").str[:5].str.join("/")
    return df


def visualize_selections(df):
    """
    Function to visualize donor and acceptor atom selections in PyMOL and then delete them to avoid clutter.

    :param df: DataFrame containing the donor and acceptor selections.
    """
    acc_selection = " + ".join(df["acc_sel"].tolist())
    dono_selection = " + ".join(df["dono_sel"].tolist())

    cmd.select("all_acceptors", acc_selection)
    cmd.delete("all_acceptors")

    cmd.select("all_donors", dono_selection)
    cmd.delete("all_donors")


def group_and_finalize(df, pdbID, incr):
    """
    Function to group hydrogen bonds and finalize their visualization in PyMOL.

    :param df: DataFrame containing bond names and selection data.
    :param pdbID: PDB ID of the structure being processed.
    :param incr: Incremental number to ensure unique naming for groups.
    :return: Tuple containing the list of bond names, the incremented PDB ID, and the final group name.
    """
    cmd.group(f"HBonds_{pdbID}_{incr}", " ".join(df["bond_name"].tolist()))
    cmd.hide("labels", f"HBonds_{pdbID}_{incr}")

    cmd.group(f"HB_{pdbID}_{incr}", pdbID)
    cmd.group(f"HB_{pdbID}_{incr}", f"HBonds_{pdbID}_{incr}")

    cmd.set_name(pdbID, f"{pdbID}_{incr}")

    return df["bond_name"].tolist(), f"{pdbID}_{incr}", f"HB_{pdbID}_{incr}"


def pymolDisplay(df, pdb_id):
    """
    Main function to display hydrogen bonds in PyMOL. It generates bond names, executes distance commands,
    creates selections, visualizes them, and finalizes groupings.

    :param df: DataFrame containing hydrogen bond data.
    :param pdb_id: PDB ID of the structure being visualized.
    :return: Group information for the hydrogen bond visualization in PyMOL.
    """
    group_date = get_date().replace("-", "")[4:]
    incr = 1

    while f"{pdb_id}_{str(incr)}" in cmd.get_names("objects"):
        incr += 1

    # Generate bond names
    df = generate_bond_names(df, pdb_id, incr)

    # Execute distance commands
    execute_distance_commands(df)

    # Create selections for acceptors and donors
    df = create_selections(df)

    # Visualize the selections
    visualize_selections(df)

    # Group and finalize
    return group_and_finalize(df, pdb_id, incr)


# In[234]:


def check_hb_save(root_path: str, hb_save_dir: str, pdb_id: str):
    """
    Checks the hydrogen bond save location, setting it to a default directory if not provided,
    and appends the PDB ID and current date. If the directory does not exist, it creates it.

    :param root_path: Root directory where the hydrogen bond results will be saved.
    :param pdb_str: Directory path to save hydrogen bond results. If not provided, a new directory is created.
    :return: Final path to the directory where hydrogen bond results will be saved.
    """

    if not hb_save_dir:
        # If pdb_str is not provided, set it to a new directory path under the root path
        hb_save_dir = os.path.join(root_path, "hb_results", pdb_id + "_" + get_date())

        # Check if the directory already exists, if not, create it
        if not os.path.exists(hb_save_dir):
            os.makedirs(hb_save_dir)

            return hb_save_dir
    else:
        # If pdb_str is provided, append the pdb_id and date to it
        hb_save_dir = os.path.join(pdb_str, "hb_results", pdb_id + "_" + get_date())

        # Check if the directory already exists, if not, create it
        if not os.path.exists(hb_save_dir):
            os.makedirs(hb_save_dir)
            return hb_save_dir


# In[235]:


def hbsearch(
    pdb_str: str,
    pdb_save_dir: str = None,
    hb_save_dir: str = None,
    hb_file: str = os.path.join(info_class.get_dict()["root_path"], "hb-define.txt"),
    pse_file: str = os.path.join(
        info_class.get_dict()["root_path"], "period-table-info.txt"
    ),
    solvent_key: str = "NONE",
    connections: str = "0",
    remove_pdb: bool = 0,
):
    """
    Executing hb_search with set parameters and extracting HBOND entries from the output.
    The function also checks for the presence of the PDB file and manages saving locations for the results.

    :param pdb_str: The PDB file path or PDB ID to be processed.
    :param pdb_save_loc: Directory path where the PDB file should be saved if fetched.
    :param hb_file: Path to the hydrogen bond definition file.
    :param pse_file: Path to the periodic table info file.
    :param solvent_key: Solvent key for hydrogen bond identification.
    :param connections: Connection settings for HB search.
    :param remove_pdb: If set to 1 the saved pdb will be removed at the end of the function.
    """

    root_path = info_class.get_dict()["root_path"]

    check_directory(root_path)

    if pdb_save_dir is None:
        pdb_save_dir = os.path.join(root_path, "pdb_files/")

    # Check if PDB string is fetchable by pymol, already in pymol session, or available if a custom path was provided
    pdb_id, pdb_path = check_pdb(pdb_str, pdb_save_dir)

    # Check if save location exists, otherwise create
    hb_save_dir = check_hb_save(root_path, hb_save_dir, pdb_id)

    # Run HBsearch
    hbs = startHBsearch(pdb_path, hb_file, pse_file, solvent_key, connections)

    with contextlib.chdir(hb_save_dir):
        with open(f"{pdb_id}.hb", "w") as file:
            # Write the content string to the file
            file.write(hbs)

    # read fromatted output from hb-searcch into dataframe
    df = readInHBS(pdb_id, hbs)

    bondlist, newPDBid, group_name = pymolDisplay(df, pdb_id)

    info_class.update_dict(
        {
            newPDBid: {
                "hb_save_dir": hb_save_dir,
                "pdb_path": pdb_path,
                "pdb_str": hb_save_dir,
                "hb_file_name": f"{pdb_id}.hb",
                "group_name": group_name,
            }
        }
    )

    if remove_pdb:
        os.remove(pdb_path)


# # HB Network

# In[236]:


def exec_hbnet(mol: str):
    """
    Function to execute the hb-network program for a specified molecule.
    This function checks whether the hb-search has been performed for the molecule,
    retrieves the corresponding hydrogen bond file, and runs the hb-network command.

    :param mol: The molecule identifier (PDB ID or structure name) for which the hb-network should be executed.
    :return: The output from the hb-network command.
    """
    info_dct = info_class.get_dict()
    root_path = info_dct["root_path"]

    if not mol in info_dct.keys():
        raise PipelineError(
            f"No hb-file found for {mol}. Make sure hbsearch was executed for the respective molecule"
        )

    # Determine operating system
    osys = platform.system()

    hbfile_path = os.path.join(
        info_dct[mol]["hb_save_dir"], info_dct[mol]["hb_file_name"]
    )
    info_class.update_dict(info_dct)

    try:
        # Executing hb_network
        with contextlib.chdir(info_dct[mol]["pdb_str"]):
            hbn = subprocess.run(
                os.path.normpath(f"{root_path}/{osys}/hb-network {hbfile_path}"),
                capture_output=True,
                shell=True,
                check=True,
                text=True,
            ).stdout

    except subprocess.CalledProcessError as e:
        # Raise custom error for hb-network failure
        raise SubprocessError(f"HB-network execution failed with error: {e}", e.output)

    return hbn


# In[237]:


def eval_cluster(mol: str):
    """
    Function to evaluate and process cluster files for a given molecule.
    It reads the cluster files, deletes empty files, transforms the donor and acceptor information,
    and stores the processed data in a dictionary, which is then updated in the info_class.

    :param mol: The molecule identifier (PDB ID or structure name) for which clusters are being evaluated.
    :return: None (updates the info_class with cluster information).
    """

    info_dct = info_class.get_dict()
    directory = info_dct[mol]["pdb_str"] + "/CLUSTER"
    file_dict = {}
    headers = ["SPEC", "ACC", "SEP", "DONO"]

    # Loop through all files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Check if it is a file
        if os.path.isfile(file_path):
            # Check if the file has size 0 and delete if so
            if os.path.getsize(file_path) == 0:
                os.remove(file_path)
                print(f"Deleted empty file: {filename}")
            else:
                with open(file_path, "r") as file:
                    for_df = [
                        line.split()
                        for line in file.read().strip().split("\n")
                        if line.strip()
                    ]
                    df = pd.DataFrame(for_df, columns=headers)

                    # Apply transformations using form4cmd and formID
                    df["ACC"] = df["ACC"].apply(lambda x: form4cmd(mol, x))
                    df["DONO"] = df["DONO"].apply(lambda x: form4cmd(mol, x))

                    df["ACC_ID"] = df["ACC"].apply(lambda x: formID(mol, x))
                    df["DONO_ID"] = df["DONO"].apply(lambda x: formID(mol, x))
                    file_dict[filename.split(".")[0]] = df

    info_dct[mol]["CLUSTERS"] = file_dict
    info_class.update_dict(info_dct)


# In[238]:


def generate_bond_names_and_commands(df, mol, clus_num):
    """
    Function to generate bond names and distance commands for a cluster.

    :param df: DataFrame containing hydrogen bond information.
    :param mol: The molecule identifier.
    :param clus_num: The cluster number.
    :return: Updated DataFrame with bond_name and distance_cmd columns.
    """
    df["bond_name"] = df["ACC_ID"] + "-" + df["DONO_ID"] + f"_clus-{clus_num}_{mol}"
    df["distance_cmd"] = (
        "distance " + df["bond_name"] + ", " + df["ACC"] + ", " + df["DONO"]
    )
    return df


def generate_selections(df):
    """
    Function to generate donor and acceptor selections for PyMOL visualization.

    :param df: DataFrame containing hydrogen bond information.
    :return: Tuple of acceptor and donor selections.
    """
    df["dono_sel"] = df["DONO"].str.split("/").str[:5].str.join("/")
    df["acc_sel"] = df["ACC"].str.split("/").str[:5].str.join("/")
    acc_selection = " + ".join(df["acc_sel"].tolist())
    dono_selection = " + ".join(df["dono_sel"].tolist())
    return acc_selection, dono_selection


def execute_pymol_commands(
    cluster_name, acc_selection, dono_selection, df, mol, group_name, color, col_bool
):
    """
    Function to execute PyMOL commands for visualizing and grouping clusters.

    :param cluster_name: Name of the cluster.
    :param acc_selection: Acceptor selection command.
    :param dono_selection: Donor selection command.
    :param df: DataFrame containing bond names.
    :param mol: The molecule identifier.
    :param group_name: The name of the group in PyMOL.
    :param color: The color to apply to the cluster.
    """
    pymol_commands = [
        f"select all_acceptors, {acc_selection}",
        "show sticks, all_acceptors",
        f"color {color}, all_acceptors",
        "delete all_acceptors",
        f"select all_donors, {dono_selection}",
        "show sticks, all_donors",
        f"color {color}, all_donors",
        "delete all_donors",
        f'group {cluster_name}, {" ".join(df["bond_name"].tolist())}',
        f"hide labels, {cluster_name}",
        f"color {color}, {cluster_name}",
        f"group Clusters_{mol}, {cluster_name}",
        f"group {group_name}, Clusters_{mol}",
    ]
    if not col_bool:
        pymol_commands = [i for i in pymol_commands if "color" not in i]

    cmd.do(";\n".join(pymol_commands))


def pymol_clusters(mol: str, clus_size: int, col_bool: bool):
    """
    Function to visualize hydrogen bond clusters in PyMOL for a given molecule.

    :param mol: The molecule identifier (PDB ID or structure name).
    """
    info_dct = info_class.get_dict()
    clusters = info_dct[mol]["CLUSTERS"]
    group_name = info_dct[mol]["group_name"]

    # Define a color scheme
    colors = [
        "red",
        "blue",
        "yellow",
        "cyan",
        "orange",
        "purple",
        "teal",
        "brown",
        "lime",
        "grey",
    ]
    color_index = 0

    # Process each cluster
    for key, df in clusters.items():
        cluster_name = f"{key}_{mol}"

        clus_num = key.split("-")[1]

        # Generate bond names and distance commands
        df = generate_bond_names_and_commands(df, mol, clus_num)
        cmp = len(list(set(df.ACC.values.tolist() + df.DONO.values.tolist())))

        # Continue if the unique number of cluster atoms is too small
        if cmp and clus_size and cmp < int(clus_size):
            print(
                f"Number of unique atoms in {cluster_name} smaller than cutoff( {clus_size} ). Skipped! "
            )
            continue

        # Execute distance commands in PyMOL
        execute_distance_commands(df)

        # Generate acceptor and donor selections
        acc_selection, dono_selection = generate_selections(df)

        # Get the color for this cluster, cycle if there are more clusters than colors
        color = colors[color_index % len(colors)]
        color_index += 1

        # Execute batch PyMOL commands for visualization, grouping, and coloring
        execute_pymol_commands(
            cluster_name,
            acc_selection,
            dono_selection,
            df,
            mol,
            group_name,
            color,
            col_bool,
        )


# In[239]:


def hbnet(mol: str, clus_size: int = None, color: bool = 1):
    """
    Function to execute the hydrogen bond network analysis and visualize the results in PyMOL.

    :param mol: The molecule identifier (e.g., PDB ID or structure name).
    :param clus_size: Optional parameter to specify the cluster size (default is None).
    :return: None
    """
    hbn = exec_hbnet(mol)  # Execute the hydrogen bond network analysis for the molecule
    eval_cluster(mol)  # Evaluate and process the clusters for the molecule
    pymol_clusters(mol, clus_size, color)  # Visualize the clusters in PyMOL


# # Show Network

# In[240]:


def format_selection(selection_str):
    """
    Function to format a PyMOL selection string into a truncated form for cross-checking with a DataFrame.

    :param selection_str: PyMOL selection string (e.g., /4akr_1/A/A/SER`3/OG or /4akr_1/A/A/3/OG).
    :return: Tuple containing the PDB object and the formatted string as '/pdb_object//chain/resi_number/atom_name'.
    :raises ValueError: If the selection string format is invalid.
    """

    # Split the selection string by '/'
    parts = selection_str.split("/")

    # Extract the necessary components from the selection string
    try:
        pdb_object = parts[1]  # The PDB object, e.g., 4akr_1
        chain = parts[2]  # Chain identifier, e.g., A

        # Handle cases where the chain may be missing
        if not chain:
            chain = parts[3]

        resi_part = parts[4]  # Residue part, could be SER`3 or just 3
        atom_name = parts[5]  # Atom identifier, e.g., OG

        # Check if the residue part contains a backtick (`), which separates residue name and number
        if "`" in resi_part:
            # Split residue part into residue name and residue number (e.g., SER`3 -> 3)
            resi_number = resi_part.split("`")[1]
        else:
            # If there's no backtick, assume it's already a residue number
            resi_number = resi_part

        # Format the string as '/pdb_object//chain/resi_number/atom_name'
        formatted_str = f"/{pdb_object}//{chain}/{resi_number}/{atom_name}"

        return pdb_object, formatted_str

    except IndexError:
        raise ValueError(f"Invalid selection string format: {selection_str}")


def flatten(lst):
    """
    Function to flatten a list of lists into a single list.

    :param lst: List of lists to be flattened.
    :return: A flattened list containing all elements from the sublists.
    """
    return [item for sublst in lst for item in sublst]


def select_atoms_from_list(atom_list, selection_name="selected_atoms"):
    """
    Function to create a PyMOL selection from a list of atom selection strings.

    :param atom_list: A list of PyMOL atom selections (e.g., ['/4akr_1//A/6/OE1', ...]).
    :param selection_name: The name of the new selection in PyMOL (default is 'selected_atoms').
    :return: The name of the created selection.
    """

    # Join all the selections in the list using 'or' operator for PyMOL
    combined_selection = " or ".join(atom_list)

    # Use cmd.select to create a new selection
    cmd.select(selection_name, combined_selection)
    cmd.zoom(selection_name, 4)
    print(f"Selection '{selection_name}' created with atoms: {atom_list}")

    return selection_name


def shownet(sel_str):
    """
    Function to show the hydrogen bond network related to the selected atom in PyMOL.

    :param sel_str: PyMOL selection string (e.g., '/4akr_1//A/6/OE1').
    :return: Tuple containing the PDB object and the formatted selection string.
    :raises InvalidSelectionError: If the selection string is invalid or no atom is found.
    """

    # Check if there is an object (characters) between the first and second slash
    parts = sel_str.split("/")

    if len(parts) < 3 or not parts[1]:
        raise InvalidSelectionError(
            f"Invalid selection string format: {sel_str}. Please define underlying PDB after first slash."
        )

    model = cmd.get_model(sel_str)

    # Check if the model contains any atoms
    if not model.atom:
        raise InvalidSelectionError(
            f"No atom present according to selection string: {sel_str}"
        )

    elif len(model.atom) > 1:
        raise InvalidSelectionError(
            f"More than one atom selected with: {sel_str}. Please pass string providing single selection"
        )

    # Format the selection string and retrieve the PDB object
    pdb_obj, frm_str = format_selection(sel_str)

    sel_lst = []
    clust_lst = []
    info_dct = info_class.get_dict()

    # Check if the atom is part of any hydrogen bond network
    for clust_num, df in info_dct[pdb_obj]["CLUSTERS"].items():
        acc_lst = df.ACC.values.tolist()
        dono_lst = df.DONO.values.tolist()
        if frm_str in acc_lst or frm_str in dono_lst:
            sel_lst.append(acc_lst + dono_lst)
            clust_lst.append(clust_num)

    # Flatten the list and remove duplicates
    sel_lst = list(set(flatten(sel_lst)))
    # clust_lst = list(set((clust_lst)))

    # Raise an error if the selection is not part of any hydrogen bond network
    if not sel_lst:
        raise InvalidSelectionError(
            f"{sel_str} is not part of a hydrogen bond network."
        )

    # Select atoms from the list in PyMOL
    select_atoms_from_list(sel_lst)
    print(70*'=')
    print('The selection {} is part of following clusters \n {}.'.format(sel_str, clust_lst))
    print(70*'=')
    
    return pdb_obj, frm_str


# ## Extension of pymol commands

# In[194]:


cmd.extend("hbsearch", hbsearch)
cmd.extend("hbnet", hbnet)
cmd.extend("shownet", shownet)


# In[241]:


# # Example usage within jupyter notebook

# hbsearch('4akr')
# hbnet('4akr_1', clus_size=8)
# shownet('/4akr_1/C/C/ASN`261/ND2')
