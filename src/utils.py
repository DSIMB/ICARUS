"""
 This module is composed of few useful functions, used by multiple classes.
"""

import os
import platform
import random
import shutil
import sys
from collections import OrderedDict


def get_random_name():
    """
    Used to generate random name for intermediate protein structures
    """
    pool = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    return "".join(random.choices(pool, k=10))

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK_DIR = os.path.join(os.getcwd(), "icarus_output")
plt = platform.system()
# Use the native tmpfs partition on linux systems to boost i/o perfs 
if plt == "Linux":
    TMP_DIR = os.path.join("/dev/shm", get_random_name())
# Use default /tmp directory because not tmpfs partition on Mac or other OS
elif plt == "Darwin" or not os.path.exists("/dev/shm"):
    TMP_DIR = os.path.join("/tmp", get_random_name())

def clean():
    """
    Clean working dir in /dev/shm/<random> or /tmp/<random> folder.
    """
    try:
        if os.path.exists(TMP_DIR):
            shutil.rmtree(TMP_DIR)
    except OSError as e:
        print("Failed to delete tmp dir '%s'. Reason: %s" % (TMP_DIR, e))


def product(liste):
    """
    Returns the product of the items of a given list
    """
    prod = liste[0]
    for i in liste[1:]:
        prod *= i
    return prod


def get_confirm():
    """
    Not used currently, asks user for confirmation to compute a long seg_level.

    Returns
    -------
        - bool
    """
    c = "a"
    while c not in ["y", "n", ""]:
        c = str(input())
    if c in ["y", ""]:
        return True
    return False


def reformat_atom(atom, new_n, new_res, new_chain=None):
    """
    Given an ATOM record from a pdb file as a string, change the atom number
    and res number

    ALL CREDIT GOES TO PIERRE POULAIN https://cupnet.net/pdb-format/

    Parameters
    ----------
        - atom: str
            ATOM record from a pdb file
        - new_n: int
            the new atom number that will be set for this record
        - new_res: int
            the new res number that will be set for this record
    Returns
    -------
        - : str
            updated ATOM record as a string
    """
    a = [x for x in range(15)]
    a[0] = atom[0:6]  # "ATOM  " or "HETATM"
    a[1] = int(atom[6:11])  # atom serial number
    a[2] = atom[12:16]  # atom name
    a[3] = atom[16:17]  # alternate location indicator
    a[4] = atom[17:20]  # residue name
    a[5] = atom[21:22]  # chain identifier
    a[6] = int(atom[22:26])  # residue sequence number
    a[7] = atom[26:27]  # code for insertion of residues
    a[8] = float(atom[30:38])  # orthogonal coordinates for X (in Angstroms)
    a[9] = float(atom[38:46])  # orthogonal coordinates for Y (in Angstroms)
    a[10] = float(atom[46:54])  # orthogonal coordinates for Z (in Angstroms)
    a[11] = str(atom[54:60])  # occupancy
    a[12] = str(atom[60:66])  # temperature factor
    a[13] = atom[76:78]  # element symbol
    a[14] = atom[78:80]  # charge on the atom
    a[1] = new_n
    a[6] = new_res
    if new_chain is not None:
        a[5] = new_chain
    return "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   "\
           "{:8.3f}{:8.3f}{:8.3f}{:6.2s}{:6.2s}         "\
           "{:>2s}{:2s}\n".format(a[0], a[1], a[2], a[3], a[4], a[5], a[6],
                                  a[7], a[8], a[9], a[10], a[11], a[12], a[13],
                                  a[14])


def reformat_struct(path, ori_res_num_and_chain, chain, new_path=""):
    """
    Reset atoms and residues number of a newly merged pdb file hence ordering
    the records properly according to the alignments.
    Also remove unsolved residues (ones that have no CA atom).
    Changes are applied directly in the file indicated by path, or, if defined, to new_path.

    Args:
        - path (str): path to the structural file
        - new_path (str): if defined, path to output clean pdb file
        - ori_res_num_and_chain (dict): mapping of old to new residues numerotation
        - chain (str): chain of the PDB to use

    Returns:
        None
    """
    three_to_one = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M"
    }
    # initialize new atom number at 0
    new_atm_num = 0
    # initialize new res number at 0
    new_res_num = 0
    # initialize prev_res to keep track of res number
    prev_res_num = None
    # initialize prev_new_res_num to keep track of the new numerotation
    prev_new_res_num = 0
    # Default chain for calculations
    new_chain_name = "A"
    # Use OrderedDict() objects to remember the order of residues insertion in the dictionary
    resnum_to_res = OrderedDict()
    with open(path, "r") as filin:
        for line in filin:
            if line.startswith("ATOM"):
                new_atm_num += 1
                old_res_num = int(line[22:26])
                atm_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_name = line[21:22].strip()
                if chain_name == "":
                    chain_name = "A" 
                alter_pos = line[16:17]
                if (alter_pos == " " or alter_pos == "A") and chain_name == chain:
                    try:
                        three_to_one[res_name]
                    except KeyError:
                        sys.exit(f"Wrong input PDB file format: {path}.\nPosition {old_res_num}: \"{res_name}\" is not an amino acid.")
                    if prev_res_num != old_res_num:
                        prev_res_num = old_res_num
                        new_res_num += 1
                    # Check that the previous residu had a CA, otherwise remove it and reset good number
                    if new_res_num != prev_new_res_num and prev_new_res_num != 0:
                        atom_names = [atom_infos[3] for atom_infos in resnum_to_res[prev_new_res_num] if len(atom_infos) > 1]
                        if "CA" not in atom_names:
                            print(f"\nWARNING: Skipped unsolved residue {prev_new_res_num} of {path}. No CA atom detected.")
                            new_res_num = prev_new_res_num
                            new_atm_num -= len(resnum_to_res[prev_new_res_num])
                            resnum_to_res[prev_new_res_num] = []
                    # Store all atom names of each residu number
                    try:
                        resnum_to_res[new_res_num].append([line, new_atm_num, new_res_num, atm_name, new_chain_name])
                    except KeyError:
                        resnum_to_res[new_res_num] = [[line, new_atm_num, new_res_num, atm_name, new_chain_name]]
                    finally:
                        ori_res_num_and_chain[new_res_num] = (old_res_num, chain_name)
                    prev_new_res_num = new_res_num
                elif chain_name != chain:
                    sys.exit(f"Error: the chain '{chain}' was not found in {path}. Check the chain requested. Also, please make sure there is only 1 model and 1 chain only in the PDB.")
                else:
                    continue
            elif line.startswith("TER") and chain_name == chain:
                resnum_to_res[new_res_num].append(["TER"])
                ori_res_num_and_chain[new_res_num] + (line,)
                break
            elif line.startswith("END"):
                break

    # Write the clean PDB file
    # If asked, write to new_path, otherwise, overwrite input PDB file
    if new_path:
        path = new_path
    with open(path, "w") as filout:
        for _, atoms_infos in resnum_to_res.items():
            for atom_infos in atoms_infos:
                # End of chain detected, but for calculations we consider all chains as one "A"
                if len(atom_infos) == 1:
                    continue
                else:
                    atom = atom_infos[0]
                    new_atm_num = atom_infos[1]
                    new_res_num = atom_infos[2]
                    new_chain_name = atom_infos[4]
                    filout.write(reformat_atom(atom, new_atm_num, new_res_num, new_chain=new_chain_name).strip() + "\n")
    return len(resnum_to_res)
                



def renum_ori_pdb_resnum(path, ori_res_num_and_chain, new_path=None):
    """
    Write the original residue number of each residue in the PDB file.

    Args:
        - path (str): path to the PDB file
        - ori_res_num_and_chain (dict): mapping of new to old residues numerotation

    Returns:
        None
    """
    # Get file in memory because we will override the file
    with open(path, "r") as filin:
        atoms = filin.readlines()

    # Do not overwrite the file if new_path is defined
    if new_path is not None:
        path = new_path
    new_n = 0  # initialize atom number at 0
    with open(path, "w") as filout:
        for atom in atoms:
            if atom.startswith("ATOM"):
                new_n += 1
                new_res = int(atom[22:26])
                ori_res = ori_res_num_and_chain[new_res][0]
                ori_chain = ori_res_num_and_chain[new_res][1]
                filout.write(reformat_atom(atom, new_n, ori_res, new_chain=ori_chain).strip() + "\n")
                # This PDB has more than one chain so we need to write down the "TER" line
                # at the end of the actual chain
                if len(ori_res_num_and_chain[new_res]) > 2:
                    filout.write(ori_res_num_and_chain[new_res][2].strip() + "\n")
            else:
                filout.write(atom)


def renum_ori_pdb_resnum_kpax(path, ori_res_num_and_chain1, ori_res_num_and_chain2):
    """
    Write the original residue number of each residue for both query and target in the KPAX result output.

    Args:
        - path (str): path to the PDB file
        - ori_res_num_and_chain (dict): mapping of new to old residues numerotation

    Returns:
        None
    """
    # Get file in memory because we will override the file
    with open(path, "r") as filin:
        atoms = filin.readlines()

    flag = False
    new_n = 0  # initialize atom number at 0
    with open(path, "w") as filout:
        for atom in atoms:
            # Detect the beginning of the second chain
            if atom.startswith("TER"):
                flag = True
                filout.write(atom)
            # Write the first chain
            elif atom.startswith("ATOM") and not flag:
                new_n += 1
                new_res = int(atom[22:26])
                ori_res = ori_res_num_and_chain1[new_res][0]
                ori_chain = ori_res_num_and_chain1[new_res][1]
                filout.write(reformat_atom(atom, new_n, ori_res, new_chain=ori_chain).strip() + "\n")
                # This PDB has ore than one chain so we need to write down the "TER" line
                # at the end of the actual chain
                if len(ori_res_num_and_chain1[new_res]) > 2:
                    filout.write(ori_res_num_and_chain1[new_res][2].strip() + "\n")
            # Write the second chain
            elif atom.startswith("ATOM") and flag:
                new_n += 1
                new_res = int(atom[22:26])
                ori_res = ori_res_num_and_chain2[new_res][0]
                ori_chain = ori_res_num_and_chain2[new_res][1]
                filout.write(reformat_atom(atom, new_n, ori_res, new_chain=ori_chain).strip() + "\n")
                # This PDB has ore than one chain so we need to write down the "TER" line
                # at the end of the actual chain
                if len(ori_res_num_and_chain2[new_res]) > 2:
                    filout.write(ori_res_num_and_chain2[new_res][2].strip() + "\n")
            # Write the rest of the lines of the PDB files
            else:
                filout.write(atom)


if __name__ == "__main__":
    pass
