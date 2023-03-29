"""
This module implements the Protein class.
"""

import os
import re
import shutil
import subprocess
import sys
import signal
import numpy as np
from collections import OrderedDict
from prody import *
confProDy(verbosity='none')

import src.utils as utils

# When icarus.py is executed, this equals: /path/to/icarus
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWORD = os.path.join(PROJECT_DIR, "bin", "SWORD")
WORK_DIR = os.path.join(os.getcwd(), "icarus_output")
PDB_CLEAN_DIR = os.path.join(WORK_DIR, "PDBs_Clean")
PDB_STAND_DIR = os.path.join(WORK_DIR, "PDBs_Stand")
TMP_DIR = utils.TMP_DIR

class Protein:
    """
    Used to manage data related to a protein. Initialized by giving the path to
    the protein structure.
    Protein object can either undergo protein peeling or not.
    """
    seg_size = 15  # min PU size by default

    def __init__(self, path, peel=True, ori_path=None, min_len_pu=15):
        """
        Creates a Protein object instance based on a 3D structural file.

        Attributes:
            - name (str): name of the protein.
            - path (str): path to structural file of the protein.
            - length (int): length of the protein.
            - seq (str): 1D sequence of the protein.
            - PUs (list): list of list of PU instances. Each sublist correspond to a partitioning level.
            - success (bool): Was the peeling successfull ?

        Arguments :
            - path (str) : path to the pdb file of the structure
        """

        __slots__ = ["name", "path", "length", "residues_num", "seq", "start", "end", "PUs_per_level", "PUs"]

        self.name = None
        self.path = None
        self.ori_path = ori_path
        self.length = None
        self.residues_num = None
        self.seq = None
        self.start = None
        self.end = None
        self.PUs_per_level = {}

        

        if not os.path.exists(path):
            raise OSError("'{}': No such file or directory".format(path))
        if peel:
            os.makedirs(PDB_STAND_DIR, exist_ok=True)
            os.makedirs(PDB_CLEAN_DIR, exist_ok=True)
            shutil.copy2(path, WORK_DIR)
            # get the file name, remove rest of path
            filename = os.path.basename(path)
            self.name = filename
            dest_dir = os.path.join(PDB_CLEAN_DIR, self.name)
            # Regenerate SWORD/Peeling outputs at each run
            # to avoid strange surprises for user
            if os.path.exists(dest_dir):
                shutil.rmtree(dest_dir)
            subprocess.run([SWORD, "-i", self.name, "--dir", WORK_DIR, "--max", "15",
                            "--size", str(min_len_pu)], capture_output=False)
            if not os.path.exists(dest_dir):
                with os.scandir(os.path.join(PDB_CLEAN_DIR)) as dir:
                    for file in dir:
                        if self.name == file:
                            os.rename(os.path.join(PDB_CLEAN_DIR, file), dest_dir)
                            break  # file name is change, no need to go further
            self.path = os.path.join(PDB_STAND_DIR, self.name)
            shutil.move(os.path.join(WORK_DIR, self.name), self.path)
        else:
            self.path = path
            # get last part of path i.e name
            self.name = os.path.basename(os.path.splitext(path)[0])
        self._set_length()
        if peel:
            PUs = self._set_PUs()
            for level, pus in enumerate(PUs):
                id_pu = 1
                self.PUs_per_level[level] = {}
                for pu in pus:
                    self.PUs_per_level[level][id_pu] = PU(pu)
                    id_pu += 1
            del PUs
        else:
            self.PUs_per_level = {}

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def _set_length(self):
        """
        Sets the length variable of a prot or PU object to the number of residues
        in the pdb file associated with the object.
        Also save the residues numbers in self.residues_num

        Args:
            - path (str): path to the pdb file

        Returns:
            - None
        """
        residues = set()
        with open(self.path, "r") as filin:
            for line in filin:
                if line.startswith("ATOM"):
                    residues.add(int(line[22:26]))
        # count unique elements of res_nb
        self.length = len(residues)
        self.residues_num = sorted(residues)

    def _set_PUs(self):
        """
        Given the path to the protein peeling output, return a list of list
        of path to the protein units. Each sublist correspond to a segmentation
        level.

        Args:
            - filename (str): name of the subdirectory containing the protein
                             peeling output

        Returns:
            - PUs (list): List of sublist of paths to Protein Units
        """
        pdbs = []
        # Get peeling output pdb filenames
        for file in os.listdir(os.path.join(PDB_CLEAN_DIR, self.name, "Peeling")):
            if file[-4:] == ".pdb":
                pdbs.append(file)
        # The following complex expression simply does one thing:
        # sort peeling output filenames by segmentation level AND PU order:
        # Ex: [ d1qasa2.ent_2_1_17.pdb, d1qasa2.ent_2_18_44.pdb, d1qasa2.ent_2_45_59.pdb, d1qasa2.ent_2_60_96.pdb ]
        sorted_pdbs = sorted(pdbs, key=lambda x: (int(x.split("_")[-3]),
                                                  int(x.split("_")[-2])))
        PUs = []
        seg_levels = []
        # Build sublist of PUs for each segmentation level
        for pdb in sorted_pdbs:
            # first number encountered is the segmentation level
            search = re.search(r"_(\d+)_\d+_\d+", pdb)
            seg_levels.append(int(search.group(1)))
        for level in set(seg_levels):
            # add empty sub list for each segmentation level
            PUs.append([])
        # Fill sublists with full paths of PUs
        for i, level in enumerate(seg_levels):
            PUs[level - 1].append(os.path.join(PDB_CLEAN_DIR, self.name, "Peeling", sorted_pdbs[i]))
        return PUs

    def set_1d_seq(self):
        """
        Set the 1d sequence of the protein from its PDB file

        Attributes:
            - set self.seq to seq (str): the 1d sequence of the protein.
        """
        three_to_one_d = {
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
        seq = ""
        seen = set()
        with open(self.path, "r") as filin:
            for line in filin:
                if line.startswith("ATOM"):
                    if int(line[22:26]) in seen:
                        continue
                    seen.add(int(line[22:26]))
                    try:
                        # add 1 char name of res to seq
                        seq += three_to_one_d[line[17:20]]
                    except KeyError:
                        seq += "X"
            self.seq = seq

    def reformat_struct(path):
        """
        Reset atoms and residues number of a newly merged pdb file hence ordering
        the records properly according to the alignments.
        Changes are applied directly in the file indicated by path.

        Args:
            - path (str): path to the structural file

        Returns:
            None
        """
        # Get in memory because we will override the file
        with open(path, "r") as filin:
            atoms = filin.readlines()  # fetch atoms in old pdb

        new_n = 0  # initialize new atom number at 0
        new_res = 0  # initialize new res number at 0
        prev_res = 0  # initialize prev_res to keep track of res number
        with open(path, "w") as filout:
            for atom in atoms:
                if atom.startswith("ATOM"):
                    new_n += 1
                    old_res = int(atom[22:26])
                    if prev_res != old_res:
                        prev_res = old_res
                        new_res += 1
                    filout.write(utils.reformat_atom(atom, new_n, new_res))
                elif atom.startswith("TER"):
                    filout.write("TER\n")
                elif atom.startswith("HETATM"):
                    continue
                else:
                    filout.write(atom)


class PU(Protein):
    """
    The PU class inherits Protein's class attributes and functions.
    Used to store and manipulate data associated with a Protein Unit.
    """

    def __init__(self, path):
        """
        Creates a PU object, very similar to a Protein objet,
        the main difference being it does not have a PU list, and PUs have
        delimitations.
        """
        if not os.path.exists(path):
            raise OSError("'{}': No such file or directory".format(path))
        self.path = path
        self.name = path.split("/")
        self.name = "{}".format(self.name[-1][:-4])
        # Get PU start and end from Protein Peeling output
        search = re.search(r"_(\d+)_(\d+)$", self.name)
        if search:
            self.start = int(search.group(1))
            self.end = int(search.group(2))

        self._set_length()  # set length
        self.set_1d_seq()  # set 1D sequence