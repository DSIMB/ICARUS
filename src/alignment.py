"""
This module implements the Alignment class.
"""

import os
import re
import subprocess
import signal
import sys
import shutil
import numpy as np
from pathlib import Path

import src.utils as utils

from .protein import PU

# When icarus.py is executed, this equals: /path/to/icarus
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KPAX = os.path.join(PROJECT_DIR, "bin", "kpax", "bin", "kpax5.1.3.x64")
WORK_DIR = os.path.join(os.getcwd(), "icarus_output")
RESULTS_DIR = os.path.join(WORK_DIR, "results")
TMP_DIR = utils.TMP_DIR

class Alignment:
    """
    The purpose of this class in storing paths and informations relative to an
    alignment realised with the TM-align program.
    """

    def __init__(self, query, target, opt_prune, seed_alignment, save_output=False, keep_ori_resnum=False):
        """
        alignment instance creator

        Args:
            - query (Protein or PU): The query protein or PU to align to the target protein
            - target (Protein or PU): The target protein on which the query will be aligned
            - opt_prune (float): TM-score threshold to prune the alignment
            - seed_alignment (bool): if True, the TM-align alignments will be fixed with a seed alignment
            - save_output (bool): if True, the TM-align output PDB will be saved in the icarus_output folder

        Attributes:
            - query (Protein or PU): passed query
            - target (Protein or PU): passed target
            - opt_prune (float): TM-score threshold to prune the alignment
            - seed_alignment (bool): if True, the TM-align alignments will be fixed with a seed alignment
            - path (str): path to the generated alignement folder
            - score (float): TM-score associated with the alignment
            - success (bool): True if the alignment was successful
            - new_query (str): string corresponding to the modifieds ATOMs of the query file after being subject to
                               rotatory constrains during alignment
            - new_target (str): string corresponding to the modifieds ATOMs of the target file after being subject to
                                rotatory constrains during alignment
            - updated_target (str): path to updated target file, returned by the update_target function
        """

        # slots allow faster attribute access and space savings in memory
        __slots__ = ["query", "target", "path", "score", "new_query", "new_target", "updated_target", "__dict__"]

        self.query = query
        self.target = target
        self.opt_prune = opt_prune
        self.seed_alignment = seed_alignment
        self.success = None
        self.path = None  # initialized by _align
        self.score = None  # initialized by _align
        self.all_aligned = None  # initialized by _get_match
        self.new_query = None  # initialized by _get_structures()
        self.new_target = None  # initialized by _get_structures()
        self.new_query_dict = None  # initialized by _get_structures()
        self.updated_target = None  # initialized by _update_target()
        self.save_output = save_output
        self.keep_ori_resnum = keep_ori_resnum
        self._align(query, target, opt_prune)
        if self.success:
            self._get_match()
            self._get_structures()
            self.name = self.path.split("/")[-2]
            self._update_target()

    def __repr__(self):
        """
        objects are represented by the name of the alignment folder.
        """
        return self.name

    def _align(self, query, target, opt_prune):
        """
        Align a query pdb file onto a target pdb file using KPAX
        and set the associated TM-score normalized by the shortest sequence.
        KPAX normalizes by the shortest sequence by default.
        KPAX aligns the 2nd file on the 1st (which stays rigid)

        Args:
            - query, str: path to a query pdb file to superimpose on target file
            - target, str: path to a target pdb file

        Returns:
            - None

        Attributes:
            - set self.success
            - set self.path
           
            - set self.score
        """
        ali = os.path.join(TMP_DIR, f"{query.name}-on-{target.name}")
        os.makedirs(ali, exist_ok=True)

        alignment = f"{ali}/kpax_results/{query.name}_{target.name}_flex.pdb"
        # Superimpose query to target with KPAX (1st structure is rigid and kpax aligns the 2nd onto the 1st one).
        opt = ["-flex", "-nosubdirs", "-nohex", "-novmd", "-nojmol", "-nomatrix",
               "-nosse", "-nofasta", "-nopir", "-nokrmsd", "-noprofit", "-nohits",
               "-notops", "-norank", "-norainbow", "-pdb", target.path, query.path]
        output = subprocess.run([KPAX] + opt, cwd=ali, capture_output=True)
        res = output.stdout.decode("utf-8").split("\n")
        if self.save_output:
            kpax_result_filename = query.name + "_on_" + target.name + "_kpax.pdb"
            kpax_result_path = os.path.join(RESULTS_DIR, query.name + "_on_" + target.name, "result_PDBs", kpax_result_filename)
            # Create the directory because in special cases it is not yet created
            os.makedirs(os.path.dirname(kpax_result_path), exist_ok=True)
            shutil.copy2(f"{alignment}", kpax_result_path)
            self.kpax_result_path = kpax_result_path
        # to_remove = ("*pml", "*lig", "*all", "*pdb")
        # files = []
        # for files_pattern in to_remove:
        #     for files in Path(ali).glob(files_pattern):
        #         files.unlink()
        tm_score = None
        for line in res:
            # Get TM_score normalized by user-input (shortest length)
            tm_score_found = re.search(r"^\s+1\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+).*$", line)
            if tm_score_found:
                tm_score = round(float(tm_score_found.group(1)), 3)
                break
        if tm_score < opt_prune or tm_score is None:
            self.success = False
        else:
            self.success = True
            self.path = ali+"/kpax_results"
            self.score = tm_score

    def _get_core_and_all_matching_residues(self):
        """
        Get the core matching residues between the query and the target from the KPAX output,
        as well as all the aligned residues.
        Also get the target sequence involved in the alignment with the PU.

        Args:
            - aa_dict (dict): dictionary of amino acids

        Returns:
            - core_pu_pos (list): list of the core PU matching residues in the query
            - core_target_pos (list): list of the core matching residues in the target
            - all_target_pos (list): list of all the matching residues in the target
        """
        core_pu_pos = []
        core_target_pos = []
        target_pos = []
        pu_pos = [] # CA only
        target_pos = []  # CA only
        pu_seq = ""
        target_seq = ""
        with open(f"{self.path}/{self.query.name}_{self.target.name}_flex.wpairs", "r") as filin:
            for line in filin:
                line = line.strip()
                if not line.startswith("#"):
                    target, query, star = line.split("|")
                    # Aligned residues
                    if star.strip() == "*":
                        q_chain, q_resnum, q_resname = query.strip().split(":")
                        t_chain, t_resnum, t_resname = target.strip().split(":")
                        core_pu_pos.append(int(q_resnum))
                        core_target_pos.append(int(t_resnum))
                        pu_pos.append(int(q_resnum))
                        target_pos.append(int(t_resnum))
                        target_seq += t_resname
                        pu_seq += q_resname
                    # Not aligned residues
                    # else:
                    #     if query.strip() != "-":
                    #         q_chain, q_resnum, q_resname = query.strip().split(":")
                    #         pu_pos.append(int(q_resnum))
                    #     if target.strip() != "-":
                    #         t_chain, t_resnum, t_resname = target.strip().split(":")
                    #         target_pos.append(int(t_resnum))
        return core_pu_pos, core_target_pos, target_pos, pu_pos, pu_seq, target_seq

    def _get_match(self):
        """
        Create and assign to self.all_aligned a list of tuple representing the
        match between query and target.
        Create and assign to self.all_aligned a dict of lists representing
        the best matches of the alignment.
        i.e tuple (1, 4) res 1 of query matched with res 4 of target.

        Attributes:
            - set self.all_aligned
        """
        # Get the best core matching residues between query and target according to KPAX
        # and get all the atom coordinates of the *best* aligned residues in the query and target
        core_pu_pos, core_target_pos, target_pos, pu_pos, pu_seq, target_seq = self._get_core_and_all_matching_residues() 
        self.all_aligned = {"core_pu_aligned_positions": core_pu_pos,
                            "core_target_aligned_positions": core_target_pos,
                            "all_pu_aligned_positions": pu_pos,
                            "all_target_aligned_positions": target_pos,
                            "pu_seq": "".join(pu_seq),
                            "target_seq": "".join(target_seq)}

    def _get_structures(self):
        """
        Fetch structural conformation of PU or target that have been subject to
        rotational translation during KPAX. The query stays rigid.
        Stores informations in Alignment.new_query and Alignment.new_target
        variables.

        KPAX flexible alignment means that the query is flexibly aligned onto the rigid target,
        in the output files we get this naming convention:
            <query_name>_<target_name>_flex.pdb and the rigid query file <query_name>_query.pdb

        Attributes:
            - set self.new_query
            - set self.new_target
        """
        new_query = ""
        new_query_dict = {}
        new_target = ""
        with open(f"{self.path}/{self.query.name}_{self.target.name}_flex.pdb", 'r') as filin:
            for line in filin:
                line = line.strip()
                if line.startswith("ATOM"):
                    # This is to have the PU as PDB to align
                    new_query += line[:-1] + ' ' * (81 - len(line)) + "\n"
                    # This is to have the correspondance between the residue number and its PDB line
                    try:
                        new_query_dict[int(line[22:26].strip())].append(line[:-1] + "\n")
                    except KeyError:
                        new_query_dict[int(line[22:26].strip())] = [line[:-1] + "\n"]
        with open(f"{self.path}/{self.target.name}_query.pdb", 'r') as filin:
            for line in filin:
                line = line.strip()
                if line.startswith("ATOM"):
                    # This is to have the PU as PDB to align
                    new_target += line[:-1] + ' ' * (81 - len(line)) + "\n"
        self.new_query = new_query
        self.new_query_dict = new_query_dict
        self.new_target = new_target

    def _update_target(self):
        """
        Update the target structure by removing the residues that matched the
        query structure.

        Args:
            - self : Alignment class instance

        Returns:
            - None

        Attributes:
            - Set self.updated_target as the path to the updated structure.
        """
        res_to_del = [res for res in self.all_aligned["core_target_aligned_positions"]]
        old_pdb = self.target
        new_pdb = os.path.join(TMP_DIR, utils.get_random_name() + ".pdb")
        with open(old_pdb.path, "r") as filin, open(new_pdb, 'w') as filout:
            for atom in filin:
                res_nb = int(atom[22:26])
                if res_nb not in res_to_del:
                    filout.write(atom)
        self.updated_target = new_pdb

    def get_new_alis(self, id_pu, alis):
        """
        Re-align PU of alignment which are not already aligned to the
        updated target.

        Args:
            - id_pu (int): id of the PU to remove because already aligned
            - alis (list of Alignment): list of new alignments of PU against a target protein

        Returns:
            - alis (list of Alignment): Realined PUs
        """
        if len(alis) == 1:
            return []  # get index of ali
        new_alis = alis.copy()
        # new_target as PU to realise alignment
        new_target = PU(self.updated_target)
        del new_alis[id_pu]  # remove aligned PU
        for i, ali in new_alis.items():
            new_alis[i] = Alignment(ali.query, new_target, self.opt_prune, self.seed_alignment)
        return new_alis
        
    def smooth_positions(self, min_contiguous=5, n_gaps=2):
        """
        Remove alignment positions that are isolated.
        --> if less than `min_contiguous` contiguous positions or separated by `n_gaps` or more positions

        Args:
            - min_contiguous (int): minimum number of contiguous positions
            - n_gaps (int): maximum number of gaps authorized between contiguous positions

        Returns:
            - None

        Attributes:
            - The best_path but with smoothed positions:
                self.all_aligned["all_target_aligned_positions"]
                self.all_aligned["all_pu_aligned_positions"]
                self.all_aligned["pu_seq"]
                self.all_aligned["target_seq"] 
        """
        # Best path in the graph == best PUs aligned against the target
        smoothed_target_positions = []
        smoothed_pu_positions = []
        smoothed_target_seq = ""
        smoothed_pu_seq = ""
        gaps_vect = np.array([], dtype=int)
        for i in range(len(self.all_aligned["all_target_aligned_positions"]) - 1):
            t_pos = self.all_aligned["all_target_aligned_positions"][i]
            t_next_pos = self.all_aligned["all_target_aligned_positions"][i+1]
            # Vector of nb of gaps between each positions: [1, 2, 4] --> [0, 1]
            gaps_vect = np.append(gaps_vect, (t_next_pos - t_pos - 1))
        # Get indexes of where there are gaps > `n_gaps` positions
        idx_of_large_gaps = np.argwhere(gaps_vect > n_gaps).flatten()
        # Smooth only when there are gaps, otherwise keep the original positions
        if len(idx_of_large_gaps) != 0:
            # Take care of the first gap first
            t_slice = self.all_aligned["all_target_aligned_positions"][:idx_of_large_gaps[0] + 1]
            # Keep only groups of positions of size `min_contiguous` or more
            if len(t_slice) >= min_contiguous:
                smoothed_target_positions += t_slice
            if len(idx_of_large_gaps) > 1:
                # Take care of middle range gaps
                for j in range(len(idx_of_large_gaps) - 1):
                    idx = idx_of_large_gaps[j]
                    next_idx = idx_of_large_gaps[j+1]
                    t_slice = self.all_aligned["all_target_aligned_positions"][idx+1:next_idx + 1]
                    if len(t_slice) >= min_contiguous:
                        smoothed_target_positions += t_slice
            # Take care of the last gap
            # If the gap occurs between the `min_contiguousth` before last and last positions
            # then skip since the group won't be `min_contiguous` long
            if idx_of_large_gaps[-1] != len(self.all_aligned["all_target_aligned_positions"]) - min_contiguous:
                t_slice = self.all_aligned["all_target_aligned_positions"][idx_of_large_gaps[-1]+1:]
                if len(t_slice) >= min_contiguous:
                    smoothed_target_positions += t_slice
            # Smooth the corresponding PUs and sequences positions
            for i, pos in enumerate(self.all_aligned["all_target_aligned_positions"]):
                if pos in smoothed_target_positions:
                    smoothed_pu_positions.append(self.all_aligned["all_pu_aligned_positions"][i])
                    smoothed_target_seq += self.all_aligned["target_seq"][i]
                    smoothed_pu_seq += self.all_aligned["pu_seq"][i]
            self.all_aligned["smoothed_all_target_aligned_positions"] = smoothed_target_positions
            self.all_aligned["smoothed_all_pu_aligned_positions"] = smoothed_pu_positions
            self.all_aligned["smoothed_pu_seq"] = smoothed_pu_seq
            self.all_aligned["smoothed_target_seq"] = smoothed_target_seq
            self.all_aligned["gaps_vect"] = gaps_vect
        else:
            self.all_aligned["smoothed_all_target_aligned_positions"] = self.all_aligned["all_target_aligned_positions"]
            self.all_aligned["smoothed_all_pu_aligned_positions"] = self.all_aligned["all_pu_aligned_positions"]
            self.all_aligned["smoothed_pu_seq"] = self.all_aligned["pu_seq"]
            self.all_aligned["smoothed_target_seq"] = self.all_aligned["target_seq"]
            self.all_aligned["gaps_vect"] = gaps_vect


    def get_smoothed_target_positions(self, min_contiguous=5, n_gaps=2):
        """
        Remove only **target** positions that are isolated.
        --> if less than `min_contiguous` contiguous positions or separated by `n_gaps` or more positions

        Args:
            - min_contiguous (int): minimum number of contiguous positions
            - n_gaps (int): maximum number of gaps authorized between contiguous positions

        Returns:
            - self.all_aligned["smoothed_core_target_aligned_positions"] (list): list of smoothed target positions
                or the original list if no smoothing was required
        """
        # Best path in the graph == best PUs aligned against the target
        smoothed_target_positions = []
        gaps_vect = np.array([])
        for i in range(len(self.all_aligned["core_target_aligned_positions"]) - 1):
            t_pos = self.all_aligned["core_target_aligned_positions"][i]
            t_next_pos = self.all_aligned["core_target_aligned_positions"][i+1]
            # Vector of gaps between positions: [1, 2, 4] --> [0, 1]
            gaps_vect = np.append(gaps_vect, (t_next_pos - t_pos - 1))
        # Get indexes of where there are gaps > `n_gaps` positions
        idx_of_large_gaps = np.argwhere(gaps_vect > n_gaps).flatten()
        # Smooth only when there are gaps, otherwise keep the original positions
        if len(idx_of_large_gaps) != 0:
            # Take care of the first gap first
            t_slice = self.all_aligned["core_target_aligned_positions"][:idx_of_large_gaps[0] + 1]
            # Keep only groups of positions of size `min_contiguous` or more
            if len(t_slice) >= min_contiguous:
                smoothed_target_positions += t_slice
            if len(idx_of_large_gaps) > 1:
                # Take care of middle range gaps
                for j in range(len(idx_of_large_gaps) - 1):
                    idx = idx_of_large_gaps[j]
                    next_idx = idx_of_large_gaps[j+1]
                    t_slice = self.all_aligned["core_target_aligned_positions"][idx+1:next_idx + 1]
                    if len(t_slice) >= min_contiguous:
                        smoothed_target_positions += t_slice
            # Take care of the last gap
            # If the gap occurs between the `min_contiguous` before last and last positions
            # then skip since the group won't be `min_contiguous` long
            if idx_of_large_gaps[-1] != len(self.all_aligned["core_target_aligned_positions"]) - min_contiguous:
                t_slice = self.all_aligned["core_target_aligned_positions"][idx_of_large_gaps[-1]+1:]
                if len(t_slice) >= min_contiguous:
                    smoothed_target_positions += t_slice
            self.all_aligned["smoothed_core_target_aligned_positions"] = smoothed_target_positions
            # Case when there are gaps, but when removing the gaps, there is nothing left.
            # Keep the core positions anyway.
            if len(smoothed_target_positions) == 0:
                self.all_aligned["smoothed_core_target_aligned_positions"] = self.all_aligned["core_target_aligned_positions"]
        else:
            self.all_aligned["smoothed_core_target_aligned_positions"] = self.all_aligned["core_target_aligned_positions"]
        return self.all_aligned["smoothed_core_target_aligned_positions"]

    @staticmethod
    def multiple_alignment(query_prot, target_prot, opt_prune, seed_alignment, seg_level=0):
        """
        Align each PU from the query_prot Protein object to the target_prot protein.
        PUs from the segmentation level seg_level are align, default 0.

        Args:
            - query_prot (Protein): Protein object which PUs, will be aligned.
            - target_prot (Protein): Protein object where PUs will be aligned to.

        Returns:
            - alis (list of Alignment): List of alignment objects resulting from PU vs prot alignment
        """
        if seg_level > len(query_prot.PUs_per_level) or (seg_level < 0):
            print("Error: Requested segmentation level isn't available")
            return None
        if seg_level == 0:
            return Alignment(query_prot, target_prot, opt_prune, seed_alignment)  # full_length alignment

        alis = {}
        for id_pu, pu in query_prot.PUs_per_level[seg_level - 1].items():
            alis[id_pu] = Alignment(pu, target_prot, opt_prune, seed_alignment)
        return alis