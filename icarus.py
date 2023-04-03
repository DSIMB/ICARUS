#!/usr/bin/env python3
"""
    Icarus (flexIble struCtural Alignment based on pRotein UnitS)
    is a program for flexible structural alignement.
    This is the main script.
"""
import argparse
import gc
import math
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import textwrap
import time
from gc import get_referents
from multiprocessing import cpu_count
from types import FunctionType, ModuleType

import numpy as np

import src.graphpu as g
import src.protein as p
import src.utils as utils
from src.alignment import Alignment

# Start monitoring program runtime
start_time = time.time()

TMP_DIR = utils.TMP_DIR
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
WORK_DIR = os.path.join(os.getcwd(), "icarus_output")
RESULTS_DIR = os.path.join(WORK_DIR, "results")
PDB_CLEAN_DIR = os.path.join(WORK_DIR, "PDBs_Clean")
PDB_STAND_DIR = os.path.join(WORK_DIR, "PDBs_Stand")
GDT = os.path.join(PROJECT_DIR, "bin", "gdt2.pl")
for path in [WORK_DIR, RESULTS_DIR, TMP_DIR]:
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
# Keys=explore level ; values=number of PUs to consider in graph
INTERVALS = {0: [0], 1: [2, 3], 2: [4, 5], 3: [6], 4: [7], 5: [8]}
NB_PUS_2_EXPLORE_LEVEL = {0: 0, 2: 1, 3: 1, 4: 2, 5: 2, 6: 3, 7: 4, 8: 5}


def signal_handler(signal, handler):
    """
    Catch CTRL+C signal
    Clean the tmp directory before stopping
    """
    if os.path.exists(TMP_DIR):
        shutil.rmtree(TMP_DIR, ignore_errors=True)
        print("\nQuitting gracefully, bye !")
    sys.exit(0)

# Catch CTRL+C signal and quit gracefully by cleaning traces
signal.signal(signal.SIGINT, signal_handler)

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType
def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: ' +
                        str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


def run_gdt(query, target, min_len_p1_p2, opt_prune, ori_res_num_and_chain1, ori_res_num_and_chain2, save_output=False, keep_ori_resnum=False):
    """
    Runs gdt2.pl with input query and target files.
    As gdt requires file to be in current folder, files are copied in
    current folder then the copy in deleted once gdt execution is over.

    Args:
        - query (Protein or PU): query structure used as an input for gdt.pl.
        - target (Protein or PU): target structure used as an input for gdt.pl.
        - min_len_p1_p2 (int): minimum length of the query and target sequences.
        - opt_prune (int): optional pruning by TM-score threshold.
        - save_output (bool): if True, gdt2.pl output is saved in a file.
        - keep_ori_resnum (bool): if True, original residue numbers are kept in the KPAX PDB output.

    Returns:
        - score (float): TM-score normalized by length of shortest structure.
    """
    mode = 0
    # KPAX before calculating the scores with gdt2.pl
    # use the original query and target residues numbers
    ali = Alignment(query, target, opt_prune, save_output, keep_ori_resnum)
    with open(query.name, "w") as filout:
        filout.write(ali.new_query)
    with open(target.name, "w") as filout:
        filout.write(ali.new_target)
    command = f"{GDT} -pdb '{query.name} {target.name}' -mode {mode} -len {min_len_p1_p2}"
    cmd_args = shlex.split(command)
    # execution of gdt2.pl
    output = subprocess.run(cmd_args, capture_output=True, check=True)
    if os.path.exists(query.name):
        os.remove(f"{query.name}")
    if os.path.exists(target.name):
        os.remove(f"{target.name}")
    output = output.stdout.decode("utf-8")
    if save_output:
        # For cases when we want to launch KPAX and keep the original residue numbers and chain names
        if keep_ori_resnum:
            utils.renum_ori_pdb_resnum_kpax(ali.kpax_result_path, ori_res_num_and_chain1, ori_res_num_and_chain2)
        base_path = os.path.join(RESULTS_DIR, query.name + "_on_" + target.name)
        os.makedirs(base_path, exist_ok=True)
        with open(f"{base_path}/gdt_{query.name}_vs_{target.name}.txt", "w") as filout:
            filout.write(output)
    output = output.split("\n")
    tm_scores = []
    for line in output:  # search output for TM-scores returned by gdt
        if ("TM-score" in line) and ("user length" in line):
            # retrieves scores as floats
            tm_scores.append(float(line[11:16]))
    # to return score normalized by short chain
    if target.length < query.length:
        return tm_scores[0]
    return tm_scores[1]


def main(p1,
         p2,
         min_len_p1_p2,
         exploration_level,
         exploration_level_p1,
         exploration_level_p2,
         opt_prune,
         smoothed_pu_output,
         ori_res_num_and_chain1,
         ori_res_num_and_chain2,
         sequential,
         nb_cpu=0,
         verbose=False):
    """
    For two proteins, creates all graphs associated with all exploration
    levels. Print a summary of the results and return a string corresponding
    to this summary

    Args:
        - p1 (Protein): One of the two protein to align.
        - p2 (Protein): The other protein to align.
        - min_len_p1_p2 (int): Length of the smallest sequence between query and target sequences
        - exploration_level (int): Exploration level requested by the user
        - exploration_level_p1 (int): If not None, represents the new maximum exploration level
                                      for the protein 1 because the original exploration level
                                      was not adapted (too high)
        - exploration_level_p2 (int): If not None, represents the new maximum exploration level
                                      for the protein 2 because the original exploration level
                                      was not adapted (too high)
        - opt_prune (int): optional pruning by TM-score threshold.
        - smoothed_pu_output (bool): if True, the output of the PU is smoothed
        - ori_res_num_and_chain1 (dict): original number of residues and chain(s) of protein 1
        - ori_res_num_and_chain2 (dict): original number of residues and chain(s) of protein 2
        - nb_cpu (int): number of CPU to use. If 0, all available CPUs are used.
        - sequential (bool): if True, do a sequential alignment: keep paths with consecutive PUs only.
        - verbose (bool): If true, print all solutions in the final output, else print the optimal only.
    """
    results = []
    scores = []
    max_expl_lvl = ""
    max_expl_lvl_p1 = ""
    solutions_cnt = 1
    textual_alignment_p1_vs_p2 = ""
    textual_alignment_p2_vs_p1 = ""
    # Print maximum exploration levels explored and corresponding number of PUs
    # Case when we can use the exploration level requested by user
    if exploration_level_p1 is None and exploration_level_p2 is None:
        if len(INTERVALS[exploration_level]) > 1:
            max_expl_lvl = " or ".join(str(i) for i in INTERVALS[exploration_level])
        else:
            max_expl_lvl = INTERVALS[exploration_level]
        results.append(f"\nINFO: Maximum exploration level {exploration_level} ({max_expl_lvl}) PUs\n")

    ###################### EXPLORE SOLUTION P1 vs. P2 ######################
    # The first protein is peeled into N PUs and aligned against second one

    # Print maximum exploration levels explored and corresponding number of PUs
    # Case when we need to use a new maximum exploration level because the one
    # requested by user is not adapted
    if exploration_level_p1 is not None:
        if exploration_level_p1 == 0:
            max_expl_lvl_p1 = "no"
        elif len(INTERVALS[exploration_level_p1]) > 1:  
            max_expl_lvl_p1 = " or ".join(str(i) for i in INTERVALS[exploration_level_p1])
        else:
            max_expl_lvl_p1 = INTERVALS[exploration_level_p1][0]
        results.append(f"\nINFO: The maximum exploration level possible for '{p1.name}' is {exploration_level_p1} ({max_expl_lvl_p1} PUs)\n")

    results.append("                                                    SCORES")
    results.append("                                                    ======\n")
    results.append(f" *  {p1.name} against {p2.name}")
    print(f"\n\nAligning {p1.name} against {p2.name}", end="")
    if exploration_level_p1 == 0:  # Peeling coudn't find any PU on the protein, we only do a KPAX
        start_t = time.time()
        results.append(f" └── level 0 | 0 PUs (plain KPAX): {run_gdt(p1, p2, min_len_p1_p2, opt_prune, ori_res_num_and_chain1, ori_res_num_and_chain2, keep_ori_resnum=True)}")
        textual_alignment_p1_vs_p2 = ""
        print("\n")
    else:
        results.append(f" ├── level 0 | 0 PUs (plain KPAX): {run_gdt(p1, p2, min_len_p1_p2, opt_prune, ori_res_num_and_chain1, ori_res_num_and_chain2, keep_ori_resnum=True)}")
        start_t = time.time()
        for level in range(g.GraphPU.max_seg_level_p1):
            nb_pu_at_level = len(p1.PUs_per_level[level])
            expl_level = NB_PUS_2_EXPLORE_LEVEL[nb_pu_at_level]
            print(f"\n\nSOLUTION {solutions_cnt}: {p1.name} [level {expl_level} => {nb_pu_at_level} PUs] vs {p2.name}")
            print(f"{56*'-'}")  # for aesthetic purposes
            graph = g.GraphPU(p1, p2, level + 1, expl_level, min_len_p1_p2, ori_res_num_and_chain1, ori_res_num_and_chain2, nb_cpu, opt_prune, smoothed_pu_output, sequential)
            if graph.succeeded:
                scores.append([graph.best_score, graph.pu_order_text, graph.best_ali, "1"])
                # aesthetics: last line to write
                if level == g.GraphPU.max_seg_level_p1 - 1:
                    results.append(f" └── level {expl_level} | {nb_pu_at_level} PUs: {graph.best_score}\n")
                else:
                    results.append(f" ├── level {expl_level} | {nb_pu_at_level} PUs: {graph.best_score}")
            solutions_cnt += 1
        if graph.succeeded:
            textual_alignment_p1_vs_p2 = draw_textual_alignment(p1, p2, graph, ori_res_num_and_chain1, smoothed_pu_output)
            results.append(textual_alignment_p1_vs_p2)
        del graph
        gc.collect()

    ###################### EXPLORE SOLUTION P2 vs. P1 ######################
    # The second protein is peeled into N PUs and aligned against first one

    # Print maximum exploration levels explored and corresponding number of PUs
    # Case when we need to use a new maximum exploration level because the one
    # requested by user is not adapted
    max_expl_lvl_p2 = ""
    print(f"\n\nAligning {p2.name} against {p1.name}", end="")
    if exploration_level_p2 is not None:
        if exploration_level_p2 == 0:
            max_expl_lvl_p2 = "no"
        elif len(INTERVALS[exploration_level_p2]) > 1:
            max_expl_lvl_p2 = " or ".join(str(i) for i in INTERVALS[exploration_level_p2])
        else:
            max_expl_lvl_p2 = INTERVALS[exploration_level_p2][0]
        results.append(f"\nINFO: The maximum exploration level possible for '{p1.name}' is {exploration_level_p2} ({max_expl_lvl_p2} PUs)\n")
    results.append("                                                    SCORES")
    results.append("                                                    ======\n")
    results.append(f" *  {p2.name} against {p1.name}")
    if exploration_level_p2 == 0:  # Peeling coudn't find any PU on the protein, we only do a KPAX
        start_t = time.time()
        results.append(f" └── level 0 | 0 PUs (plain KPAX): {run_gdt(p2, p1, min_len_p1_p2, opt_prune, ori_res_num_and_chain2, ori_res_num_and_chain1, keep_ori_resnum=True)}")
        textual_alignment_p2_vs_p1 = ""
    else:
        results.append(f" ├── level 0 | 0 PUs (plain KPAX): {run_gdt(p2, p1, min_len_p1_p2, opt_prune, ori_res_num_and_chain2, ori_res_num_and_chain1, keep_ori_resnum=True)}")
        start_t = time.time()
        for level in range(g.GraphPU.max_seg_level_p2):
            nb_pu_at_level = len(p2.PUs_per_level[level])
            expl_level = NB_PUS_2_EXPLORE_LEVEL[nb_pu_at_level]
            print(f"\n\nSOLUTION {solutions_cnt}: {p2.name} [level {expl_level} => {nb_pu_at_level} PUs] vs {p1.name}")
            print(f"{56*'-'}")  # for aesthetic purposes
            graph = g.GraphPU(p2, p1, level + 1, expl_level, min_len_p1_p2, ori_res_num_and_chain2, ori_res_num_and_chain1, nb_cpu, opt_prune, smoothed_pu_output, sequential)
            if graph.succeeded:
                scores.append([graph.best_score, graph.pu_order_text, graph.best_ali, "2"])
                # aesthetics: last line to write
                if level == g.GraphPU.max_seg_level_p2 - 1:
                    results.append(f" └── level {expl_level} | {nb_pu_at_level} PUs: {graph.best_score}\n")
                else:
                    results.append(f" ├── level {expl_level} | {nb_pu_at_level} PUs: {graph.best_score}")
            solutions_cnt += 1
        if graph.succeeded:
            textual_alignment_p2_vs_p1 = draw_textual_alignment(p2, p1, graph, ori_res_num_and_chain2, smoothed_pu_output)
            results.append(textual_alignment_p2_vs_p1)
    results += "\n"
    os.makedirs(base_path, exist_ok=True)

    print("\n\n")
    print("**************************************************************************************************************")
    print("*************************************************** RESULTS **************************************************")
    print("**************************************************************************************************************")

    # Get the best result scores
    best_scores, best_score_pos, which_alignment = get_best_scores(scores)
    # As textual output
    best, best_for_file, best_path = get_textual_scores(scores, best_score_pos,
                                                        which_alignment,
                                                        textual_alignment_p1_vs_p2,
                                                        textual_alignment_p2_vs_p1)

    solutions_txt = ""
    # Identify the PDB of the best solution(s)
    solutions_txt += f"Best solution(s):\n"
    for i, path in enumerate(best_path):
        res = re.search(r"/(((.*)-level.+)(-on-.+))/", path)
        res1 = os.path.basename(res.group(1))
        res2 = os.path.basename(res.group(2))
        res3 = os.path.basename(res.group(3))
        res4 = os.path.basename(res.group(4))
        res4 = re.sub("-on-", "_on_", res4)
        shutil.copy2(f"{RESULTS_DIR}/{res3}{res4}/result_PDBs/{res1}.pdb", f"{base_path}/solution_{best_score_pos[i]+1}_{res1}.pdb")
        shutil.copy2(f"{RESULTS_DIR}/{res3}{res4}/result_PDBs/{res1}_renum.pdb", f"{base_path}/solution_{best_score_pos[i]+1}_{res1}_renum.pdb")
        shutil.copy2(f"{RESULTS_DIR}/{res3}{res4}/result_PDBs/{res2}.pdb", f"{base_path}/solution_{best_score_pos[i]+1}_{res2}.pdb")
        shutil.copy2(f"{RESULTS_DIR}/{res3}{res4}/result_PDBs/{res2}_renum.pdb", f"{base_path}/solution_{best_score_pos[i]+1}_{res2}_renum.pdb")
        solutions_txt += f"--> {base_path}/solution_{best_score_pos[i]+1}_{res1}.pdb\n"
    
    if verbose:
        print("\n\n\033[93mVERBOSE MODE: show all intermediate exploration levels scores and best alignments for P1 vs P2 and P2 vs P1\033[0m")
        print("              \033[93mBest overall results are shown at the end.\033[0m\n")
        print("\n".join(results))
        # Saving output of simple KPAX alignment
        run_gdt(p1, p2, min_len_p1_p2, opt_prune, ori_res_num_and_chain1, ori_res_num_and_chain2, save_output=True, keep_ori_resnum=True)
        run_gdt(p2, p1, min_len_p1_p2, opt_prune, ori_res_num_and_chain2, ori_res_num_and_chain1, save_output=True, keep_ori_resnum=True)
        dest1 = os.path.join(base_path, f"{p1.name}_on_{p2.name}")
        dest2 = os.path.join(base_path, f"{p2.name}_on_{p1.name}")
        if os.path.exists(dest1):
            print(f"INFO: Overwriting {p1.name} vs. {p2.name} results")
            shutil.rmtree(dest1)
        if os.path.exists(dest2):
            print(f"INFO: Overwriting {p2.name} vs. {p1.name} results\n")
            shutil.rmtree(dest2)
        if os.path.exists(f"{RESULTS_DIR}/{p1.name}_on_{p2.name}"):
            shutil.move(f"{RESULTS_DIR}/{p1.name}_on_{p2.name}", dest1)
        if os.path.exists(f"{RESULTS_DIR}/{p2.name}_on_{p1.name}"):
            shutil.move(f"{RESULTS_DIR}/{p2.name}_on_{p1.name}", dest2)
        print("\n\n\n\n\n\t\tGLOBAL BEST\n\n"+best)
        with open(f"{base_path}/summary.txt", "a") as filin:
            filin.write("\n".join(results))
            filin.write("\n\n\n\n\n\t\tGLOBAL BEST\n\n"+best_for_file)
            filin.write(f"Results can be found here:\n--> {dest1}\n--> {dest2}\n")
        print(f"Results can be found here:\n--> {dest1}\n--> {dest2}\n")
    else:
        print("\n")
        print(best)
        with open(f"{base_path}/summary.txt", "a") as filin:
            filin.write(best_for_file)
        if os.path.exists(f"{RESULTS_DIR}/{p1.name}_on_{p2.name}"):
            shutil.rmtree(f"{RESULTS_DIR}/{p1.name}_on_{p2.name}")
        if os.path.exists(f"{RESULTS_DIR}/{p2.name}_on_{p1.name}"):
            shutil.rmtree(f"{RESULTS_DIR}/{p2.name}_on_{p1.name}")
        if os.path.exists(PDB_CLEAN_DIR):
            shutil.rmtree(PDB_CLEAN_DIR)
        if os.path.exists(PDB_STAND_DIR):
            shutil.rmtree(PDB_STAND_DIR)
    print(solutions_txt)


def get_best_scores(scores):
    """
    Retrieves best score(s) from the scores of the explore_all_graphs.

    Args:
        - scores (list of lists): Each list contains:
                                [graph.best_score, graph.pu_order, graph.best_ali, which_alignment]

    Returns:
        - best_scores (list): best tm-scores scores of best alignments
        - best_score_pos (list): positions of best scores in the total `scores` data structure
        - which_alignment (int): "1" = P1 vs P2 and "2" = P2 vs P1
    """
    best_scores = []
    maxi = -1
    best_score_pos = []
    # which alignment: "1" = P1 vs P2 and "2" = P2 vs P1
    # There may be more than one optimal solutions,
    # therefore which_alignment is a python set which can contain "1" and "2"
    which_alignment = set()
    for i, score in enumerate(scores):
        if score[0] > maxi:
            maxi = score[0]
            best_score_pos = [i]
            which_alignment = set(score[3])
        elif score[0] == maxi:
            best_score_pos.append(i)
            which_alignment.add(score[3])
    # Retrieve best scores !
    for i in best_score_pos:
        best_scores.append(scores[i][0])
    
    return best_scores, best_score_pos, which_alignment


def get_textual_scores(scores, best_score_pos, which_alignment,
                       textual_alignment_p1_vs_p2, textual_alignment_p2_vs_p1):
    """
    Output best scores textually and visually arranged.

    Args:
        - best_scores (list): best tm-scores scores of best alignments
        - best_score_pos (list): positions of best scores in the total `scores` data structure
        - which_alignment (int): "1" = P1 vs P2 and "2" = P2 vs P1
        - textual_alignment_p1_vs_p2
        - textual_alignment_p2_vs_p1

    Returns:
        - output (str): colored terminal output
        - simple_output (str): simple text output
        - best_path (str): Path of the best alignment
    """
    # Following code is printing results when program is finished
    best_path = []
    output = ""
    if len(best_score_pos) > 1:
        output += f"There are {len(best_score_pos)} optimal solutions for this alignment:\n"
    else:
        output += "There is only 1 optimal solution for this alignment:\n"
    simple_output = ""
    for i in best_score_pos:
        simple_output = output  # version without color to write in file later
        output += f"\n\n                                                 SOLUTION {i + 1}\n"
        output += f"                                                 ==========\n"
        simple_output += f"\n\n                                                 SOLUTION {i + 1}\n"
        simple_output += f"                                                 ==========\n"
        output += f"\n *  Score: \033[92m{scores[i][0]}\033[0m\n"
        simple_output += f"\n *  Score: {scores[i][0]}\n"
        output += f"{scores[i][1]}\n\n"
        simple_output += f"{scores[i][1]}\n\n"
        if "1" in which_alignment:
            simple_output += textual_alignment_p1_vs_p2
            output += textual_alignment_p1_vs_p2
        if "2" in which_alignment:
            simple_output += textual_alignment_p2_vs_p1
            output += textual_alignment_p2_vs_p1
        best_path.append((scores[i][2]))
    return output, simple_output, best_path

def extract_query_target_and_dist_from_gdt(best_gdt2_output):
    """
    Extracts the query, the target and the distances between 
    respective residues from the GDT output.

    Args:
        - best_gdt2_output (str): GDT2 output of the best solution

    Returns:
        - query (str): query sequence
        - target (str): target sequence
        - dist (str): distances between residues as a string
        - distances (list): list of distances between residues as integers
    """
    # Extract truncated query and truncated target sequences
    gdt_query = ""
    gdt_target = ""
    distances = ""
    for line in best_gdt2_output:
        chain1 = re.search(r"CHAIN_1:(.*)", line)
        chain2 = re.search(r"CHAIN_2:(.*)", line)
        dist = re.search(r"DIST\s*:(.*)", line)
        if chain1:
            gdt_query = chain1.group(1)
        if chain2:
            gdt_target = chain2.group(1)
        if dist:
            # 'x' means distance >= 10 A
            distances = [int(d) if d not in [" ", "x"] else " " for d in dist.group(1)]
    gdt_dist = "".join([str(d) for d in distances])
    return gdt_query, gdt_target, gdt_dist, distances

def get_pu_delimitations(graph, target, pu_order_ascending_target_pos, smoothed_pu_output):
    """
    According to ascending query positions, retrieve the PU delimitations.

    Args:
        - graph (Graph): Graph object
        - target (Protein): target protein
        - pu_order_ascending_target_pos (list): order of PUs according to ascending target positions
        - smoothed_pu_output (bool): if True, use the smoothed PU delimitations 

    Returns:
        - PUs_order_to_delim (dict): keys are the PU order and values are the delimitations
        - PUs_delim_to_keep (list): list of PU delimitations to keep (was it able to align to the target ?)
        - PUs_not_shown (list): list of discarded PU delimitations (query probably too long for the target)
    """
    PUs_order_to_delim = {}
    PUs_delim_to_keep = []
    PUs_not_shown = []
    for i in pu_order_ascending_target_pos:
        ali = graph.best_path[i]
        # User requests smoothed PUs in output
        if smoothed_pu_output:
            # Smoothing can result in empty positions for aligned PU positions
            if ali.all_aligned["smoothed_all_pu_aligned_positions"] == []:
                continue
            else:
                start = ali.all_aligned["smoothed_all_pu_aligned_positions"][0]
                end = ali.all_aligned["smoothed_all_pu_aligned_positions"][-1]
        else:
            start = ali.all_aligned["all_pu_aligned_positions"][0]
            end = ali.all_aligned["all_pu_aligned_positions"][-1]
        # print only PUs that can be aligned to the target
        #if end <= len(target.seq):
            PUs_delim_to_keep.append(start)
            PUs_delim_to_keep.append(end)
            # key is the PU number and value is a tuple of start & end positions
            # of a PU X. Ex:
            #                PU2
            #      start +--------+ end
            PUs_order_to_delim[i] = (start, end)
        #else:
        #    PUs_not_shown.append(i + 1)
    return PUs_order_to_delim, PUs_delim_to_keep, PUs_not_shown

def get_pu_lengths(PUs_delim_to_keep):
    """
    Calculates the lengths of the PUs
    according to their respective delimitations.

    Args:
        - PUs_delim_to_keep (list): list of PU delimitations to keep

    Returns:
        - PUs_lengths (list): list of PU lengths
    """
    pu_lengths = []
    for id in range(len(PUs_delim_to_keep)-1):
        if id % 2 == 0:
            pu_lengths.append(PUs_delim_to_keep[id + 1] - PUs_delim_to_keep[id] + 1)
    return pu_lengths

def get_pus_connections_and_new_pu_delims(line_CHAIN1, PUs_delim_to_keep, pu_lengths):
    """
    Set PUs connections symbols:  +--------++-----------+
    Set the new PU delimitations by taking into account the
    different gaps in the alignment (begining gaps and intra sequence gaps)

    Args:
        - line_CHAIN1 (str): line of the CHAIN1 sequence
        - PUs_delim_to_keep (list): list of PU delimitations to keep
        - pu_lengths (list): list of PU lengths

    Returns:
        - line_PUs_connections (str): line of the PUs connections
        - new_delims (list): list of new PU delimitations
    """
    line_PU_connection = ""
    i = 0
    # Shift the beginning of line if the alignment starts by gaps
    nb_begining_gaps = 0
    while line_CHAIN1[i] == "-":
        line_PU_connection += " "
        nb_begining_gaps += 1
        i += 1
    new_delims = {}
    # Go through the query after removing flanking gaps "-"
    # Calculate the new delimitations of PUs taking into account gaps
    idx = 0
    pu_idx = 0
    pu_len = 1
    nb_intra_gaps = 0
    for i, res in enumerate(line_CHAIN1.strip("-"), start=1):
        # Stop if we reach the maximum position of the available PUs
        # before the end of the full alignmnent
        if (i - nb_intra_gaps) > max(PUs_delim_to_keep):
            break
        # This is a gap
        if res == "-":
            line_PU_connection += "-"
            nb_intra_gaps += 1
        else:
            # We traced the full length of a PU
            if pu_len == pu_lengths[idx]:
                new_delims[PUs_delim_to_keep[pu_idx]] = i + nb_begining_gaps
                line_PU_connection += "+"
                pu_len = 1
                idx += 1
                pu_idx += 1
            # This is the beginning of a new traced PU
            elif pu_len == 1:
                new_delims[PUs_delim_to_keep[pu_idx]] = i + nb_begining_gaps
                line_PU_connection += "+"
                pu_len += 1
                pu_idx += 1
            # This is a residue of the current PU
            else:
                line_PU_connection += "-"
                pu_len += 1
    # Case when the end of the alignement is only gaps so we cannot keep
    # some of the PUs that were supposed to be there
    if len(new_delims) % 2 != 0:
       del new_delims[max(new_delims.keys())]
       PUs_delim_to_keep = [delim for delim in PUs_delim_to_keep if delim in list(new_delims.keys())]
    return line_PU_connection, new_delims, PUs_delim_to_keep

def get_smoothed_pus_connections_and_new_pu_delims(graph, line_CHAIN1, pu_order_ascending_target_pos, PUs_order_to_delim_ascending_target_pos, PUs_delim_to_keep_ascending_target_pos):
    """
    Set PUs connections symbols:  +--------++-----------+
    Set the new PU delimitations by taking into account the
    different gaps in the alignment (begining gaps and intra sequence gaps)

    Args:
        - graph (Graph): Graph object
        - line_CHAIN1 (str): line of the CHAIN1 sequence
        - pu_order_ascending_target_pos (list): list of PU order according to ascending target positions
        - PUs_delim_to_keep_ascending_target_pos (list): list of PU delimitations according to ascending target positions

    Returns:
        - line_PUs_connections (str): line of the PUs connections
        - new_delims (list): list of new PU delimitations
    """
    i = 0
    # Count the number of begining gaps
    nb_begining_gaps = 0
    while line_CHAIN1[i] == "-":
        nb_begining_gaps += 1
        i += 1

    line_CHAIN1 = line_CHAIN1.upper()
    # Remember the positions of the gaps in the line_CHAIN1
    gaps_pos_in_line_CHAIN1 = [i for i, pos in enumerate(line_CHAIN1) if pos == "-"]
    # Remove the gaps from the line_CHAIN1
    no_gaps_line_CHAIN1 = line_CHAIN1.replace("-", "")
    line_PU_connection = " " * len(no_gaps_line_CHAIN1)

    # Draw the +-------+ delimitating the PUs
    # and register the new delimitations
    new_pu_delim_pos = []
    for i in pu_order_ascending_target_pos:
        ali = graph.best_path[i]
        match = re.search(ali.all_aligned["smoothed_pu_seq"], no_gaps_line_CHAIN1)
        start = match.start()
        end = match.end()
        new_pu_delim_pos.append(start + 1)
        new_pu_delim_pos.append(end)
        line_PU_connection = line_PU_connection[:start] + "+" + "-" * len(line_PU_connection[start + 1:end - 1]) + "+" + line_PU_connection[end:]

    # Put back the gaps dynamically to the string
    i = 0
    while i < len(line_CHAIN1):
        if i in gaps_pos_in_line_CHAIN1:
            if i == 0:
                line_PU_connection = "-" + line_PU_connection
            else:
                line_PU_connection = line_PU_connection[:i] + "-" + line_PU_connection[i:]
        i += 1

    # Add beginning gaps
    line_PU_connection = " " * nb_begining_gaps + line_PU_connection[nb_begining_gaps:]

    # Associate old delimitation positions with new ones
    new_delims = {}
    cnt = 0
    for i, pos in enumerate(line_PU_connection, start=1):
        if pos == "+" and len(new_delims) < len(PUs_delim_to_keep_ascending_target_pos):
            new_delims[PUs_delim_to_keep_ascending_target_pos[cnt]] = i
            cnt += 1
    # Replace inter-PUs gaps by " " to go from "+---+-- -+-------+"" to "+---+    +-------+""
    for i, delim in enumerate(PUs_delim_to_keep_ascending_target_pos[:-1]):
        if i % 2 != 0:
            line_PU_connection = line_PU_connection[:new_delims[delim]] + \
                                 " " * (new_delims[PUs_delim_to_keep_ascending_target_pos[i + 1]] - new_delims[delim] - 1) + \
                                 line_PU_connection[new_delims[PUs_delim_to_keep_ascending_target_pos[i + 1]] - 1:]
    line_PU_connection = line_PU_connection[:new_delims[PUs_delim_to_keep_ascending_target_pos[-1]]]
    return line_PU_connection, new_delims

def draw_inter_pus_as_blank(line_PU_connection, PUs_delim_to_keep, new_delims):
    """
                                             ↓↓↓
    Set the inter PUs as blank: :  +--------+   +-----------+

    Args:
        - line_PU_connection (str): line of the PUs connections
        - PUs_delim_to_keep (list): list of PU delimitations to keep
        - new_delims (dict): keys are the original delimitations and values are the new delimitations

    Returns:
        - line_PU_connection (str): line of the PUs connections
    """
    for i in range(len(PUs_delim_to_keep) - 1):
        if PUs_delim_to_keep[i] in new_delims and PUs_delim_to_keep[i + 1] in new_delims:
            if i % 2 != 0:
                begin = new_delims[PUs_delim_to_keep[i]]
                stop = new_delims[PUs_delim_to_keep[i + 1]]
                # Check that there is something to remove between 2 PUs
                if (stop - begin) > 1:
                    nb_to_replace = stop - begin - 1
                    line_PU_connection = line_PU_connection[:begin] + \
                                            " " * nb_to_replace + \
                                            line_PU_connection[stop-1:]
    return line_PU_connection

def draw_pu_ranges(PUs_delim_to_keep, new_delims, line_PU_connection, ori_res_query):
    """
    Set the PUs ranges: ┌12    21┐┌22       34┐  ┌36                  56┐
    Positions correspond to the original positions of the PU in the query.

    Args:
        - PUs_delim_to_keep (list): list of PU delimitations to keep
        - new_delims (dict): keys are the original delimitations and values are the new delimitations
        - line_PU_positions (str): line of the PUs positions

    Returns:
        - line_PU_positions (str): line of the PUs positions
    """
    line_PU_positions = " " * len(line_PU_connection)
    for i in range(len(PUs_delim_to_keep) - 1):
        if PUs_delim_to_keep[i] in new_delims and PUs_delim_to_keep[i + 1] in new_delims:
            if i % 2 == 0:
                begin = new_delims[PUs_delim_to_keep[i]]
                stop = new_delims[PUs_delim_to_keep[i + 1]]
                delim1 = "┌" + str(ori_res_query[PUs_delim_to_keep[i]][0])
                delim2 = str(ori_res_query[PUs_delim_to_keep[i+1]][0]) + "┐"
                line_PU_positions = line_PU_positions[:begin-1] + \
                                    delim1 + \
                                    line_PU_positions[begin + len(delim1):stop - len(delim2) + 1] + \
                                    delim2 + \
                                    line_PU_positions[stop:]
    return line_PU_positions

def draw_pu_number(line_CHAIN1, pu_order_ascending_target_pos, PUs_order_to_delim, new_delims):
    """
    Draw the PUs numbers: →      PU1         PU2                 PU3
                              ┌12    21┐┌22       34┐  ┌36                  56┐

    Args:
        - line_CHAIN1 (str): line of the CHAIN1 sequence
        - pu_order_ascending_target_pos (list): list of PUs order ascending by target position
        - PUs_order_to_delim (dict): keys are the PUs order and values are the delimitations
        - new_delims (dict): keys are the original delimitations and values are the new delimitations
    
    Returns:
        - line_PU (str): line of the PUs numbers
    """
    line_PU = " " * len(line_CHAIN1.rstrip("-"))
    for pu_i in pu_order_ascending_target_pos:
        # if the PU is not out of the scope of the target sequence length
        if pu_i in PUs_order_to_delim:
            start, end = PUs_order_to_delim[pu_i]
            if start in new_delims and end in new_delims:
                new_start = new_delims[start] - 1
                new_end = new_delims[end] - 1
                # Place the string in the middle of the PU range
                index = (new_end + new_start) // 2
                line_PU = line_PU[:index - 2] + "PU" + str(pu_i + 1) + line_PU[index + 1:]
    return line_PU

def draw_matching_symbols(line_DIST, distances):
    """
    Draw the matching symbols: 
            VKSII-TLDG-GALVQVQKWD----GKSTTIK----VGVGFATRKVAGMAK
       →    ||||| .|   .|||||||:.    .||||||    .|||||||||||||.
            NTFTIPKTDYDNFLMAHLINEKDGETFQLMGLYGREPDLSSDIKERFAQLC

    Aligned distance (match <=> dist): '|' <= 1 Å
                                   ':' <= 2 Å
                                   '.' <= 3 Å

    Args:
        - line_DIST (str): line of the distances
        - distances (list): list of distances

    Returns:
        - line_DIST (str): line of the distances
    """
    line_MATCH = ""
    # Avoid beginning gaps
    for i, dist in enumerate(line_DIST):
        if dist != " ":
            if distances[i] == " ":
                line_MATCH += " "
            elif distances[i] <= 1:
                line_MATCH += "|"
            elif distances[i] <= 2:
                line_MATCH += ":"
            elif distances[i] <= 4:
                line_MATCH += "."
            else:
                line_MATCH += " "
        else:
            line_MATCH += " "
    return line_MATCH

def draw_aligned_pus_title():
    """
    Draw the aligned PUs title

    Returns:
        - The aligned PUs title
    """
    aligned_pu_title = " " * 48 + "ALIGNED PU(S)\n"
    aligned_pu_title += " " * 48 + "=============\n\n"
    return aligned_pu_title

def draw_best_alignment_title():
    """
    Draw the best alignment title

    Returns:
        - The best alignment title
    """
    alignment_title = " " * 48 + "BEST ALIGNMENT\n"
    alignment_title += " " * 48 + "==============\n\n"
    return alignment_title

def draw_single_pus(line_CHAIN1, line_MATCH, line_CHAIN2, pu_order_ascending_query_pos, PUs_order_to_delim, PUs_delim_to_keep, new_delims, ori_res_query):
    """
    Draw the single PUs details.
    One PU at a time.
    Example:
                PU 1     :CDAFVGTWKLVSSEN--F-DD---YMKE
                          ::|||||||||||||    |
                TARGET   :VEKINGEWHTIILASDKREKIEDNGNFR
                ali. pos. 12                         39
                ori. pos. 1                          22

                PU 2     :PNMIISVNGDLVTIR
                            |||||||||||||
                TARGET   :FLEQIHVLENSLVLK
                ali. pos. 41            55
                ori. pos. 23            37

    Args:
        - line_CHAIN1 (str): line of the CHAIN1 sequence
        - line_MATCH (str): line of the matching symbols
        - line_CHAIN2 (str): line of the CHAIN2 sequence
        - pu_order_ascending_query_pos (list): list of PUs order ascending by query position
        - PUs_order_to_delim (dict): keys are the PUs order and values are the delimitations
        - PUs_delim_to_keep (list): list of PU delimitations to keep
        - new_delims (dict): keys are the original delimitations and values are the new delimitations
        - ori_res_query (dict): mapping of the original query residues to the new clean query residues

    Returns:
        - single_pus (str): textual alignment of the single PUs
    """
    single_pus = ""
    for i, pu_i in enumerate(pu_order_ascending_query_pos):
        if pu_i in PUs_order_to_delim:
            start, end = PUs_order_to_delim[pu_i]
            if start in new_delims and end in new_delims:
                new_start = int(new_delims[start])
                new_end = int(new_delims[end])
                pu_section = line_CHAIN1[new_start - 1:new_end]
                matching_section = line_MATCH[new_start - 1:new_end]
                target_section = line_CHAIN2[new_start - 1:new_end]
                # Remove gaps of target's ending sequence where some PUs
                # could potentially not be aligned because they are too long for target
                target_section = target_section[:max(PUs_delim_to_keep) + 1] + target_section[max(PUs_delim_to_keep):].replace("-", "")
                delim = "ali. pos. {:<4}".format(str(new_start)) + " " * (len(target_section) - 5) + str(new_end)
                ori_delim = "ori. pos. {:<4}".format(str(ori_res_query[PUs_order_to_delim[pu_i][0]][0])) + " " * (len(target_section) - 5) + str(ori_res_query[PUs_order_to_delim[pu_i][1]][0])
                single_pus += "PU {:<2}    :{}".format(str(i + 1), pu_section.upper()) + "\n"
                single_pus += "          " + matching_section + "\n"
                single_pus += "TARGET   :" + target_section.upper() + "\n"
                single_pus += delim + "\n"
                single_pus += ori_delim + "\n\n"
    return single_pus

def hide_beginning_gaps(line_CHAIN1):
    """
    Hide the beginning gaps of the (query) sequence by replacing
    them by " ".

    Args:
        - line_CHAIN1 (str): line of the CHAIN1 sequence

    Returns:
        - new_line_CHAIN1 (str): line of the CHAIN1 sequence with the beginning gaps hidden
    """
    i = 0
    new_line_CHAIN1 = ""
    while line_CHAIN1[i] == "-":
        new_line_CHAIN1 += " "
        i += 1
    new_line_CHAIN1 = new_line_CHAIN1 + line_CHAIN1.strip("-")
    return new_line_CHAIN1

def split_ali_if_too_wide(line_PU, line_PU_connection, line_CHAIN1, line_MATCH, line_CHAIN2, line_DIST, line_PU_positions, line_POS):
    """
    Split the textual alignment if it is too wide.
    Show n = 100 chars of alignment per line.

    Args:
        - line_PU (str): line of the PU sequence
        - line_PU_connection (str): line of the PU connection
        - line_CHAIN1 (str): line of the CHAIN1 sequence
        - line_MATCH (str): line of the matching symbols
        - line_CHAIN2 (str): line of the CHAIN2 sequence
        - line_DIST (str): line of the distances
        - line_PU_positions (list): list of PUs positions
        - line_POS (str): line of the positions


    Returns:
        - splitted_text_alignment (dict): keys are the types of lines of the alignment
                                          and values are the actual lines of the alignment
    """
    n = 100
    nb_splits = math.ceil(len(line_CHAIN1) / n)
        # Keys = line type, values = corresponding chunks of size <= n
    splitted_text_alignment = {
            "PUs": [line_PU[i:i + n] for i in range(0, len(line_PU), n)],
            "ranges": [line_PU_positions[i:i + n] for i in range(0, len(line_PU_positions), n)],
            "connect": [line_PU_connection[i:i + n] for i in range(0, len(line_PU_connection), n)],
            "query": [line_CHAIN1[i:i + n] for i in range(0, len(line_CHAIN1), n)],
            "match": [line_MATCH[i:i + n] for i in range(0, len(line_MATCH), n)],
            "target": [line_CHAIN2[i:i + n] for i in range(0, len(line_CHAIN2), n)],
            "dist": [line_DIST[i:i + n] for i in range(0, len(line_DIST), n)],
            "pos": [line_POS[i:i + n] for i in range(0, len(line_POS), n)]
        }
    return splitted_text_alignment, nb_splits

def draw_ali_pos_landmark(line_CHAIN2):
    """
    Draw the alignment positions landmarks:
            ali. pos.:         10        20        30        40 ...

    Args:
        - line_CHAIN2 (str): line of the CHAIN2 sequence
    
    Returns:
        - line_ALI_POS (str): line of the ali. pos. landmarks
    """
    line_POS = ""
    iter_target = enumerate(line_CHAIN2, start=1)
    for i, _ in iter_target:
            # Don't show position of hundreds (100, 200, etc...) as it cannot be printed
        if i % 10 == 0 and str(i)[1:] != "00":
            line_POS += str(i)
            [next(iter_target) for _ in range(len(str(i)) - 1) if i < len(line_CHAIN2) - 2]
        else:
            line_POS += " "
    return line_POS

def draw_splitted_ali(PUs_not_shown, splitted_text_alignment, nb_splits):
    """
    Consider cases when just few PUs are aligned and therefore the 
    lines like "PUs" or "connect" do not spread accross several chunks

    Args:
        - PUs_not_shown (list): list of PUs not shown
        - splitted_text_alignment (dict): keys are the types of lines of the alignment
                                            and values are the actual lines of the alignment
        - nb_splits (int): number of chunks of size <= n

    Returns:
        - trunc_textual_alignment (str): The textual alignment with the partial info not shown
    """
    trunc_textual_alignment = ""
    for i in range(nb_splits):
        # Consider cases when just few PUs are aligned and therefore the
        # lines like "PUs" or "connect" do not spread accross several chunks
        if i < len(splitted_text_alignment["PUs"]):
            trunc_textual_alignment += "PUs      :" + splitted_text_alignment["PUs"][i] + "\n"
        else:
            trunc_textual_alignment += "PUs      :" + "\n"
        if i < len(splitted_text_alignment["connect"]):
            # Add indication of which PUs are not shown because for
            # example the query sequence is longer than target sequence
            if len(PUs_not_shown) != 0 and i == (nb_splits - 1):
                trunc_textual_alignment += "ori. pos.:" + splitted_text_alignment["ranges"][i] + "\n"
                trunc_textual_alignment += "connect  :" + splitted_text_alignment["connect"][i] + " (PUs not aligned: " + str(PUs_not_shown) + ")\n"
            else:
                trunc_textual_alignment += "ori. pos.:" + splitted_text_alignment["ranges"][i] + "\n"
                trunc_textual_alignment += "connect  :" + splitted_text_alignment["connect"][i] + "\n"
        else:
            if len(PUs_not_shown) != 0 and i == (nb_splits - 1):
                trunc_textual_alignment += "ori. pos.:\n"
                trunc_textual_alignment += "connect  : (PUs not aligned: " + str(PUs_not_shown) + ")\n"
            else:
                trunc_textual_alignment += "ori. pos.:\n"
                trunc_textual_alignment += "connect  :\n"
        trunc_textual_alignment += "QUERY    :" + splitted_text_alignment["query"][i].upper() + "\n"
        if i < len(splitted_text_alignment["match"]):
            trunc_textual_alignment += "match    :" + splitted_text_alignment["match"][i] + "\n"
        else:
            trunc_textual_alignment += "match    :" + "\n"
        if i < len(splitted_text_alignment["target"]):
            trunc_textual_alignment += "TARGET   :" + splitted_text_alignment["target"][i].upper() + "\n"
        else:
            trunc_textual_alignment += "TARGET   :" + "\n"
        if i < len(splitted_text_alignment["dist"]):
            trunc_textual_alignment += "dist     :" + splitted_text_alignment["dist"][i] + "\n"
        else:
            trunc_textual_alignment += "dist     :" + "\n"
        if i < len(splitted_text_alignment["pos"]):
            trunc_textual_alignment += "ali. pos.:" + splitted_text_alignment["pos"][i] + "\n"
        else:
            trunc_textual_alignment += "ali. pos.:" + "\n"
        trunc_textual_alignment += "\n"
    return trunc_textual_alignment

def draw_single_ali(PUs_not_shown, line_PU, line_PU_connection, line_CHAIN1, line_MATCH, line_CHAIN2, line_DIST, line_PU_positions, line_POS):
    """
    Draw the alignment of a single alignment

    Args:
        - PUs_not_shown (list): list of PUs not shown
        - line_PU (str): line of the PUs
        - line_PU_connection (str): line of the PUs connection
        - line_CHAIN1 (str): line of the CHAIN1 sequence
        - line_MATCH (str): line of the MATCH sequence
        - line_CHAIN2 (str): line of the CHAIN2 sequence
        - line_DIST (str): line of the DIST sequence
        - line_PU_positions (str): line of the PUs positions
        - line_POS (str): line of the ali. pos. landmarks
        
    Returns:
        - single_line_textual_alignment (str): The textual alignment of a single alignment
    """
    single_line_textual_ali = ""
    single_line_textual_ali += "PUs      :" + line_PU + "\n"
    single_line_textual_ali += "ori. pos.:" + line_PU_positions + "\n"
    if PUs_not_shown != []:
        single_line_textual_ali += "connect  :" + line_PU_connection + "  (PUs not aligned: " + str(PUs_not_shown) + ")\n"
    else:
        single_line_textual_ali += "connect  :" + line_PU_connection + "\n"
    single_line_textual_ali += "QUERY    :" + line_CHAIN1.upper() + "\n"
    single_line_textual_ali += "match    :" + line_MATCH + "\n"
    single_line_textual_ali += "TARGET   :" + line_CHAIN2.upper() + "\n"
    single_line_textual_ali += "dist     :" + line_DIST + "\n"
    single_line_textual_ali += "ali. pos.:" + line_POS + "\n"
    single_line_textual_ali += "\n\n"
    return single_line_textual_ali

def draw_textual_alignment(query, target, graph, ori_res_query, smoothed_pu_output):
    """
    This functions prints the best alignement of a query to a target from
    a given Graph object.
    It represents full query sequence, full target sequence and full PUs sequences.
    This means that all the gaps will be represented in the alignment.

    Args:
        query (str): The "query" as Protein object
        target (str): The "target" as Protein object
        graph (Graph): a graph.Graph object
        ori_res_query (dict): the original residue number of the query
        ori_res_target (dict): the original residue number of the target

    Returns:
        The best textual alignment
    """
    line_PU = ""
    line_PU_connection = ""
    line_CHAIN1 = ""
    line_MATCH = ""
    line_CHAIN2 = ""
    line_DIST = ""
    textual_alignment = ""

    # Set full target sequence
    target.set_1d_seq()
    query.set_1d_seq()

    # Get the order of PUs aligned against target according to ascending target positions
    pu_order_ascending_target_pos = graph.get_pu_order_from_target_ascending_pos(graph.best_path)
    # Get the order of PUs aligned against target according to ascending query positions
    pu_order_ascending_query_pos = graph._get_pu_order_from_query_ascending_positions()
    # Extract the alignment resulting from the gdt2.pl output of the best solution (sequences and distances)
    line_CHAIN1, line_CHAIN2, line_DIST, distances = extract_query_target_and_dist_from_gdt(graph.best_gdt2_output)
    # Get PU delimitations according to ascending target positions
    PUs_order_to_delim_ascending_target_pos, PUs_delim_to_keep_ascending_target_pos, PUs_not_shown = get_pu_delimitations(graph, target, pu_order_ascending_target_pos, smoothed_pu_output)
    # Get the length of the PUs in the right order
    pu_lengths = get_pu_lengths(PUs_delim_to_keep_ascending_target_pos)
    # In some cases, the smallest PU is already bigger than the target,
    # so we cannot use any of them.
    if PUs_delim_to_keep_ascending_target_pos != []:
        # Set PUs connections symbols:  +--------++-----------+
        if smoothed_pu_output:
            line_PU_connection, new_delims = get_smoothed_pus_connections_and_new_pu_delims(graph, line_CHAIN1, pu_order_ascending_target_pos, PUs_order_to_delim_ascending_target_pos, PUs_delim_to_keep_ascending_target_pos)
            #PUs_delim_to_keep_ascending_target_pos = list(new_delims.keys())
            #pu_order_ascending_target_pos = list(new_delims.values())
        else:
            line_PU_connection, new_delims, PUs_delim_to_keep_ascending_target_pos = get_pus_connections_and_new_pu_delims(line_CHAIN1, PUs_delim_to_keep_ascending_target_pos, pu_lengths)
            # Replace inter-PUs connections by blanks: +--------+    +-----------+
            line_PU_connection = draw_inter_pus_as_blank(line_PU_connection, PUs_delim_to_keep_ascending_target_pos, new_delims)
        # If some PU delimitations could not align, we need to remove the corresponding PUs from the alignment
        nb_pus_left = int((len(PUs_delim_to_keep_ascending_target_pos)/2))
        pu_order_ascending_target_pos = pu_order_ascending_target_pos[:nb_pus_left]
        pu_order_ascending_query_pos = pu_order_ascending_query_pos[:nb_pus_left]
        # Set PUs position range: ┌12    21┐  ┌22       34┐
        #                         +--------+  +-----------+
        line_PU_positions = draw_pu_ranges(PUs_delim_to_keep_ascending_target_pos, new_delims, line_PU_connection, ori_res_query)                    
        # Set order of PUs number for "PUs" text line:       PU1           PU2
        #                                                 ┌12    21┐  ┌22       34┐
        #                                                 +--------+  +-----------+
        line_PU = draw_pu_number(line_CHAIN1, pu_order_ascending_target_pos, PUs_order_to_delim_ascending_target_pos, new_delims)
        # Set the "matching symbols" line
        line_MATCH = draw_matching_symbols(line_DIST, distances)
        # Textual output for single PU / line matched to the
        # corresponding region on the target
        textual_alignment += draw_aligned_pus_title()
        # Draw each PU details in the textual alignment
        # One pU at a time
        textual_alignment += draw_single_pus(line_CHAIN1, line_MATCH, line_CHAIN2, pu_order_ascending_query_pos, PUs_order_to_delim_ascending_target_pos, PUs_delim_to_keep_ascending_target_pos, new_delims, ori_res_query)
        # Replace beginning gaps of the alignment by blanks
        line_CHAIN1 = hide_beginning_gaps(line_CHAIN1)
        # Remove gaps of query's and target's ending sequences because some PUs of the query
        # could potentially not be aligned because they are too long for target
        #if max(PUs_delim_to_keep_ascending_target_pos) in new_delims:
        #    line_CHAIN2 = line_CHAIN2[:new_delims[max(PUs_delim_to_keep_ascending_target_pos)]] + line_CHAIN2[new_delims[max(PUs_delim_to_keep_ascending_target_pos)]:].replace("-", "")
        #    line_CHAIN1 = line_CHAIN1[:new_delims[max(PUs_delim_to_keep_ascending_target_pos)]] + line_CHAIN1[new_delims[max(PUs_delim_to_keep_ascending_target_pos)]:].replace("-", "")
        # Line of positions for the textual output
        # Set a position every 10 characters
        line_POS = draw_ali_pos_landmark(line_CHAIN2)

        # Some PUs may not have been aligned, so we shorten the lines 
        # to remove the corresponding positions (MATCH and DIST) 
        line_MATCH = line_MATCH[:line_PU_connection.rindex("+") + 1]
        line_DIST = line_DIST[:line_PU_connection.rindex("+") + 1]

        n = 100
        # Split the textual alignment to print a maximum of 100 caracters / line
        if len(line_CHAIN1) > n:
            splitted_text_alignment, nb_splits = split_ali_if_too_wide(line_PU, line_PU_connection, line_CHAIN1, line_MATCH, line_CHAIN2, line_DIST, line_PU_positions, line_POS)
            textual_alignment += draw_best_alignment_title()
            textual_alignment += draw_splitted_ali(PUs_not_shown, splitted_text_alignment, nb_splits)
        else:  # If the alignment is small enough, we can print it in one line
            #textual_alignment = ""
            textual_alignment += draw_best_alignment_title()
            textual_alignment += draw_single_ali(PUs_not_shown, line_PU, line_PU_connection, line_CHAIN1, line_MATCH, line_CHAIN2, line_DIST, line_PU_positions, line_POS)
        textual_alignment += "Aligned distance (match <=> dist): '|' <= 1 Å\n"+35*" "+"':' <= 2 Å\n"+35*" "+"'.' <= 3 Å\n"
    # Some PUs are not shown because they are too long for the target
    if PUs_not_shown != []:
        textual_alignment += "\n\nINFO: These PUs could not be aligned because they are too long for the target: " + str(PUs_not_shown) + "\n"
    return textual_alignment

def set_exploration_limits(prot, nb_pus_requested, exploration_level):
    """
    This function sets the limits of the exploration of the graph.

    Args:
        prot (Protein): the protein for which we set exploration limits
        nb_pus_requested (int): the number of PUs requested
        exploration_level (int): the level of exploration given by the user

    Returns:
        - exploration_level (int): the level of exploration determined after checking the limits
        - max_seg_level (int): the maximum level of segmentation for Peeling
    """
    max_seg_level = None
    if len(prot.PUs_per_level) == 0:
        print(f"\n* {prot.name}")
        print(f"   - Protein Peeling program could not find any PU in {prot.name}. ICARUS will only perform a simple KPAX when {prot.name} is used as the query")
        max_seg_level = 0
        exploration_level = 0
    else:
        for level, PUs in prot.PUs_per_level.items():
            level += 1
            if len(PUs) in nb_pus_requested:
                # Stop when segmentation level gives the number of PUs
                # requested by the user (if it can it will set the highest possible)
                max_seg_level = level
            elif len(PUs) > nb_pus_requested[-1] and max_seg_level is None:
                max_seg_level = level - 1
                new_exploration_level = exploration_level - 1
                print(f"\n* {prot.name}")
                print(f"   - The protein {prot.name} could not be segmented in {nb_pus_requested[-1]} PUs.")
                print(f"     However, it seems that there are other exploration levels available.")
                print(f"     --> For this job, we are setting the new exploration level to the level just below the one requested: {exploration_level} -> {new_exploration_level}")
                exploration_level = new_exploration_level
                break
            # The number of PUs requested is higher than the number of PUs available
            # So we set to max level possible for according max number of PUs
            if level == len(prot.PUs_per_level) and \
                    max_seg_level is None:
                max_nb_pus = max(len(PUs) for _, PUs in prot.PUs_per_level.items())
                new_exploration_level = NB_PUS_2_EXPLORE_LEVEL[max_nb_pus]
                print(f"\n* {prot.name}")
                print(f"   - Maximum number of PUs: {max_nb_pus}")
                print(f"     --> Setting new exploration level to maximum possible: {exploration_level} -> {new_exploration_level}")
                # The max number of PUs is necessarily in the last level of segmentation
                max_seg_level = len(prot.PUs_per_level)
                exploration_level = new_exploration_level
    return exploration_level, max_seg_level
    
def large_job_warning(max_len_p1_p2, exploration_level_p1, exploration_level_p2, force):
    """
    Print a warning if the job is large.
    Inform about CPU and memory usage.

    Args:
        max_len_p1_p2 (int): maximum length of the P1 and P2 sequences
        exploration_level_p1 (int): exploration level for P1
        exploration_level_p2 (int): exploration level for P2
        force (bool): if True, bypass the warning
    """
    if max_len_p1_p2 > 200 and \
        (exploration_level_p1 is not None and exploration_level_p1 >= 4 or
        exploration_level_p2 is not None and exploration_level_p2 >= 4 or
        exploration_level >= 4):
            print('\n\033[1;93mWarning\033[0m: jobs with proteins >200 residues '
                'and exploration level >= 4\n\tusually require a high amount '
                'of computational resources (RAM and disk memory).')
            if not force:
                answer = ""
                while answer not in ["y", "n"]:
                    answer = input('\nPlease acknowledge the warning.'
                                '\nUse option --force to bypass this step.'
                                '\nProceed to calculations anyway ? [y/n]').lower()
                if answer == "n":
                    sys.exit()
            else:
                print("\nUser is bypassing the warning")

def clean_input_pdb_files(path1, path2, chain1, chain2):
    """
    Clean the input PDB files.

    Args:
        path1 (str): path to the PDB file of P1
        path2 (str): path to the PDB file of P2
        chain1 (str): pdb chain to parse in the query
        chain2 (str): pdb chain to parse in the target
    
    Returns:
        new_path1 (str): path to the cleaned PDB file of P1
        new_path2 (str): path to the cleaned PDB file of P2
        ori_res_num_and_chain1 (dict): original residue numbers and chain(s) of P1
        ori_res_num_and_chain2 (dict): original residue numbers and chain(s) of P2
    """
    print("\nClean input PDB files ... ", end="")
    # Clean input PDB files:
    # - Reindex
    # - Remove unsolved residues (ones which have no CA atom)
    # - Check if it is a legitimate amino acids PDB
    ori_res_num_and_chain1 = {}
    ori_res_num_and_chain2 = {}
    new_path1 = os.path.join(TMP_DIR, os.path.basename(os.path.splitext(path1)[0]))
    new_path2 = os.path.join(TMP_DIR, os.path.basename(os.path.splitext(path2)[0]))
    len1 = utils.reformat_struct(path1, ori_res_num_and_chain1, chain1, new_path = new_path1)
    len2 = utils.reformat_struct(path2, ori_res_num_and_chain2, chain2, new_path = new_path2)
    print("done\n")
    return new_path1, new_path2, ori_res_num_and_chain1, ori_res_num_and_chain2, len1, len2


def parse_arguments():
    """
    This functions parses the arguments given to the script.
    It returns the arguments as an object.

    Args:
        None

    Returns:
        args: The arguments as an object
    """
    def check_path(path):
        """
        Check the input path sanity
        """
        if os.path.isfile(path):
            return path
        raise argparse.ArgumentTypeError(f"Error: file:{path} does not exist")

    def check_cpu(nb_cpu):
        """
        Check if the user input CPU nb is valid
        """
        try:
            nb_cpu = int(nb_cpu)
        except ValueError as e:
            print("Unable to cast ", type(nb_cpu), " into integer (int)")
            raise argparse.ArgumentTypeError(
                "Error option -c/--cpu: please input an integer or string integer"
            ) from e
        if 0 <= nb_cpu <= cpu_count():
            return nb_cpu
        raise argparse.ArgumentTypeError(
            f"Error option -c/--cpu: nb_cpu should be 0 <= nb_cpu <= {cpu_count()}")

    def check_chain(chain):
        """
        Check if the user input chain is valid
        """
        if (chain.isalpha() or chain.isdigit()) and len(chain) == 1:
            return chain
        raise argparse.ArgumentTypeError(
            "Error option -a/--chain: PDB chain should be a single alphabetic or digit character [a-zA-Z0-9]")

    def check_pu_size(size):
        """
        Check if the user input PU size is valid
        """
        try:
            size = int(size)
        except ValueError as e:
            print("Unable to cast ", type(size), " into integer (int)")
            raise argparse.ArgumentTypeError(
                "Error option -m: please input an integer or string integer"
            ) from e
        if 15 <= size <= 99:
            return size
        raise argparse.ArgumentTypeError(
            "Error option -m: PU size should be 15 <= min <= 99")

    def check_prune(score):
        """
        Check if the user input pruning TM-score threshold is valid
        """
        try:
            score = float(score)
        except ValueError as e:
            print("Unable to cast ", type(score), " into float")
            raise argparse.ArgumentTypeError(
                "Error option -p: please input a float or string float"
            ) from e
        if 0.0 <= score <= 1.0:
            return score
        raise argparse.ArgumentTypeError(
            "Error option -p: score should be 0.0 <= min <= 1.0")

    def check_exploration_level(lvl):
        """
        Check if the user input refinement level is valid.
        It should be 1 <= lvl <= 5
        """
        try:
            lvl = int(lvl)
        except ValueError as e:
            raise argparse.ArgumentTypeError(
                "Error option -l: please input an integer") from e
        if 1 <= lvl <= 5:
            return lvl
        raise argparse.ArgumentTypeError(
            "Error option -l: refinement level should be 1 <= lvl <= 5")

    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''\
                    Icarus is a flexible structural alignment program.
                    It takes as input 2 pdb files and returns the optimal
                    alignment based on different protein exploration levels.'''),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
                Explanations of ICARUS output.

                ICARUS generates results in a directory named 'icarus_output'
                in the directory from which you launch ICARUS.

                Non-verbose output
                ------------------
                ./icarus_output/results/query_and_target:
                 ├── solution1_query-on-target-level_X_N_PUs.pdb
                 └── summary.txt (terminal textual output)

                Verbose output
                --------------
                See the README for detailed information on verbose output'''))
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-p1",
                        "--protein1",
                        help="Path to the first protein to align",
                        type=check_path,
                        required=True)
    required.add_argument("-p2",
                        "--protein2",
                        help="Path to the second protein to align",
                        type=check_path,
                        required=True)
    optional.add_argument("-m",
                        "--min-size",
                        help=textwrap.dedent('''\
                            Minimum size of Protein Units (PUs).
                            Must be 15 <= min <= 99, default 15'''),
                        default=15,
                        type=check_pu_size)
    optional.add_argument("-c1",
                        "--chain1",
                        help=textwrap.dedent('''\
                            PDB Chain of the query.
                            Single alphabetical character [A-B].
                            Default is chain A.'''),
                        default="A",
                        type=check_chain)
    optional.add_argument("-c2",
                        "--chain2",
                        help=textwrap.dedent('''\
                            PDB Chain of the target.
                            Single alphabetical character [A-B].
                            Default is chain A.'''),
                        default="A",
                        type=check_chain)
    optional.add_argument(
        "-l",
        "--exploration-level",
        help=textwrap.dedent('''\
                              The exploration level determines up to how many PUs to
                              consider to build the graph of solutions.
                              A high exploration level will consider more PUs,
                              therefore explore more possibilities and
                              potentially find better results, but it will
                              also increase exponentially the complexity of
                              calculations and runtime.
                              A low exploration level can potentially miss good results.
                              A good trade-off is to set level 2 (default) or 3.
                              Attention !! Levels 4 and 5 may require a high amount 
                              of memory and will run longer.
                              Available levels and corresponding number of PUs (up to):
                                1 -> [2, 3],
                                2 -> [4, 5],
                                3 -> 6,
                                4 -> 7,
                                5 -> 8'''),
        default=2,
        type=check_exploration_level)
    optional.add_argument(
        "-p",
        "--prune",
        help=textwrap.dedent('''\
                            The pruning threshold corresponds to a TM-score value
                            used to filter the KPAX between Protein Units and the target protein.
                            If an alignment is below this threshold, it will be pruned from the graph.
                            A high pruning threshold (TM-score value) will filter out more solutions
                            whereas a low one will prune less solutions.
                            Default value is no pruning (0.).'''),
        default=0.0,
        type=check_prune)
    group = optional.add_mutually_exclusive_group()
    optional.add_argument(
        "-t",
        "--smoothed-pu-output",
        help=textwrap.dedent('''\
                            If this option is set, the Protein Units (PUs) that were aligned to the target protein
                            are smoothed / trimmed to keep only the "core" aligned position.
                            This option only changes the final presentation of the textual alignment for "visual"
                            purposes. The output PDB files contain all the residues. Default is not set'''),
        action="store_true",
        default=False)
    optional.add_argument(
        "-e",
        "--sequential",
        help=textwrap.dedent('''\
                            If this option is set, the program considers only solutions with consecutive Protein Units.
                            Either ascending or descending order.'''),
        action="store_true",
        default=False)
    optional.add_argument(
        "-f",
        "--force",
        help="Bypass asking user confirmation for exploration level >= 4",
        action="store_true",
        default=False)
    optional.add_argument(
        "-c",
        "--cpu",
        help=f"How many CPUs to use. Default all (0). Max on this computer is: {cpu_count()}",
        default=0,
        type=check_cpu)
    optional.add_argument(
        "-v",
        "--verbose",
        help=textwrap.dedent('''\
                        Set verbose mode: print longer output and generate
                                          intermediate results and alignments'''),
        action="store_true",
        default=False)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    # Parse CLI arguments
    args = parse_arguments()

    # Set variables
    p.Protein.seg_size = args.min_size
    exploration_level = args.exploration_level
    opt_prune = args.prune
    verbose = args.verbose
    force = args.force
    nb_cpu = args.cpu
    chain1 = args.chain1
    chain2 = args.chain2
    smoothed_pu_output = args.smoothed_pu_output
    sequential = args.sequential
    # Default is 0 -> use all cores
    if nb_cpu == 0:
        nb_cpu = cpu_count()

    path1 = args.protein1
    path2 = args.protein2
    # Check, reformat and clean input PDB files
    new_path1, new_path2, ori_res_num_and_chain1, ori_res_num_and_chain2, len1, len2 = clean_input_pdb_files(path1, path2, chain1, chain2)
    nb_pus_requested = INTERVALS[exploration_level]
    min_len_pu = len1 // max(nb_pus_requested)
    min_len_pu1 = min_len_pu if 15 <= min_len_pu <= 25 else 15
    p1 = p.Protein(new_path1, ori_path=path1, min_len_pu=min_len_pu1)
    p1.set_1d_seq()
    print(f"\n    {p1.length} aa\n    Seq: {p1.seq}")
    min_len_pu = len2 // max(nb_pus_requested)
    min_len_pu2 = min_len_pu if 15 <= min_len_pu <= 25 else 15
    p2 = p.Protein(new_path2, ori_path=path2, min_len_pu=min_len_pu2)
    p2.set_1d_seq()
    print(f"\n    {p2.length} aa\n    Seq: {p2.seq}")
    min_len_p1_p2 = min(p1.length, p2.length)
    max_len_p1_p2 = max(p1.length, p2.length)
    nb_pus_requested = INTERVALS[exploration_level]
    exploration_level_p1 = None
    exploration_level_p2 = None

    base_path = os.path.join(RESULTS_DIR, p1.name + "_and_" + p2.name)
    if os.path.exists(base_path):
        print(f"\n\nWARNING: Overwriting existing results at {base_path}\n")
        shutil.rmtree(base_path)

    # Set limits for Protein 1
    exploration_level_p1, g.GraphPU.max_seg_level_p1 = set_exploration_limits(p1, nb_pus_requested, args.exploration_level)

    # Set limits for Protein 2
    exploration_level_p2, g.GraphPU.max_seg_level_p2 = set_exploration_limits(p2, nb_pus_requested, args.exploration_level)

    large_job_warning(max_len_p1_p2, exploration_level_p1, exploration_level_p2, force)

    main(p1, p2, min_len_p1_p2, exploration_level,
        exploration_level_p1, exploration_level_p2,
        opt_prune, smoothed_pu_output,
        ori_res_num_and_chain1, ori_res_num_and_chain2, sequential, nb_cpu, verbose)

    # Add program runtime to terminal output and summary.txt
    runtime = time.time() - start_time
    print("Total runtime: {:.1f} seconds".format(runtime))
    with open(os.path.join(base_path, "summary.txt"), "a") as f_out:
        f_out.write("Total runtime: {:.1f} seconds".format(runtime))
    utils.clean()
