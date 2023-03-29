"""
This module implements the GraphPU class.
"""

import multiprocessing
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
from collections import OrderedDict
from functools import partial
import itertools

import networkx as nx
from networkx.classes.function import nodes
import numpy as np
from networkx import algorithms

import src.utils as utils

from .alignment import Alignment
from .protein import Protein

# When icarus.py is executed, this equals: /path/to/icarus.py
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK_DIR = os.path.join(os.getcwd(), "icarus_output")
RESULTS_DIR = os.path.join(WORK_DIR, "results")
GDT = os.path.join(PROJECT_DIR, "bin", "gdt2.pl")
TMP_DIR = utils.TMP_DIR

class GraphPU:
    """
    Rooted Directed Graph: nodes represent an alignment between a
    Protein Unit (PU) edges between two nodes represent the alignment
    """

    max_seg_level_p1 = None
    max_seg_level_p2 = None

    def __init__(self,
                 query,
                 target,
                 seg_level,
                 expl_level,
                 min_len_p1_p2,
                 ori_res_num_and_chain_query,
                 ori_res_num_and_chain_target,
                 nb_cpu,
                 opt_prune,
                 seed_alignment,
                 smoothed_pu_output,
                 sequential):
        """
        Initialize a DiGraph which first's and only node is the initial state
        of the alignement ie naive all query PUs vs target alignments.

        Args:
            - query, Protein:
                Query protein to align to target.
            - target, Protein:
                Target protein on which PUs will be aligned.
            - seg_level, int:
                Segmentation level at which the alignment should be realised.
            - expl_level, int:
                Exploration level at which the alignment should be realised.
            - ori_res_num_and_chain_query (dict): Mapping of the original residue number and chain(s) of the query to the new numerotation and chain.
            - ori_res_num_and_chain_target (dict): Mapping of the original residue number and chain(s) of the target to the new numerotation and chain.
            - nb_cpu (int): Number of CPUs to use for multiprocessing
            - opt_prune (float): TM-score threshold to prune TM-align alignments of PUs
            - seed_alignment (bool): If True, the TM-align alignments are seeded with a
                                        fixed alignment (used when query and target are the same)
            - smoothed_pu_output (bool): If True, the PUs are trimmed/smoothed for the final textual alignment only
            - sequential (bool): If True, do a sequential alignment: keep paths with consecutive PUs only
        """

        __slots__ = ("query", "target", "seg_level", "expl_level", 
                    "min_len_p1_p2", "nb_cpu", "opt_prune", "seed_alignment"
                    "best_ali", "best_score", "final_alis", "best_path", 
                    "best_gdt2_output", "best_pu_order", "pu_order_text",
                    "pu_range", "all_alis", "query", "target", "seg_level",
                    "expl_level", "nb_PUs", "opt_prune")

        self.best_ali = None
        self.best_score = None
        self.final_alis = None  # Instantiated when self.merge_all() is called.
        self.best_path = None
        self.best_gdt2_output = None
        self.best_pu_order = None
        self.pu_order_text = None
        self.pu_range = {}
        self.all_alis = None
        self.query = query
        self.target = target
        self.seg_level = seg_level
        self.expl_level = expl_level
        self.nb_PUs = len(query.PUs_per_level[self.seg_level - 1])
        self.ori_res_num_and_chain_query = ori_res_num_and_chain_query
        self.ori_res_num_and_chain_target = ori_res_num_and_chain_target
        self.best_ori_query_renum_pdb = None
        self.opt_prune = opt_prune
        self.seed_alignment = seed_alignment
        self.smoothed_pu_output = smoothed_pu_output
        self.sequential = sequential
        self.succeeded = False
        # Graph representing all possible alignments.
        # Nodes are state where different alignments are possible on a
        # target. Alignments are stored in 'alis' parameter of the node.
        # Targets are updated after each alignment and are intrinsic to
        # alignments
        # Edges represents performed alignments, stored in attribute ali
        self.graph = nx.DiGraph()
        # First set of alignments
        self.alis = Alignment.multiple_alignment(query, target, opt_prune, seed_alignment, seg_level)
        # Number of nodes to create
        # key=nb of PUs values=nb of nodes to create accordingly
        self.edges_todo_all_d = {
            1: 1,
            2: 4,
            3: 15,
            4: 64,
            5: 325,
            6: 1956,
            7: 13699,
            8: 109600,
            9: 986409
        }
        self.build_graph(nb_cpu)
        self.merge_all(nb_cpu)
        succeeded = self.compute_scores(min_len_p1_p2, nb_cpu, opt_prune, seed_alignment)
        if succeeded:
            # Smooth the positions of the PUs in the best alignments to be able to print them correctly
            # Remove PU positions that are isolated:
            # --> if less than 5 contiguous positions or separated by 2+ positions (default values)
            if smoothed_pu_output:
                for ali in self.best_path:
                    ali.smooth_positions(min_contiguous=5, n_gaps=2)
            pu_order = self._get_pu_order_from_query_ascending_positions()
            self.write_output_PUs_query_regions(pu_order)

            # Prepare output PDB result files
            # 1   - PDB file: Best aligned query structure + target (renumbered)
            # 1.1 - PDB file: Best aligned query structure + target (original residue numbers)
            # 2   - PDB file: Best aligned query structure (renumbered)
            # 2.1 - PDB file: Best aligned query structure (original residue numbers)
            base_path = os.path.join(RESULTS_DIR, query.name + "_on_" + target.name)
            os.makedirs(f"{base_path}/result_PDBs", exist_ok=True)
            solution_pdb_renum = f"{base_path}/result_PDBs/{query.name}-level_{expl_level}_{self.nb_PUs}_PUs_renum.pdb"
            solution_pdb = f"{base_path}/result_PDBs/{query.name}-level_{expl_level}_{self.nb_PUs}_PUs.pdb"
            solution_and_target_pdb_renum = f"{base_path}/result_PDBs/{query.name}-level_{expl_level}_{self.nb_PUs}_PUs-on-{target.name}_renum.pdb"
            solution_and_target_pdb = f"{base_path}/result_PDBs/{query.name}-level_{expl_level}_{self.nb_PUs}_PUs-on-{target.name}.pdb"
            # Move the best aligned query structure to the result directory
            shutil.copy2(self.best_ori_query_renum_pdb, solution_pdb_renum)
            shutil.copy2(self.best_ali, solution_pdb)
            target_renum_path = os.path.join(TMP_DIR, target.name)
            # Copy the original target and reformat it
            shutil.copy2(target.ori_path, target_renum_path)
            Protein.reformat_struct(target_renum_path)
            # Concatenate the initial target to the pdb file containing the
            # best aligned query structure, to be able to view both structures in pymol
            with open(solution_pdb, "r") as f, open(solution_pdb_renum, "r") as f1, open(target.ori_path, "r") as f2, \
                 open(target_renum_path, "r") as f3, open(solution_and_target_pdb_renum, "w") as f_out, \
                 open(solution_and_target_pdb, "w") as f_out1 :
                # Write the QUERY
                f_out.write(f"HEADER    query_{query.name}-level_{expl_level}_{self.nb_PUs}_PUs\n")
                f_out1.write(f"HEADER    query_{query.name}-level_{expl_level}_{self.nb_PUs}_PUs\n")
                for line in f:
                    f_out1.write(line)
                for line in f1:
                    f_out.write(line)
                f_out.write("END\n")
                f_out1.write("END\n")
                # Write the TARGET
                f_out.write(f"HEADER    target_{target.name}\n")
                f_out1.write(f"HEADER    target_{target.name}\n")
                for line in f2:
                    f_out1.write(line)
                for line in f3:
                    f_out.write(line)
                f_out.write("END\n")
                f_out1.write("END\n")
            # Build directory for all intermediate results (when --verbose)
            # This is destroyed at the end when not --verbose
            os.makedirs(f"{base_path}/intermediate/alignments_level_{expl_level}_{self.nb_PUs}_PUs", exist_ok=True)
            for ali in self.best_path:
                dest = f"{base_path}/intermediate/alignments_level_{expl_level}_{self.nb_PUs}_PUs/ali_{ali}"
                if os.path.exists(dest):
                    shutil.rmtree(dest)
                shutil.copytree(f"{ali.path}", dest)
            utils.clean()
            self.succeeded = True


    def _get_graph_skeletton(self, graph_skeletton, node="T"):
        """
        This method builds recursively the skeletton (pseudo adjacency matrix)
        of the PUs graph. It is not doing the actual alignments, only generating
        a dictionnay which represents the maximum depth and size graph that can
        be generated.

        Args:
            - graph_skeletton (dict of dicts): dictionnary representing the
                                                graph, empty at the beginning,
                                                it fills up recursively during
                                                runtime
            - node (str): Name of the node currently expanded

        Returns:
            - skeletton graph (dict of dicts)
            Example:
            graph_skeletton(graph_skeletton, 1, node='T')
                level |  node  |   successors
            >>> {0:     {'T':    ['T1', 'T2', 'T3']},
                 1:     {'T1':   ['T11', 'T12', 'T13'],
                         'T2':   ['T21', 'T22', 'T23'],
                         'T3':   ['T31', 'T32', 'T33']}
                 2: ...}
        """
        # Stop condition for recursive function: peeling depth limit
        level = len(node) - 1  # T[0..n]
        if level > len(self.alis):
            return True
        # Range goes from 1-->nb_alis, the deeper we go in the graph the less
        # PUs we need to align and so the less edges we need to create
        for id_pu in self.alis:
            # Do not try to realign same PU again
            if str(id_pu) in list(node):
                continue
            new_node = f"{node}{id_pu}"
            # Keys = nodes ; values = corresponding successors in DiGraph
            try:
                graph_skeletton[level][node][id_pu] = new_node
            except KeyError:
                graph_skeletton[level][node] = {}
                graph_skeletton[level][node][id_pu] = new_node
            self._get_graph_skeletton(graph_skeletton, f"{node}{id_pu}")
        return True

    def _add_nodes(graph_edges, graph_nodes, node, successors, alis, encoded_pus_at_curr_level, tree_pu_vectors, init_target_len):
        """
        DO NOT PUT "self" TO THIS FUNCTION,
        IT SLOWS DOWN MULTIPROCESSING BY A LOT !

        For a given node N, and its associated alignments A1..Ai,
        this function builds nodes N(A1)..N(Ai) and the edges between
        N and N(A1)..N(Ai).
        The attribute of a node is the alignments left to do against the target.
        The attribute of an edge is the alignment between connected nodes.
        The nodes and edges are stored inside pythonic data structures, shared
        in the memory during multiprocessing and updated during runtime.

        Args:
            - graph_edges (list): Stores all edges as a list of 3-tuples:
                            [(u, v, attr), ...] <==> [('T1', 'T2', 'ali'=ali)]
            - graph_nodes (dict): Keep order of nodes by smallest ('T': root)
                            to longest (T112: deepest).
                            We need this order to expand the graph correctly
            - node (str): Name of the node (Ex: "T11")
            - successors (list): List of out-degree nodes connected to this node
            - alis (list[Alignment,]): List of Alignments of this node

        Returns:
            None
        """
        for id_pu, ali in alis.items():
            if not ali.success:  # if alignment has failed, don't create node
                continue  # skips node
            # First node of the graph
            if node == "T":
                # Initialize vector of size equal to the full target length, to 0 (nothing aligned)
                encoded_ali = np.zeros(shape=init_target_len, dtype=np.int8)
            else:  # We are not at the root node, we need to get the encoded vector of the previous node
                encoded_ali = tree_pu_vectors[(node[:-1], int(node[-1]))]
            # Retrieve aligned core positions of the PU on the target
            aligned_pos = np.asarray([x-1 for x in ali.all_aligned["core_target_aligned_positions"]])
            # Replace positions with id_pu for aligned positions in full vector
            encoded_ali[aligned_pos] = id_pu
            tree_pu_vectors[(node, id_pu)] = encoded_ali
            # Concatenate all
            # Use strings because comparison of strings (later) is faster than comparison of numpy arrays
            encoded_ali = ''.join([str(x) for x in encoded_ali])
            encoded_pus_at_curr_level[node] = encoded_ali
            # Add node (attribute = alis to do)
            graph_nodes[successors[id_pu]] = {"alis": ali.get_new_alis(id_pu, alis)}
            # Add edge (attribute = ali done)
            graph_edges.append((node, successors[id_pu], {"ali": ali}))
            
            

    def build_graph(self, nb_cpu):
        """
        This method first sets the theoretical (max possibilities) full graph
        architecture as a pseudo adjacency matrix. Then builds the actual graph
        by multiprocessing the calculations of all the alignments for each
        segmentation level.
        The graph is stored as a Networkx DiGraph object, where nodes represent
        the alignments to do, and edges the actual alignements between PUs and
        targets.
        When a PU is aligned on a target at node i, the next alignments to be
        done on this target at nodes i+n will use a new target that omits the
        already aligned regions.

        Args:
            nb_cpu (int): Number of CPUs to use for multiprocessing

        Returns:
            - None
        """
        # Build the skeletton of the graph:
        # --> precalculation of all the alignments that need to be done.
        # For each level of segmentation from 0 to nb of alis after 1st SWORD:
        #    keys=nodes and values=list of successor nodes.
        # Order of nodes by smallest ('T': root) to longest (T1...N: deepest)
        # We need to keep this order to expand the graph correctly afterwards
        graph_skeletton = OrderedDict(
            {i: {}
             for i in range(len(self.alis) + 1)})
        self._get_graph_skeletton(graph_skeletton)

        # Stores all edges as a list of 3-tuples:
        # [(u, v, attr), ...] <==> [('T1', 'T2', 'ali'=ali)]
        graph_edges = multiprocessing.Manager().list()
        # Stores nodes and their attribute 'alis'
        # Initialized with the root node 'T' and first N alignments
        graph_nodes = multiprocessing.Manager().dict(
            {"T": {
                "alis": self.alis
            }})

        # We need to keep track of the alignments that have been done
        # This stores for each node, a vector of zeros, and the PU id
        # at positions where it was aligned.
        # ex: { T1234: [0,0,0,1,1,1,1,1,1,1,0,0,0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,0,0,0,4,4,4,4,4],
        #       ...
        #      }
        tree_pu_vectors = multiprocessing.Manager().dict()

        # Each alignment in each node of the graph is represented
        # as a string (length of the target) of zeros and a unique identifier of each PU (from 1 to ~8-9)
        # to identify core aligned positions of the PU on the target
        # example: "000000001111111111100000000000000000000"
        # key=node, value=encoded aligned positions
        encoded_pus_at_curr_level = multiprocessing.Manager().dict()
        # use the length of the target sequence as max length to store encoded aligned positions of PUs
        init_target_len = graph_nodes["T"]["alis"][1].target.length

        removed_nodes = set()
        removed_nodes_per_level = dict()
        cnt = 0
        
        # MULTIPROCESSING
        nb_alis_done = 0
        nb_alis_todo = sum([len(graph_skeletton[level]) for level in range(len(self.alis) + 1)])
        # SWORD initially splitted the query protein into n PUs.
        # For each successive alignment of PU that has to be done onto the
        # target, this loop creates a graph where a node represents the
        # alignments that are left to do, and edges represent an alignment.
        # The graph explores all possible alignments.
        for level in range(len(self.alis) + 1):
            encoded_pus_at_curr_level = multiprocessing.Manager().dict()
            with multiprocessing.Pool(processes=nb_cpu, initializer=signal.signal, initargs=(signal.SIGINT, signal.SIG_IGN)) as p:
                try:
                    # For each node of the graph, get the alignments to do and
                    # parallelize the actual alignments calculations for successors
                    for node, successors in graph_skeletton[level].items():
                        nb_alis_done += 1
                        self.progressbar(nb_alis_done, nb_alis_todo, "Build graph")
                        # Nodes terminology represents segmentation depth
                        # Ex: T1, T2, T3 => level 1 ;
                        # T11, T12, T21, T22, T31, T32 => level 2 ...
                        # skip non-existant nodes (failed alignment probably)
                        try:
                            alis = graph_nodes[node]["alis"]
                        except KeyError:
                            continue
                        args = (graph_edges, graph_nodes, node, successors, alis, encoded_pus_at_curr_level, tree_pu_vectors, init_target_len)
                        p.apply_async(GraphPU._add_nodes, args)
                    p.close()
                    p.join()
                except KeyboardInterrupt:
                    if os.path.exists(TMP_DIR):
                        shutil.rmtree(TMP_DIR, ignore_errors=True)
                        print("\nQuitting gracefully, bye !")
                    sys.exit(0)
            # BRANCH & BOUND
            # For each node at a given level, check if the alignments of successors
            # can be skipped because the order of alignment was already done once.
            if 1 < level < len(self.alis):
                removed_nodes_per_level[level] = []
                # Generate combinations of possible PUs at this level
                pu_combinations = list(itertools.combinations("".join([str(x) for x in list(range(1, len(self.alis)))]), len(list(graph_skeletton[level])[0][1:])))
                pu_combinations = [''.join(comb) for comb in pu_combinations]
                for combination in pu_combinations:
                    # Generate all permutations of given PUs
                    perms = ["".join(sequence) for sequence in itertools.permutations(combination)]
                    len_perms = len(perms)
                    # Explore only necessary permutations (upper right non diagonal triangle matrix) 
                    for i in range(len_perms - 1):
                        done_node = "T" + perms[i]
                        try:
                            done_encoded_ali = encoded_pus_at_curr_level[done_node]
                        except KeyError:
                            continue
                        for j in range(i + 1, len_perms):
                            current_node = "T" + perms[j]
                            try:
                                current_encoded_ali = encoded_pus_at_curr_level[current_node]
                            except KeyError:
                                continue
                            # This condition is necessary because the next loop deletes some nodes of the current level
                            if current_node in graph_skeletton[level] and done_node in graph_skeletton[level]:
                                # Check if the alignment is already being explored in another branch
                                if current_encoded_ali == done_encoded_ali:
                                    # Kill children of the current node
                                    for _, successor in graph_skeletton[level][current_node].items():
                                        nb_alis_done += 1
                                        removed_nodes.add(successor)
                                        removed_nodes_per_level[level].append(successor)
                                        nlevel = len(successor) - 1
                                        graph_skeletton[nlevel].pop(successor, None)
                                        graph_nodes.pop(successor, None)
                                    cnt += 1
                                    graph_skeletton[level].pop(current_node, None)
                                    graph_nodes.pop(current_node, None)

        # Delete edges of removed nodes
        graph_edges = [(i, j, k) for i, j, k in graph_edges if i not in removed_nodes]

        # Build the DiGraph from dictionaries
        self.graph.add_nodes_from(graph_nodes)
        self.graph.add_edges_from(graph_edges)
        print("")
        # Uncomment the following lines to draw the graph in PNG file.
        # A = nx.nx_agraph.to_agraph(self.graph)
        # A.layout()
        # A.draw(f"Graph.png")

    def get_terminal_nodes(self):
        """
        Returns a list of all terminal nodes in the graph.

        Args:
            None

        Returns:
            - None
        """
        return [v for v, d in self.graph.out_degree() if d == 0]

    def get_paths(self):
        """
        Returns a list of all paths constituting a full alignment:
        All paths from root node "T" to last leafs of the graph.

        Returns:
            - None
        """
        paths = []
        terminal_nodes = self.get_terminal_nodes()
        for node in terminal_nodes:
            paths.append(algorithms.shortest_path(self.graph, "T", node))
        # Pruning the graph remove some leafs so we need to remove the corresponding incomplete paths
        max_path_len = max([len(path) for path in paths])
        paths = [path for path in paths if len(path) == max_path_len]
        return paths

    def get_alis(self, path):
        """
        Returns the list of alignments that constitute the path from root node
        "T" to a leaf of the graph.

        Args:
            None

        Returns:
            - None
        """
        alis = []
        for i in range(len(path) - 1):
            alis.append(self.graph.edges[path[i], path[i + 1]]["ali"])
        return alis

    def get_all_alis(self):
        """
        Returns a list of all chains of alignments that constitutes a full
        alignment.

        Returns:
            - None
        """
        paths = self.get_paths()
        for i, path in enumerate(paths):
            paths[i] = self.get_alis(path)
        return paths

    def calc_chunksize(self, n_workers, len_iterable, factor=4):
        """
        Calculate the optimal chunksize to send to multiprocessing workers.
        Source: https://stackoverflow.com/a/54032744/6401758

        Args:
            - n_workers (int): Number of workers (a.k.a cpu_count())
            - len_iterable (int): Length of the iterable over which
                                    multiprocessing will go through

        Returns:
            - chunksize (int): Size of chunks the iterable is chopped off
                                by, and which will be submitted to the
                                process pool as separate tasks
        """
        chunksize, extra = divmod(len_iterable, n_workers * factor)
        if extra:
            chunksize += 1
        return chunksize

    def merge_all(self, nb_cpu):
        """
        Merge all paths of alignments to pdb file. Every file represents a
        solution.

        Args:
            nb_cpu (int): Number of CPUs to use for multiprocessing

        Returns:
            - None
        """
        self.all_alis = self.get_all_alis()
        len_tot = len(self.all_alis)
        # List of tuples: (path of merged alignments, list of corresponding alignments)
        paths = multiprocessing.Manager().list()

        chunksize = self.calc_chunksize(nb_cpu, len_tot)
        with multiprocessing.Pool(processes=nb_cpu, initializer=signal.signal, initargs=(signal.SIGINT, signal.SIG_IGN)) as p:
            try:
                func = partial(GraphPU.merge_alis, paths, self.query,
                            self.target, self.sequential, self.expl_level, self.nb_PUs, self.ori_res_num_and_chain_query)
                for i, _ in enumerate(
                        p.imap_unordered(func, self.all_alis, chunksize)):
                    self.progressbar(i + 1, len_tot, "Merge alignments")
                p.close()
                p.join()
            except KeyboardInterrupt:
                if os.path.exists(TMP_DIR):
                    shutil.rmtree(TMP_DIR, ignore_errors=True)
                    print("\nQuitting gracefully, bye !")
                sys.exit(0)
        # After pruning by TM-score intermediate alignments,
        # some paths may be empty. We remove them.
        for path in paths:
            if path[2] == []:
                paths.remove(path)
        print("")
        self.final_alis = paths

    @staticmethod
    def merge_alis(paths, query, target, sequential, expl_level, nb_PUs, ori_res_num_and_chain_query, alignments):
        """
        DO NOT PUT "self" TO THIS FUNCTION,
        IT SLOWS DOWN MULTIPROCESSING BY A LOT !

        Merge multiple alignments into a single pdb file. Parts merged are taken
        from the outputs of TM-align ie query protein that were subject to
        transformation.

        Args:
            - paths (Manager().list()): List of path of PUs for a whole alignment
            - query (str): Name of the query protein
            - target (str): Name of the target protein
            - sequential (bool): if True, only keep the alignements with consecutive PUs to keep the query residues order
            - expl_level (int): Exploration level
            - nb_PUs (int): Number of PUs considered for the given exploration level
            - alignements: list of alignment objects

        Returns:
            - pu_path: str
                path to generated pdb file
        """
        # Get the order in which the PUs need to be merged to reconstruct the target
        # from 1 --> end
        pu_order = GraphPU.get_pu_order_from_target_ascending_pos(alignments)
        suffix = utils.get_random_name()
        base_path = f"{TMP_DIR}/{query.name}-level_{expl_level}_{nb_PUs}_PUs-on-{target.name}/"
        os.makedirs(base_path, exist_ok=True)
        merged_pus_path_renum = os.path.join(base_path, f"merged_pus_ali_{suffix}_renum.pdb")
        # Merge the positions of the PUs that aligned to corresponding ascending target positions
        merged_pus = []
        query_seq = ""
        for pu_idx in pu_order:
            ali = alignments[pu_idx]
            query_seq += ali.all_aligned["pu_seq"]
            for residue in ali.all_aligned["all_pu_aligned_positions"]:
                atoms = [atom for atom in ali.new_query_dict[residue]]
                merged_pus += atoms
        # Check if the PUs in the paths are consecutive ascending or descending
        # If sequential is True, we only keep the paths that are consecutive
        if sequential and query_seq != query.seq:
            pass
        else:
            merged_pus_pdb = "".join(merged_pus)
            with open(merged_pus_path_renum, "w") as filout:
                filout.write(merged_pus_pdb)
            # Renumber the merged PUs so that residues numbers correspond to the original PDB ones
            merged_pus_path = os.path.join(base_path, f"merged_pus_ali_{suffix}.pdb")
            utils.renum_ori_pdb_resnum(merged_pus_path_renum, ori_res_num_and_chain_query, new_path=merged_pus_path)
            # Renumbered PDB from 1 to length of the protein for easier PyMOL vizualisation
            Protein.reformat_struct(merged_pus_path_renum)  # reset ATOM and res sequence numbers
            paths.append((merged_pus_path, merged_pus_path_renum, alignments))
        return True

    def run_gdt2(scores, min_len_p1_p2, opt_prune, seed_alignment, target, query):
        """
        DO NOT PUT "self" TO THIS FUNCTION,
        IT SLOWS DOWN MULTIPROCESSING BY A LOT !

        Runs gdt2.pl with input query and target files.
        As gdt requires file to be in current folder, files are copied in
        current folder then the copy in deleted once gdt execution is over.

        Args:
            - query, Protein or PU:
                query structure used as an input for gdt.pl.
            - target, Protein or PU:
                target structure used as an input for gdt.pl.

        Returns:
            - score, float:
                TM-score normalized by length of shortest structure.
        """
        # Case when query and target have identical amino acid sequences
        # We need to treat the gdt alignment differently in this case
        if seed_alignment:
            mode = 1
        else:
            mode = 0
        # TM-align before calculating the scores with gdt2.pl
        ali = Alignment(query, target, opt_prune, seed_alignment)
        if ali.success:
            path_query_name = f"{TMP_DIR}/{query.name}"
            with open(path_query_name, "w") as filout:
                filout.write(ali.new_query)
            path_new_target_name = f"{TMP_DIR}/{target.name}_ali_{query.name}"
            with open(path_new_target_name, "w") as filout:
                filout.write(ali.new_target)
            command = f"{GDT} -pdb '{path_query_name} {path_new_target_name}' -mode {str(mode)} -len {min_len_p1_p2}"
            cmd_args = shlex.split(command)
            output = subprocess.run(cmd_args, capture_output=True, check=True)
            output = output.stdout.decode("utf-8")
            os.remove(path_new_target_name)
            os.remove(path_query_name)
            output = output.split("\n")
            for line in output:  # search output for TM-scores returned by gdt
                if "user length input" in line:
                    # retrieves scores as floats and return also the corresponding
                    # gdt2 output
                    scores[query.path] = (float(line[11:16]), output)
                    # skip next lines
                    break
            return True

    def compute_scores(self, min_len_p1_p2, nb_cpu, opt_prune, seed_alignment):
        """
        Computes all scores using the run_gdt2 function, select best alignment
        and its score and stores it in self.

        Args:
            - min_len_p1_p2 (int): minimum length of the two proteins
            - nb_cpu (int): number of cpu to use
            - opt_prune (float): threshold for pruning
            - seed_alignment (bool): if True, use a seed alignment

        Returns:
            - succeeded (bool): True if computation was successful, False otherwise
        """
        prots = []
        manager = multiprocessing.Manager()
        scores = manager.dict()

        # The pruning by TM-score can have removed all PU alignments, so we check that
        # there are still some alignments full paths to compute
        if len(self.final_alis) != 0:
            for path_tuple in self.final_alis:
                path = path_tuple[0]
                prots.append(Protein(path, peel=False))
            nb_prots = len(prots)

            # Multiprocessing
            chunksize = self.calc_chunksize(nb_cpu, nb_prots)
            with multiprocessing.Pool(processes=nb_cpu, initializer=signal.signal, initargs=(signal.SIGINT, signal.SIG_IGN)) as p:
                try:
                    func = partial(GraphPU.run_gdt2, scores, min_len_p1_p2, opt_prune, seed_alignment, self.target)
                    for i, _ in enumerate(p.imap_unordered(func, prots, chunksize)):
                        self.progressbar(i + 1, nb_prots, "Compute scores")
                    p.close()
                    p.join()
                except KeyboardInterrupt:
                    if os.path.exists(TMP_DIR):
                        shutil.rmtree(TMP_DIR, ignore_errors=True)
                        print("\nQuitting gracefully, bye !")
                    sys.exit(0)
            self.best_ali = max(scores.items(), key=lambda x: x[1][0])[0]
            self.best_score = scores[self.best_ali][0]
            self.best_gdt2_output = scores[self.best_ali][1]
            for tuple_paths_ali in self.final_alis:
                if self.best_ali == tuple_paths_ali[0]:
                    self.best_path = tuple_paths_ali[2]
                    self.best_ori_query_renum_pdb = tuple_paths_ali[1]
            return True
        elif len(self.final_alis) == 0 and opt_prune > 0.0:
            print(f"\nWARN: No solution found for this level after pruning alignments of PUs by TM-score of: {opt_prune}.\n==> Try with a lower threshold.")
        else:
            print(f"\nWARN: No solution found for this level")
        return False

    def write_output_PUs_query_regions(self, pu_order):
        """
        Write the output order of the query PU order and corresponding regions.
        
        Args:
            None
        
        Returns:
            None

        Parameters:
            - pu_order_text (string):
                example:
                    Target: d1mup__
                    Query: d1adl__
                    ├── PU order:     PU1 |     PU5 |     PU2 |     PU3 |     PU4 |
                    └── Regions :    1-22 |   23-37 |   38-64 |  65-105 | 106-131 |
        """
        # writing pu order
        output = f"    Query: {self.query.name}"
        output += "\n     ├── PU order: "
        for i, pu in enumerate(pu_order):
            output += "{:>10}".format(f"PU{i + 1} | ")

        # writing regions
        output += "\n     └── Regions : "
        for pu in pu_order:
            # Print the original residue numbers, not the renumerotation from 1
            output += "{:>10}".format(f"{self.ori_res_num_and_chain_query[self.best_path[pu].query.residues_num[0]][0]}-{self.ori_res_num_and_chain_query[self.best_path[pu].query.residues_num[-1]][0]} | ")

        output += f"\n    Target: {self.target.name}"
        output += f"\n     └── Sequence length : {self.target.length}"
        self.pu_order_text = output

    def _get_pu_order_from_query_ascending_positions(self):
        """
        Determine the order of PUs according to their ascending residues positions.
        """
        # Positions of the PUs
        pos = [ali.query.start for ali in self.best_path]
        # Keep in memory the order of the PUs that aligned to the target in the best_path
        pu_pos = [(pu_i, p) for pu_i, p in enumerate(pos)]
        # Sort the positions
        temp = sorted(pu_pos, key=lambda x: x[1])
        pu_order = [i[0] for i in temp]
        return pu_order

    @staticmethod
    def get_pu_order_from_target_ascending_pos(alignments):
        """
        Determine the order of PUs to write in the final PDB file, because we will need
        to write the PUs in the order of the corresponding ascending positions in the target.
        """
        # Smooth the positions of the target to avoid huge gaps messing with the right order of PUs
        smoothed_target_pos = []
        target_begin_pos = []
        for ali in alignments:
            smoothed_target_pos = ali.get_smoothed_target_positions(min_contiguous=5, n_gaps=2)
            # List of first positions of the target fragments that aligned with each PU
            target_begin_pos.append(smoothed_target_pos[0])
        # Reconstruct target from 1 -> end
        # For each alignment between a PU and the target we keep its number in
        # the order and the first position aligned
        idx_target_begin_pos = [(pu_i, p) for pu_i, p in enumerate(target_begin_pos)]
        # sort by position: second value of tuples
        sorted_target_positions = sorted(idx_target_begin_pos, key=lambda x: x[1])
        # Retrieve the pu_i of each sorted_target_positions
        # which also correspond to the index of PUs
        order = [i[0] for i in sorted_target_positions]
        return order


    def progressbar(self, n_done, n_todo, msg):
        """
        Simple progressbar to follow the graph's creation progress.

        Args:
            - n_done (int): Number of events passed
            - n_todo (int): Total number of events to pass
            - msg (str): Message to print before the progress bar

        Returns:
            - None
        """
        perc = n_done / n_todo
        done = round(30 * perc)
        print(f"{msg: <17}[{(done*'#')}{(30-done)*'.'}] [{perc:>4.0%}] \r",
              end="")

    def nb_ali(self, n):
        """
        Computes number of alignement to do if n PUs to align.

        Args:
            - n (int): Number of PUs to align.

        Returns:
            - total (int): Number of alignments to do.
        """
        total = 0
        liste = range(1, n + 1)
        for i, _ in enumerate(liste):
            total += utils.product(liste[i:])
        return total
