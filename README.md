<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://user-images.githubusercontent.com/25644865/195085417-51ecbae0-2722-49d2-8603-2d0a1cbe1f9f.png" width="300">
  <img alt="" src="https://user-images.githubusercontent.com/25644865/195085445-b3af5175-8c61-4710-847c-df8907cb7617.png" width="300">
</picture>
</p>

# Icarus: Flexible protein structural alignment based on Protein Units


[![DOI](https://zenodo.org/badge/549521192.svg)](https://zenodo.org/badge/latestdoi/549521192)

Paper: [Cretin, G., Périn, C., Zimmermann, N., Galochkina, T., & Gelly, J. C. (2023). ICARUS: flexible protein structural alignment based on Protein Units. Bioinformatics, 39(8), btad459.](https://doi.org/10.1093/bioinformatics/btad459)

Icarus is a method which uses the Protein Peeling algorithm (Gelly et al. (2006a), Gelly et al. (2006b), Gelly et al. (2011), Postic et al. (2017), Cretin et al. (2022)) to identify compact regions i.e Protein Units (PUs). PUs define rigid regions to be aligned to the target and delimit hinge positions in the structure.  
Protein Peeling allows a hierarchical segmentation of a protein into compact "independent" domains (i.e that maximise intra-domain contact while minimizing inter-domain contact).  
A protein can be divided into different exploration levels, each level containing more and more PUs as the level rises. The user can choose between different exploration levels:  
| Exploration level | Max number of Protein Units to consider (up to)  |
|:-----------------:|--------------------------------------------------|
|         1         |        2 and/or 3                                |
|       **2**       | **4 and/or 5 (default)**                         |
|         3         |            6                                     |
|         4         |            7                                     |

If user asks for an exploration of 4, the program will explore **up to** 7 PUs,  
meaning that if the best alignment is found at level 3, it will give this best alignment anyways.
You can set --verbose mode to have detailed output with all intermediate alignments.

When given a pair of proteins to align, the program will generate several outputs:  

**Default output**  

* One PDB representing the best *query* protein structure transformed after alignment.
* One PDB representing the best *query* protein structure transformed after alignment aligned against the original target protein.
* The textual output (scores and textual alignment) is saved into "summary.txt" file.

**Verbose output**  

* All intermediate alignments are kept, for each exploration level up to the one chosen by the user.
* For each exploration level, the best alignments are kept in a `results_PDB` directory, containing both the PDB files of the query alone and the query aligned against the target.


## Dependencies

First, you need to run `install.sh`.
This will essentially compile and install SWORD and its dependencies.

```bash
./install.sh
```

### Conda

Create an environment using the `environment.yml` file in the git repository which contains all dependencies.

```bash
# Create the environment
conda env create -f environment.yml
or
mamba env create -f environment.yml
# Activate the environment
conda activate icarus
```

## Usage

```
$ ./icarus.py --help
usage: icarus.py [-h] -p1 PROTEIN1 -p2 PROTEIN2 [-m MIN_SIZE] [-c1 CHAIN1] [-c2 CHAIN2] [-l EXPLORATION_LEVEL] [-p PRUNE] [-s | -n] [-t] [-e] [-f] [-c CPU] [-v] [-u]

Icarus is a flexible structural alignment program.
It takes as input 2 pdb files and returns the optimal
alignment based on different protein exploration levels.

required arguments:
  -p1 PROTEIN1, --protein1 PROTEIN1
                        Path to the first protein to align
  -p2 PROTEIN2, --protein2 PROTEIN2
                        Path to the second protein to align

optional arguments:
  -m MIN_SIZE, --min-size MIN_SIZE
                        Minimum size of Protein Units (PUs).
                        Must be 15 <= min <= 99, default 15
  -c1 CHAIN1, --chain1 CHAIN1
                        PDB Chain of the query.
                        Single alphabetical character [A-B].
                        Default is chain A.
  -c2 CHAIN2, --chain2 CHAIN2
                        PDB Chain of the target.
                        Single alphabetical character [A-B].
                        Default is chain A.
  -l EXPLORATION_LEVEL, --exploration-level EXPLORATION_LEVEL
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
                          4 -> 7
  -p PRUNE, --prune PRUNE
                        The pruning threshold corresponds to a TM-score value
                        used to filter the KPAX between Protein Units and the target protein.
                        If an alignment is below this threshold, it will be pruned from the graph.
                        A high pruning threshold (TM-score value) will filter out more solutions
                        whereas a low one will prune less solutions.
                        Default value is no pruning (0.).
  -s, --seed-alignment  In the case when both input PDBs have identical amino acid sequences
                        but differ in 3D, setting this option will make sure that
                        all alignments by KPAX will be assigned to avoid any displacement.
                        This is useful for instance in the case of an analysis of structures inside a
                        dynamics simulations. Default is not set
  -n, --no-seed-alignment
                        Force the flexible alignment of identical amino acid sequences
                        without setting a seed alignment for KPAX. Default is not set
  -t, --smoothed-pu-output
                        If this option is set, the Protein Units (PUs) that were aligned to the target protein
                        are smoothed / trimmed to keep only the "core" aligned position.
                        This option only changes the final presentation of the textual alignment for "visual"
                        purposes. The output PDB files contain all the residues. Default is not set
  -e, --sequential      If this option is set, the program considers only solutions with consecutive Protein Units.
                        Either ascending or descending order.
  -f, --force           Bypass asking user confirmation for exploration level >= 4
  -c CPU, --cpu CPU     How many CPUs to use. Default all (0). Max on this computer is: 32
  -v, --verbose         Set verbose mode: print longer output and generate
                                          intermediate results and alignments
  -u, --use-ramfs       Use the native TMPFS partition /dev/shm on linux systems to boost i/o perfs.
                        This may use a large amount of RAM if proteins are large e.g. > 200 residues

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
See the README for detailed information on verbose output
```

## Toy example

```bash
./icarus.py -p1 data/RIPC/d1adl__.ent -p2 data/RIPC/d1mup__.ent
```

```
Clean input PDB files ... done

Peel d1adl__:
    131 aa
    Seq: CDAFVGTWKLVSSENFDDYMKEVGVGFATRKVAGMAKPNMIISVNGDLVTIRSESTFKNTEISFKLGVEFDEITADDRKVKSIITLDGGALVQVQKWDGKSTTIKRKRDGDKLVVECVMKGVTSTRVYERA
Peel d1mup__:
    157 aa
    Seq: EEASSTGRNFNVEKINGEWHTIILASDKREKIEDNGNFRLFLEQIHVLENSLVLKFHTVRDEECSELSMVADKTEKAGEYSVTYDGFNTFTIPKTDYDNFLMAHLINEKDGETFQLMGLYGREPDLSSDIKERFAQLCEEHGILRENIIDLSNANRC


Aligning d1adl__ against d1mup__

SOLUTION 1: d1adl__ [level 1 => 3 PUs] vs d1mup__
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

SOLUTION 2: d1adl__ [level 2 => 4 PUs] vs d1mup__
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

SOLUTION 3: d1adl__ [level 2 => 5 PUs] vs d1mup__
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

Aligning d1mup__ against d1adl__

SOLUTION 4: d1mup__ [level 1 => 3 PUs] vs d1adl__
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

SOLUTION 5: d1mup__ [level 2 => 4 PUs] vs d1adl__
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

SOLUTION 6: d1mup__ [level 2 => 5 PUs] vs d1adl__
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]


**************************************************************************************************************
*************************************************** RESULTS **************************************************
**************************************************************************************************************


There is only 1 optimal solution for this alignment:


                                                 SOLUTION 3
                                                 ==========

 *  Score: 0.719
    Query: d1adl__
     ├── PU order:     PU1 |     PU2 |     PU3 |     PU4 |     PU5 |
     └── Regions :    1-22 |   23-37 |   38-64 |  65-105 | 106-131 |
    Target: d1mup__
     └── Sequence length : 157

                                                ALIGNED PU(S)
                                                =============

PU 1     :CDAFVGTWKLVSSEN--FDDY---MKE
          ||:||||||||||||   :.
TARGET   :VEKINGEWHTIILASDKREKIEDNGNF
ali. pos. 12                        38
ori. pos. 1                         22

PU 2     :VGVGFATRKVAGMAK
          :|||||||||||||:
TARGET   :PDLSSDIKERFAQLC
ali. pos. 130           144
ori. pos. 23            37

PU 3     :PNMIISVNGDLVTIRSES----TFKNTEISF
          ..:|||||||||||||:.    .:||||||:
TARGET   :FLEQIHVLENSLVLKFHTVRDEECSELSMVA
ali. pos. 41                            71
ori. pos. 38                            64

PU 4     :K-----LGVEFDEITADDRKVKSIITLDGG------ALVQVQKWDGKSTTIK
          .      |:||||      .||||||||::      ||||||||||||||||
TARGET   :DKTEKAGEYSVTY------DGFNTFTIPKTDYDNFLMAHLINEKDGETFQLM
ali. pos. 72                                                 123
ori. pos. 65                                                 105

PU 5     :RKRDGDKLVVECVMKGVT-STRVY---ERA
                       .: .  .||||   :|:
TARGET   :-----------EEHGILRENIIDLSNANRC
ali. pos. 145                          174
ori. pos. 106                          131

                                                BEST ALIGNMENT
                                                ==============

PUs      :                      PU1                            PU3                                      PU4
ori. pos.:           ┌1                      22┐  ┌38                         64┐┌65
connect  :           +-------------------------+  +-----------------------------++----------------------------
QUERY    :           CDAFVGTWKLVSSEN--FDDY---MKE--PNMIISVNGDLVTIRSES----TFKNTEISFK-----LGVEFDEITADDRKVKSIITLDG
match    :           ||:||||||||||||   :.         ..:|||||||||||||:.    .:||||||:.      |:||||      .||||||||:
TARGET   :EEASSTGRNFNVEKINGEWHTIILASDKREKIEDNGNFRLFLEQIHVLENSLVLKFHTVRDEECSELSMVADKTEKAGEYSVTY------DGFNTFTIPK
dist     :           112000101000001  5248   79   342111010010111024    4201110023     5121100      3000110002
ali. pos.:         10        20        30        40        50        60        70        80        90

PUs      :                                  PU2                   PU5
ori. pos.:                   105┐      ┌23         37┐┌106                      131┐
connect  :----------------------+      +-------------++----------------------------+
QUERY    :G------ALVQVQKWDGKSTTIK------VGVGFATRKVAGMAKRKRDGDKLVVECVMKGVT-STRVY---ERA
match    ::      ||||||||||||||||      :|||||||||||||:             .: .  .||||   :|:
TARGET   :TDYDNFLMAHLINEKDGETFQLMGLYGREPDLSSDIKERFAQLC-----------EEHGILRENIIDLSNANRC
dist     :2      1001100010100000      200000000001102           7532635 31111   212
ali. pos.:         110       120       130       140       150       160       170

Aligned distance (match <=> dist): '|' <= 1 Å
                                   ':' <= 2 Å
                                   '.' <= 3 Å

Best solution(s):
--> /home/republique/cretin/PROJECTS/icarus-github/icarus_output/results/d1adl___and_d1mup__/solution_3_d1adl__-level_2_5_PUs-on-d1mup__.pdb

Total runtime: 30.0 seconds
```

## ICARUS output directories/files explained:

### Non-verbose output  

Contains the best transformed query PDB structure as a PDB file  
and the textual output (best score(s) and corresponding textual alignment(s))
```
work/results/query_and_target:
 ├── solution1_query-on-target-level_X_N_PUs.pdb
 └── summary.txt (terminal textual output)
```

### Verbose output  

Contains results for protein1 vs. protein2 and vice versa,
and summary.txt is the terminal textual output.
```
work/results/d1adl___and_d1mup__:
 ├── d1adl___on_d1mup__/
 ├── d1mup___on_d1adl__/
 └── summary.txt
```

Contains all intermediate alignments with KPAX outputs,
the ICARUS final result PDB and the simple protein1 vs. protein2
KPAX output for information.
```
work/results/d1adl___and_d1mup__/d1adl__-on-d1mup__:
 ├── intermediate/
 ├── result_PDBs/
 └── gdt_d1adl___vs_d1mup__.txt
```

Contains the best intermediate alignments that were done by ICARUS
during runtime while graph exploration.
```
work/results/d1adl___and_d1mup__/d1adl__-on-d1mup__/intermediate:
 ├── alignements_level_1_3_PUs/
 ├── alignements_level_2_4_PUs/
 └── alignements_level_2_5_PUs/
```

Each directory contains the TM-sup and TM-sup_all_atm output PDB files of
the best KPAXs done between best alignment PUs and the target to which
previously aligned portions of PUs were removed.
```
work/results/d1ggga__and_d1wdna_/d1adl__-on-d1mup__/intermediate/alignements_level_1_3_PUs:
 ├── ali_d1adl__-d1adl__.ent_1_106_131-on-d1mup__/
 ├── ali_d1adl__-d1adl__.ent_1_1_64-on--4SF3zUEzYS/
 └── ali_d1adl__-d1adl__.ent_1_65_105-on--cRHbY4G4OX/
```

Contains all the solutions found by ICARUS at each exploration level
and for each number of PUs of these levels.
Are given:
1 - query-exploration_level_and_number_of_PUs-on-target.pdb
2 - query-exploration_level_and_number_of_PUs.pdb
```
work/results/d1ggga__and_d1wdna_/d1adl__-on-d1mup__/intermediate/result_PDBs:
 ├── d1adl__-level_1_3_PUs-on-d1mup__.pdb
 ├── d1adl__-level_1_3_PUs.pdb
 ├── d1adl__-level_2_4_PUs-on-d1mup__.pdb
 ├── d1adl__-level_2_4_PUs.pdb
 ├── d1adl__-level_2_5_PUs-on-d1mup__.pdb
 ├── d1adl__-level_2_5_PUs.pdb
 └── d1adl___on_d1mup___kpax.pdb
```
