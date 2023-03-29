<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://user-images.githubusercontent.com/25644865/195085417-51ecbae0-2722-49d2-8603-2d0a1cbe1f9f.png" width="300">
  <img alt="" src="https://user-images.githubusercontent.com/25644865/195085445-b3af5175-8c61-4710-847c-df8907cb7617.png" width="300">
</picture>
</p>

# Icarus: Flexible protein structural alignment based on Protein Units




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


## Docker (recommended) - Linux, MacOS (Intel) and Windows

Icarus is available as a Docker image (143.42 MB compressed, 391MB MB on disk).  
You can either pull the latest image from Dockerhub:  
```
docker pull dsimb/icarus
```

Run the image:  
```
# Set internal /dev/shm size as your RAM size, it is used to write tmp files for faster processing.
# The space used depends on several parameters, but a good value generally is >=1 GB.
# If an error occures about lack of space on device, increase this value.
# For convenience the following command lines (Linux or Mac) retrieve the max value of RAM available
# For windows users, please set shm variable manually: shm!<int>
shm=$(free -g | awk '/^Mem:/{print $2}')gb # linux
shm=$(system_profiler SPHardwareDataType | grep "Memory:" | awk '{print $2}')gb # MacOS
mkdir icarus_output
docker run -it --shm-size $shm -v $(pwd)/icarus_output:/icarus/icarus_output -v ./data:/data dsimb/icarus -p1 /data/RIPC/d1adl__.pdb -p2 /data/RIPC/d1mup__.pdb

# Show help
docker run dsimb/icarus
```

### Docker - MacOS (Apple Silicon M1)  

Docker Desktop >= 3.3.1 is required.
For the ARM64 architecture of Apple Silicon M1(+) CPU, there is an additional option to add to all docker command: `--platform linux/arm64`:  

```
# Build the image and tag it icarus
docker build --platform linux/arm64 -t icarus .

# Run the image
docker run --platform linux/arm64 -it --shm-size $shm -v $(pwd)/icarus_output:/icarus/icarus_output -v ./data:/data dsimb/icarus -p1 /data/RIPC/d1adl__.pdb -p2 /data/RIPC/d1mup__.pdb
```

## Dependencies

First, you need to run the installer in order to deploy the workspace.
This will essentially compile and install SWORD and its dependencies,  
TMalign, and create a few directories to work.

```bash
./install.sh
```

### Conda

If you are familiar with conda, you should create an environment using the `environment.yml` file in the git repository which contains all dependencies.

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
usage: icarus.py [-h] -p1 PROTEIN1 -p2 PROTEIN2 [-s MIN_SIZE] [-l EXPLORATION_LEVEL] [-f] [-v]

Icarus is a flexible structural alignment program.
It takes as input 2 pdb files and returns the optimal
alignment based on different protein exploration levels.

optional arguments:
  -h, --help            show this help message and exit
  -p1 PROTEIN1, --protein1 PROTEIN1
                        Path to the first protein to align
  -p2 PROTEIN2, --protein2 PROTEIN2
                        Path to the second protein to align
  -s MIN_SIZE, --min-size MIN_SIZE
                        Minimum size of Protein Units (PUs).
                        Must be 15 <= min <= 99, default 15
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
                          4 -> 7,
                          5 -> 8
  -f, --force           Bypass asking user confirmation for exploration level >= 4
  -v, --verbose         Set verbose mode: print longer output and generate
                                          intermediate results and alignments

Explanations of ICARUS output:

Non-verbose output
------------------
work/results/query_and_target:
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

Peel d1adl__.ent:
    131 aa
    Seq: CDAFVGTWKLVSSENFDDYMKEVGVGFATRKVAGMAKPNMIISVNGDLVTIRSESTFKNTEISFKLGVEFDEITADDRKVKSIITLDGGALVQVQKWDGKSTTIKRKRDGDKLVVECVMKGVTSTRVYERA
Peel d1mup__.ent:
    157 aa
    Seq: EEASSTGRNFNVEKINGEWHTIILASDKREKIEDNGNFRLFLEQIHVLENSLVLKFHTVRDEECSELSMVADKTEKAGEYSVTYDGFNTFTIPKTDYDNFLMAHLINEKDGETFQLMGLYGREPDLSSDIKERFAQLCEEHGILRENIIDLSNANRC




d1adl__ level 1 vs d1mup__, 3 PUs
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

d1adl__ level 2 vs d1mup__, 4 PUs
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

d1adl__ level 2 vs d1mup__, 5 PUs
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

d1mup__ level 1 vs d1adl__, 3 PUs
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

d1mup__ level 2 vs d1adl__, 4 PUs
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

d1mup__ level 2 vs d1adl__, 5 PUs
--------------------------------------------------------
Build graph      [##############################] [100%]
Merge alignments [##############################] [100%]
Compute scores   [##############################] [100%]

INFO: Overwriting existing results at /home/republique/cretin/PROJECTS/icarus/icarus_output/results/d1adl___and_d1mup__




**************************************************************************************************************
*************************************************** RESULTS **************************************************
**************************************************************************************************************


There is only 1 optimal solution for this alignment:


                                                 SOLUTION 3
                                                 **********

 *  Score: 0.721
    Query: d1adl__
     ├── PU order:     PU1 |     PU2 |     PU3 |     PU4 |     PU5 |
     └── Regions :    1-22 |   23-37 |   38-64 |  65-105 | 106-131 |
    Target: d1mup__
     └── Sequence length : 157

                                                ALIGNED PU(S)
                                                -------------

PU 1     :CDAFVGTWKLVSSEN--F-DD---YMKE
          :||||||||||||||    |
TARGET   :VEKINGEWHTIILASDKREKIEDNGNFR
ali. pos. 12                         39
ori. pos. 1                          22

PU 2     :VGVGFATRKVAGMAK
          .|||||||||||||.
TARGET   :PDLSSDIKERFAQLC
ali. pos. 126           140
ori. pos. 23            37

PU 3     :PNMIISVNGDLVTIRSES----TFKNTEISF
            ||||||||||||||.      .||||||:
TARGET   :FLEQIHVLENSLVLKFHTVRDEECSELSMVA
ali. pos. 41                            71
ori. pos. 38                            64

PU 4     :K--LGVEFDEITADDRKVKSII-TLDG-GALVQVQKWD----GKSTTIK
               |||||:.  ::|||||| .|.  .|||||||:.    .||||||
TARGET   :KTEKAGEYSVTY--DGFNTFTIPKTDYDNFLMAHLINEKDGETFQLMGL
ali. pos. 73                                              121
ori. pos. 65                                              105

PU 5     :RKRDGDKLVVECVMKGVT--STRV--YERA
                       .|.|.  |::|  :.||
TARGET   :-----------EEHGILRENIIDLSNANRC
ali. pos. 141                          170
ori. pos. 106                          131

                                                BEST ALIGNMENT
                                                --------------

PUs      :                      PU1                            PU3                                      PU4
ori. pos.:           ┌1                       22┐ ┌38                         64┐ ┌65
connect  :           +--------------------------+ +-----------------------------+ +---------------------------
QUERY    :           CDAFVGTWKLVSSEN--F-DD---YMKE-PNMIISVNGDLVTIRSES----TFKNTEISF-K--LGVEFDEITADDRKVKSII-TLDG-
match    :           :||||||||||||||    |           ||||||||||||||.      .||||||:      |||||:.  ::|||||| .|.
TARGET   :EEASSTGRNFNVEKINGEWHTIILASDKREKIEDNGNFRLFLEQIHVLENSLVLKFHTVRDEECSELSMVADKTEKAGEYSVTY--DGFNTFTIPKTDYD
dist     :           211011000000000  6 07   968  550000110010000146    531000002 7   60100123  22111001 3048
ali. pos.:         10        20        30        40        50        60        70        80        90

PUs      :                              PU2                   PU5
ori. pos.:                 105┐    ┌23         37┐┌106                      131┐
connect  :--------------------+    +-------------++----------------------------+
QUERY    :GALVQVQKWD----GKSTTIK----VGVGFATRKVAGMAKRKRDGDKLVVECVMKGVT--STRV--YERA
match    :.|||||||:.    .||||||    .|||||||||||||.             .|.|.  |::|  :.||
TARGET   :NFLMAHLINEKDGETFQLMGLYGREPDLSSDIKERFAQLC-----------EEHGILRENIIDLSNANRC
dist     :3100000023    3100010    311000000000003           9731413  0221  2310
ali. pos.:         110       120       130       140       150       160       170

Aligned distance (match <=> dist): '|' <= 1 Å
                                   ':' <= 2 Å
                                   '.' <= 3 Å

Best solution(s):
--> /home/republique/cretin/PROJECTS/icarus/icarus_output/results/d1adl___and_d1mup__/solution_1_d1adl__-level_2_5_PUs-on-d1mup__.pdb

Total runtime: 26.3 seconds
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
