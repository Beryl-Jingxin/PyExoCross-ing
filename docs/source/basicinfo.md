# Basic information

Basic information includes data source and file path.

## Data source

In data source section, please provide:

1. The name of database, molecule, isotopologue and dataset.
2. The molecule ID and isotopologue ID.

**For ExoMol**

The name of `Database`, `Molecule`, `Isotopologue` and `Dataset` are necessary.

The molecule and isotopologue ID `MolIsoID` can be set as '0' or any other integers.

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
MolIsoID                                501
```

**For HITRAN**

The `Database` name, `Molecule` name and the molecule and isotopologue ID `MolIsoID` are necessary. The first two digits of `MolIsoID` are molcule ID and the third digit is isotopologue ID and there is no blank between molecule ID and isotopologue ID.

The name of `Isotopologue` and `Dataset` can be set as 'none' or any other words.

*Example*

```bash
# Data source #
Database                                HITRAN
Molecule                                CO2
Isotopologue                            none
Dataset                                 none
MolIsoID                                21
```

```bash
# Data source #
Database                                HITRAN
Molecule                                C2H2
Isotopologue                            none
Dataset                                 none
MolIsoID                                261
```

## File path

File path section records the file path for both reading and saving.

**For ExoMol**

`ReadPath` is the input files' folder path when the input line list files path is stored in the following format.

```
/FolderPath/molecule/iso-slug/dataset/iso-slug__dataset.states.bz2

/mnt/data/exomol/exomol3_data/MgH/24Mg-1H/XAB/24Mg-1H__XAB.states.bz2
```

```
└── exomol3_data
           ├── C2H2
           ├── CO2
           │     ├── 12C-16O2
           │     ├── 13C-16O2
           │     ├── ...
           │     ├── 12C-16O2__air.broad
           │     └── 12C-16O2__self.broad
           ├── MgH
           │     ├── 24Mg-1H
           │     │       ├── Yadin
           │     │       └── XAB
           │     │            ├── 24Mg-1H__XAB.def
           │     │            ├── 24Mg-1H__XAB.pf
           │     │            ├── 24Mg-1H__XAB.states.bz2
           │     │            └── 24Mg-1H__XAB.trans.bz2
           │     ├── 25Mg-1H
           │     │       ├── Yadin
           │     │       └── XAB
           │     │            ├── 25Mg-1H__XAB.def
           │     │            ├── 25Mg-1H__XAB.pf
           │     │            ├── 25Mg-1H__XAB.states.bz2
           │     │            └── 25Mg-1H__XAB.trans.bz2
           │     └── 26Mg-1H
           ├── ...
           │
```

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

*Example*

```bash
# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
```

**For HITRAN**

`ReadPath` is the file path of input line list (.par) file .

```
/home/username/data/hitran/CO2.par
```

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

*Example*

```bash
# File path #
ReadPath                                /home/jingxin/data/HITRAN/CO2.par
SavePath                                /home/jingxin/data/pyexocross/
```

## Functions

In current version, *PyExoCross* can convert data format between the ExoMol and the HITRAN formats. *PyExoCross* also implements the computations of other useful functions including partition functions, specific heats, cooling functions, radiative lifetimes, stick spectra and cross sections.

Use this function or not:

`0` means no

`1` means yes

If the value of a function's second column is `0`, then there is no need to do any changes in this function's own section, the program won't process data with this function. Although this function won't be used by users, please don't delete this function's own section, otherwise, the program cannot run.

*Example*

```bash
# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           1
CoolingFunctions                        0
Lifetimes                               0
StickSpectra                            0
CrossSections                           1
```

---

***Note***

1. Just change the information which you will use, please do not delete other information. Please do not change the first column strings.
2. Cooling functions' temperature starts from T = 200 K, the temperature of the other functions (partition functions, specific heats and lifetimes) starts from T = 1 K.
