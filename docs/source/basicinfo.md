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
            ├── CO2
            │    ├── 12C-16O2
            │    ├── 13C-16O2
            │    ├── ...
            │    ├── 12C-16O2__air.broad
            │    └── 12C-16O2__self.broad
            ├── C2H2
            ├── MgH
            │    ├── 24Mg-1H
            │    │      ├── Yadin
            │    │      └── XAB
            │    │         ├── 24Mg-1H__XAB.def
            │    │         ├── 24Mg-1H__XAB.pf
            │    │         ├── 24Mg-1H__XAB.states.bz2
            │    │         └── 24Mg-1H__XAB.trans.bz2
            │    ├── 25Mg-1H
            │    │      ├── Yadin
            │    │      └── XAB
            │    │         ├── 25Mg-1H__XAB.def
            │    │         ├── 25Mg-1H__XAB.pf
            │    │         ├── 25Mg-1H__XAB.states.bz2
            │    │         └── 25Mg-1H__XAB.trans.bz2
            │    └── 26Mg-1H
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
