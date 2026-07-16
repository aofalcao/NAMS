# NAMS
Non-Contiguous atom molecular structural similarity

NAMS is a molecular similarity computation program. It needs that the molecular data is stored in a specific formmat. This is accomplihed by using the supplied `makenamsdb.py` utility that encapsulates all the chemical knowledge about a molecule or molecule database

## Running `makenamsdb`

makenamsdb is a sophisticated program that allows many possibilities but its main goal is to create files for `nams`
running it asking for help gives the following output

```
python makenamsdb.py --help
usage: makenamsdb.py [-h] [--input-format INPUT_FORMAT] [--binary] [--stereo] [--explicit-hydrogens]
                     [--fragment-policy {legacy-first,largest}] [--delimiter DELIMITER]
                     [--smiles-column SMILES_COLUMN] [--id-column ID_COLUMN] [--skip-header] [--write-images]
                     [--svg-directory SVG_DIRECTORY] [--strict] [--quiet] [--manifest MANIFEST]
                     input output

Build a NAMS text or legacy binary database using Open Babel 3.

positional arguments:
  input                 Delimited input molecular file
  output                Output NAMS database

options:
  -h, --help            show this help message and exit
  --input-format INPUT_FORMAT
                        Open Babel input format for the molecular column (default: smi)
  --binary              Write legacy compact binary
  --stereo              Encode atom chirality and double-bond E/Z (text output recommended)
  --explicit-hydrogens  Add explicit hydrogens to bond environments
  --fragment-policy {legacy-first,largest}
                        Disconnected-component policy (default: legacy-first)
  --delimiter DELIMITER
                        Column delimiter; escape sequences accepted, e.g. '\t'
  --smiles-column SMILES_COLUMN
                        One-based molecular-text column (default: 1)
  --id-column ID_COLUMN
                        One-based integer-ID column (default: 2)
  --skip-header
  --write-images        Write indexed SVG files
  --svg-directory SVG_DIRECTORY
  --strict              Stop on the first invalid record
  --quiet
  --manifest MANIFEST   Manifest path (default: OUTPUT.manifest.json)

```

Please make sure [OpenBabel](http://openbabel.org/docs/Installation/install.html) is properly installed


## Creating a simple database

ani SMILEs file (SMILES encoding and Numerical ID) can be transformed into a NAMS file simply by running makenamsdb without any other options

```
$ python makenamsdb.py data\aminoacids.smi data\aminoacids.nams
$ python makenamsdb.py data\ibuprofen.smi data\ibuprofen.nams
$ python makenamsdb.py data\aspirin.smi data\aspirin.nams
```
The above will create 3 databases, 2 with 1 molecule (for aspirin and ibuprofen and the other for the 20  regular amino acids

## Using NAMS

### NAMS Options

calling the help file will give all the available NAMS parameters:

```
$ nams -help 
NAMS - Non-contiguous atom matching structural similarity function
nams version - 0.99.20260716
This program will compare one or several molecules against others producing
a similarity score for each pair of molecules
Basic use:
nams -mode [mode] -in [mol_filename] -db [mod db filename]
NOTE: -mode, and -db are required options
Details
    -in [filename] filename is the name of a file in the NAMS format with one
       or more molecules
    -db [filename] filename is the name of a file in the NAMS format with one
        or more molecules
    -mode: Operation mode - must follow a value: 1,2,3 or 4
        1 - Compares one molecule (given with the -in option) with a full
            database (-db)
        2 - Compares a database with itself. All molecules of db are compared
            to all the others in the database
        3 - Compares a list of molecules (-in) with all the molecules in (-db)
    -opts [AJLMS] Output options. Each of these 5 letters or combinations may
            be used. If options A or M are absent results are provided in a table
            [default -JLS - if -opts present all is reset]
        M - Presents the atom matching scores for all atom pairs for both
            molecules
        A - Gives the optimum matching between atoms of both molecules.
            Each atom pair is identified by their atom ID and the
            corresponding matching score
        L - Self similarity scores of both molecules [default]
        S - Similarity score between both molecules [default]
        J - Jaccard score between both molecules [default]
Filter parameters:
    -tj [value] - presents only the results above a given Jaccard threshold
    -fdb [filename] -defines a filename that includes the molids to search within
             the database, as specified by -db
    -fin [filename] -defines a filename that includes the molids to search within
             the NAMS input file, as specified by -in
    -fmxmw [value] - the maximum molecular %weight difference to search within the
             database [default=9999]
    -fmnmw [value] - the minimum molecular %weight difference to search within the
             database [default=100]
NAMS Heuristic parameters:
    -alpha [value] - alpha parameter of the model [default=0.9]
    -nrings [value] - score that acounts for how many rings an atom belongs
           [default=0.80]
    -stiso [value] - parameter for accountin chiral stereo isomerism
           [default=0.95]
    -ctiso [value] - accounts for double bond cis-trans isomerisms
           [default=0.95]
    -bring [value] - whether a bond is in a ring or not [default=0.90]
    -barom [value] - whether a bond is in an aromatic ring or not
           [default = 0.90]
    -border [value] - bond covalent order [default=0.90]
    -pen [value] - penalty factor to account for unmatched bonds
           [default= 1.0]
    -ADM [code] - atom distance matrix code. Admissible values are 0,1,2, 3 and 4
           [default=3]
(c) Andre Falcao 2013-2026
    This Program is provided free of charge and WITHOUT ANY WARRANTY
  if used please cite:
     AL Teixeira, AO Falcao. 2013. Non contiguous atom matching structural similarity
     function DOI: 10.1021/ci400324u. Journal of Chemical Information and Modeling
```
Despite its dounting options there are essentially 3 modes of Using NAMS

* `-mode 1` - Compares one molecule to a database and alseo another allowing very detailed output 
* `-mode 2` - Compares a database to itself, producing a table with the similarity of each molecule to every other one (essentiaaly constructing a similarities matrix)
* `-mode 3` - Compares 2 databases

### MODE 1 - Comparing 2 molecules

One can compute 2 molecules alone directly using option 1
```
$ nams -in data\aspirin.nams -db data\ibuprofen.nams -mode 1
```
produces the following output that can be redirected to a file
```
    molid1      molid2     ss1     ss2   sim12   Sscore
         1         2 134.209 168.889 100.774  0.4981
Number of molecular alignments performed: 1
Time required for execution: 0.002000 seconds
```
here we can see the self similarity score of 2 molecules their common similarity and the Jaccard score that specifies how similar each one is to the other. a socre of 1.0 implies identity and very close isomers. NAMS can differentiate several types of isomerism already.

We can get a more through output by using output options `-JALMS` that produces for each heavy atom of the molecule what is its similarity score to all the other atoms of the molecule, thus

```
$ nams -in data\aspirin.nams -db data\ibuprofen.nams -mode 1 -opts JALMS
Atom matching matrix
         1      2      3      4      5      6      7      8      9     10     11     12     13     14     15
  1   6.20   6.09   6.06   5.69   5.32   5.42   5.90   5.42   5.32   6.84   6.86   6.62   6.68   6.70   6.20
  2   5.77   6.75   6.63   6.40   6.38   6.27   6.62   6.27   6.38   7.32   7.96   6.39   6.46   7.26   5.77
  3   6.10   5.91   5.96   5.77   5.38   5.38   5.70   5.38   5.38   6.17   6.35   7.27   7.08   6.36   6.10
  4   5.38   6.10   7.26   7.06   7.28   7.55   7.49   7.55   7.28   7.85   6.83   6.58   6.77   6.87   5.38
  5   5.01   5.63   6.74   8.03   8.29   8.55   8.67   8.55   8.29   7.31   6.55   5.96   6.04   6.70   5.01
  6   4.81   5.46   6.55   7.80   8.29   8.55   8.08   8.55   8.29   6.98   6.03   5.54   5.53   6.16   4.81
  7   4.65   5.38   6.61   7.98   8.10   8.04   8.06   8.04   8.10   6.77   5.64   5.08   5.08   5.86   4.65
  8   4.67   5.49   6.65   7.80   7.97   8.30   8.09   8.30   7.97   6.90   5.90   5.03   5.07   6.14   4.67
  9   4.83   5.49   6.65   7.97   8.24   8.41   8.58   8.41   8.24   7.22   6.23   5.46   5.52   6.53   4.83
 10   4.89   5.67   6.72   8.08   8.15   8.41   8.57   8.41   8.15   7.70   6.74   5.78   5.85   6.58   4.89
 11   5.12   5.96   7.25   6.79   6.72   6.82   7.18   6.82   6.72   8.07   7.39   5.81   5.86   6.69   5.12
 12   5.36   6.49   6.33   6.04   5.85   6.00   6.30   6.00   5.85   6.85   7.04   6.74   6.55   7.05   5.36
 13   5.40   6.52   6.36   6.04   5.84   6.06   6.36   6.06   5.84   6.91   7.10   6.55   6.74   7.11   5.40
Atom matching and scores
  1  13 ->  6.68
  2  11 ->  7.96
  3  12 ->  7.27
  4   3 ->  7.26
  5   7 ->  8.67
  6   8 ->  8.55
  7   5 ->  8.10
  8   6 ->  8.30
  9   9 ->  8.24
 10   4 ->  8.08
 11  10 ->  8.07
 12   2 ->  6.49
 13  14 ->  7.11

Self Similarity:    1 -> 134.209
Self Similarity:    2 -> 168.889
Similarity: (1 2) -> 100.774
Jaccard Score: (1 2) -> 0.49809
Number of molecular alignments performed: 1
Time required for execution: 0.009000 seconds
```

### MODE 1 - comparing one molecule to a database

Mode 1 still allow us to get a table for comparing one molecule (ibuprofen) against all amino acids
```
Snams -in data\ibuprofen.nams -db data\aminoacids.nams -mode 1
    molid1      molid2     ss1     ss2   sim12   Sscore
         2         1 168.889  27.079  21.960  0.1262
         2         2 168.889 100.559  54.165  0.2516
         2         3 168.889  60.439  37.182  0.1935
         2         4 168.889  60.439  36.788  0.1911
         2         5 168.889  36.939  25.018  0.1384
         2         6 168.889  73.139  45.570  0.2320
         2         7 168.889  73.139  45.208  0.2297
         2         8 168.889  18.319  14.133  0.0817
         2         9 168.889  97.169  58.974  0.2848
         2        10 168.889  61.039  44.712  0.2414
         2        11 168.889  60.439  44.370  0.2399
         2        12 168.889  72.179  50.166  0.2628
         2        13 168.889  59.579  35.774  0.1857
         2        14 168.889 113.769 101.187  0.5576
         2        15 168.889  55.179  36.906  0.1972
         2        16 168.889  36.939  25.298  0.1401
         2        17 168.889  48.579  32.818  0.1777
         2        18 168.889 183.359 120.515  0.5201
         2        19 168.889 131.169 109.468  0.5744
         2        20 168.889  48.579  35.813  0.1971
Number of molecular alignments performed: 20
Time required for execution: 0.017000 seconds
```
we can see that the closest symilarity is to amino acid 19 (Typrosine) as expected, even though Phenylalanine and Tryptophan are also close


# Older Tutorials using R and RStudio

Please check the and the original source code release, along with some user tutorials here:

http://www.di.fc.ul.pt/~afalcao/NAMS.html



# Reference

For further information on what is NAMS and how to use it, check the original publication. and cite it if you find it useful 

https://pubs.acs.org/doi/abs/10.1021/ci400324u 

