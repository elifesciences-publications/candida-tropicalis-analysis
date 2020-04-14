# Analysis of _Candida tropicalis_ 3C-seq data
This repository includes codes used for analyzing 3C-seq data of _Candida tropicalis_ type strain MYA-3404. Manuscript available at https://www.biorxiv.org/content/10.1101/2020.02.07.938175v1. The codes can be applied for Hi-C/3C-seq data of other organisms as well.
### Dependencies
To use these scripts, the following Python libraries and their respective dependencies are required.
- numpy
- scipy
- matplotlib
- mirnylib
- hiclib
### Dowload
```
git clone https://github.com/Yao-Chen/candida-tropicalis-analysis.git
cd candida-tropicalis-analysis
```
### Input data
All the scripts are based on the matrix **HDF5** files exported by function `export` in class `binnedData` of `hiclib` package.

Prior to the matrix export, data processing by `hiclib` package, including iterative mapping of Hi-C/3C-seq raw reads, fragment assignment and filtering, binning and bin-level filtering and iterative bias correction, are recommended.

The convertion of genome-wide **interaction matrix** to **contact probability matrix** can be performed by normalizing the sum of entire matrix to the number of rows/columns, so that the sum of each row/column approximates 1 and each value _Cij_ represents the probability of contacts between bin _i_ and bin _j_.
```
>>> matrix = BD.dataDict[NAME]
>>> totalsum = sum(sum(matrix))
>>> BD.dataDict[NAME] = matrix / float(totalsum) * len(matrix)
```
Next, genome-wide and chromosome-wide contact probability matrices can be exported to HDF5 format and used as input for the codes in this repository.
```
>>> BD.export(NAME, GENOME_HDF5) # export genome-wide matrix
>>> BD.export(NAME, CHROMOSOME_HDF5, byChromosome="all") # export chromosome-wide matrices
```
### Usage
Please refer to the usage description in each script.
