# Analysis of Candida tropicalis 3C-seq data
This repository includes codes used for analyzing 3C-seq data of Candida tropicalis type strain MYA-3404. Manuscript available at https://www.biorxiv.org/content/10.1101/2020.02.07.938175v1.

**Input data**
All the codes are based on the matrix hdf5 files exported from Hiclib package - binned data - export (https://mirnylab.bitbucket.io/hiclib/binneddata.html#hiclib.binnedData.binnedData.export)

Prior to the matrix export, data processing by the same package (https://mirnylab.bitbucket.io/hiclib/index.html), including iterative mapping of Hi-C/3C-seq raw reads, fragment assignment and filtering, binning and bin-level filtering and iterative bias correction, are recommended.

The convertion of genome-wide interaction matrix to contact probability matrix can be performed by normalizing the sum of entire matrix to the number of rows/columns, so that the sum of each row/column approximates 1 and each value Cij represents the probability of contacts between bin i and bin j.
> matrix = BD.dataDict[NAME]
> totalsum = sum(sum(matrix))
> BD.dataDict[NAME] = matrix / float(totalsum) * len(matrix)

Next, genome-wide and chromosome-wide contact probability matrices can be exported to HDF5 format and used as input for this script.
> BD.export(NAME, GENOME\_HDF5) # export genome-wide matrix
> BD.export(NAME, CHROMOSOME\_HDF5, byChromosome="all") # export chromosome-wide matrices
