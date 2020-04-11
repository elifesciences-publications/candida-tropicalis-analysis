#! /usr/bin/python

"""
<Description>
This script generates a Hi-C/3C-seq contact heatmap from the genome-wide/chromosome-wide contact probability matrices.

<Input data>
The hdf5 file of genome-wide and/or chromosome-wide contact probability matrices.

<Usage>
python heatmap_generation.py [-h] [--genomeHDF5 GENOMEHDF5]
                             [--chromosomeHDF5 CHROMOSOMEHDF5] [--plotGenome]
                             [--plotChromosome PLOTCHROMOSOME [PLOTCHROMOSOME ...]]
                             [--colorMap COLORMAP] [--vmin VMIN] [--vmax VMAX]
                             [--chromosomeGrid CHROMOSOMEGRID]

<Examples>
1) Plotting a genome-wide heatmap
$ python 01_heatmap_generation.py --genomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb.hdf5 --plotGenome

2) Plotting a chromosome-wide heatmap
$ python 01_heatmap_generation.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotChromosome 1

3) Plotting a heatmap of more than one chromosomes
$ python 01_heatmap_generation.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotChromosome 1 2

4) Plotting both genome-wide and chromosome-wide heatmaps at the same time
$ python 01_heatmap_generation.py --genomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb.hdf5 --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotGenome --plotChromosome 1 2

5) Customizing heatmap
$ python 01_heatmap_generation.py --genomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb.hdf5 --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotGenome --plotChromosome 1 2 --colorMap Blues --vmin -12 --vmax -5

<Output>
Genome-wide/chromosome-wide contact probability heatmaps in PNG format.
"""

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mirnylib import h5dict, plotting

# Main function
def main():
	GENOME_HDF5, CHROMOSOME_HDF5, PLOT_GENOME, PLOT_CHROMOSOME, COLORMAP, VMIN, VMAX, CHROMOSOME_GRID, HEATMAP_GENOME, HEATMAP_CHROMOSOMES = parseInput()
	if PLOT_GENOME:
		matrix, resolution, grids, chrmLabels = getGenomeMatrix(GENOME_HDF5)
		plotHeatmap(matrix, resolution, grids, chrmLabels, HEATMAP_GENOME, COLORMAP, VMIN, VMAX, CHROMOSOME_GRID)
	if PLOT_CHROMOSOME:
		matrix, resolution, grids, chrmLabels = getChromosomesMatrix(CHROMOSOME_HDF5, PLOT_CHROMOSOME)
		plotHeatmap(matrix, resolution, grids, chrmLabels, HEATMAP_CHROMOSOMES, COLORMAP, VMIN, VMAX, CHROMOSOME_GRID)

# Parse input arguments
def parseInput():
	# Define arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--genomeHDF5', dest = 'genomeHDF5', default = None, help = 'Input HDF5 file of Genome-wide contact probability matrix. [Required when --plotGenome is speficied]')
	parser.add_argument('--chromosomeHDF5', dest = 'chromosomeHDF5', default = None, help = 'Input HDF5 file of chromosome-wide contact probability matrices. [Required when --plotChromosome is specified]')
	parser.add_argument('--plotGenome', dest = 'plotGenome', action='store_true', help = 'Plot genome-wide contact probability heatmap. [Optional]')
	parser.add_argument('--plotChromosome', dest = 'plotChromosome', nargs='+', default = 0, help = 'Plot heatmap of one or more chromosomes. For single chromosome-wide heatmap, chromosome number is required (e.g. use "--plotChromosome 1"). For heatmap of multiple chromosomes, list of chromosomes is required where chromosome numbers are separared by space (e.g. use "--plotChromosome 1 2" to plot a heatmap of chromosome 1 and 2). [Optional]')
	parser.add_argument('--colorMap', dest = 'colorMap', default = 'hot_r', help = 'Built-in colormaps in Matplotlib package (https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html) [Optional; Default = hot_r]')
	parser.add_argument('--vmin', dest = 'vmin', type = float, default = None, help = 'minimum value of color scale. [Optional; Default = None]')
	parser.add_argument('--vmax', dest = 'vmax', type = float, default = None, help = 'maximum value of color scale. [Optional; Default = None]')
	parser.add_argument('--chromosomeGrid', dest = 'chromosomeGrid', type = bool, default = True, help = 'Add grid lines to start and end of each chromosome. [Optional; Default = True]')

	# Print help information if no argument is given
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	info = parser.parse_args()
	GENOME_HDF5 = info.genomeHDF5
	CHROMOSOME_HDF5 = info.chromosomeHDF5
	PLOT_GENOME = info.plotGenome
	PLOT_CHROMOSOME = info.plotChromosome
	COLORMAP = info.colorMap
	VMIN = info.vmin
	VMAX = info.vmax
	CHROMOSOME_GRID = info.chromosomeGrid

	# Check arguments
	if GENOME_HDF5 == None and  PLOT_GENOME:
		print "[Error] --genomeHDF5 is required when --plotGenome is specified."
		sys.exit(1)
	if CHROMOSOME_HDF5 == None and PLOT_CHROMOSOME:
		print "[Error] --chromosomeHDF5 is required when --plotChromosome is specified."
		sys.exit(1)

	# Name output figure
	tmp = "_" + COLORMAP
	if VMIN != None:
		tmp += "_vmin" + str(VMIN)
	if VMAX != None:
		tmp += "_vmax" + str(VMAX)
	if CHROMOSOME_GRID:
		tmp += "_grid"
	HEATMAP_GENOME = GENOME_HDF5.split("/")[-1][:-5] + tmp + ".png"
	HEATMAP_CHROMOSOMES = ""
	if PLOT_CHROMOSOME:
		HEATMAP_CHROMOSOMES = CHROMOSOME_HDF5.split("/")[-1][:-5] + "_chr" + "_chr".join(PLOT_CHROMOSOME) + tmp + ".png"

	# Return
	return GENOME_HDF5, CHROMOSOME_HDF5, PLOT_GENOME, PLOT_CHROMOSOME, COLORMAP, VMIN, VMAX, CHROMOSOME_GRID, HEATMAP_GENOME, HEATMAP_CHROMOSOMES

# Get genome-wide matrix
def getGenomeMatrix(GENOME_HDF5):
	# Read hdf5
	f = h5dict.h5dict(GENOME_HDF5, mode='r')
	matrix = f["heatmap"]
	resolution = f["resolution"]
	chromosomeStarts = f["chromosomeStarts"]
	binNumber = f["binNumber"]

	# Get grid location and chromosome labels
	grids = list(chromosomeStarts) + [binNumber]
	genomeIdxToLabel = f['genomeIdxToLabel']
	chrmLabels = genomeIdxToLabel.values()

	# Return
	return matrix, resolution, grids, chrmLabels

# Get chromosome-wide matrix
def getChromosomesMatrix(CHROMOSOME_HDF5, PLOT_CHROMOSOME):
	# Read hdf5
        f = h5dict.h5dict(CHROMOSOME_HDF5, mode="r")
	resolution = f["resolution"]
        genomeIdxToLabel = f["genomeIdxToLabel"]
        chromosomeStarts = f["chromosomeStarts"]
        binNumber = f["binNumber"]

	# Get grid location and chromosome labels
        chrmIdx = []
        grids = [0]
	chrmLabels = []
        for i in range(len(genomeIdxToLabel)):
                if genomeIdxToLabel[i] in PLOT_CHROMOSOME:
                        chrmIdx.append(i)
			chrmLabels.append(genomeIdxToLabel[i])
                        if i == len(genomeIdxToLabel) - 1:
                                size = binNumber - chromosomeStarts[i]
                        else:
                                size = chromosomeStarts[i+1] - chromosomeStarts[i]
                        grids.append(grids[-1] + size)
	
	# Get chromosome-wide matrix
        slices = []
        for i in chrmIdx:
                m = []
                for j in chrmIdx:
                        key = str(i) + " " + str(j)
                        m.append(f[key])
                slices.append(np.concatenate(m,axis=1))
        matrix = np.concatenate(slices,axis=0)

	# Return
	return matrix, resolution, grids, chrmLabels

# Plot heatmap from matrix
def plotHeatmap(matrix, resolution, grids, chrmLabels, HEATMAP, COLORMAP, VMIN, VMAX, CHROMOSOME_GRID):
	# Plot heatmap
	plotting.plot_matrix(np.log2(matrix),cmap=COLORMAP,vmin=VMIN,vmax=VMAX,label="Contact probability (log2)")

        if len(chrmLabels) > 1 and CHROMOSOME_GRID:
		# Plot grid for heatmap of more than one chromosomes
              	xticks = []
                xticklabels = []
                yticks = []
                yticklabels = []
                for i in range(len(chrmLabels)):
                        start = grids[i]
			end = grids[i+1]
			for j in range(len(chrmLabels)):
				start2 = grids[j]
				end2 = grids[j+1]
	                        plt.axhline(start-0.5,xmin=(start2+15.0)/grids[-1],xmax=(end2-15.0)/grids[-1],color="k",linewidth=1)
        	                plt.axvline(start-0.5,ymin=1-(start2+15.0)/grids[-1],ymax=1-(end2-15.0)/grids[-1],color="k",linewidth=1)
			plt.axhline(grids[-1]-0.5,xmin=(start+15.0)/grids[-1],xmax=(end-15.0)/grids[-1],color="k",linewidth=1)
                        plt.axvline(grids[-1]-0.5,ymin=1-(start+15.0)/grids[-1],ymax=1-(end-15.0)/grids[-1],color="k",linewidth=1)
                        xticks.append((start-0.5+end-0.5)/2)
                        yticks.append((start-0.5+end-0.5)/2)
		
		# Set tick labels to chromosomes
                plt.xticks(xticks,chrmLabels)
                plt.yticks(yticks,chrmLabels)

		# Set x-axis and y-axis titles
		plt.xlabel("Chromosome")
		plt.ylabel("Chromosome")
        else:
		# Set x-axis and y-axis titles for heatmap with no grid
                plt.xlabel("Genomic bin number")
                plt.ylabel("Genomic bin number")
	
	# Set label location
	plt.gca().xaxis.set_tick_params(label1On=False,label2On=True,tick1On=False,tick2On=True)
        plt.gca().yaxis.set_tick_params(label1On=True,label2On=False,tick1On=True,tick2On=False)
	plt.gca().xaxis.set_label_position("top")

	# Remove figure frame
	for spine in plt.gca().spines.values():
                spine.set_visible(False)
	
	# Save figure
	plt.savefig(HEATMAP,dpi=500)
	plt.close()

# Execute main function
if __name__ == "__main__":
	main()
