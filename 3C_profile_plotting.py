#! /usr/bin/python

'''
<Description>
This script plots the 3C profile anchored on a specified genomic bin.

<Input data>
The hdf5 file of chromosome-wide contact probability matrices.

<Usage>
python 3C_profile_plotting.py [-h] --chromosomeHDF5 CHROMOSOMEHDF5
                              --anchorChromosome ANCHORCHROMOSOME --anchor
                              ANCHOR [--regionChromosome REGIONCHROMOSOME]
                              [--regionStart REGIONSTART]
                              [--regionEnd REGIONEND] [--ymax YMAX]

<Examples>
1) Plotting the 3C profile of a bin (chr1: 10000-19999) on the entire chromosome
$ python 02_3C_profile_plotting.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --anchorChromosome 1 --anchor 10000

2) Plotting the 3C profile of a bin (chr1: 10000-19999) in a specified region (chr1: 0-99999)
$ python 02_3C_profile_plotting.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNorm
Prob_byChr.hdf5 --anchorChromosome 1 --anchor 10000 --regionStart 0 --regionEnd 99999

3) Plotting the 3C profile of a bin (chr1: 10000-19999) on another chromosome
$ python 02_3C_profile_plotting.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNorm
Prob_byChr.hdf5 --anchorChromosome 1 --anchor 10000 --regionChromosome 2 --regionStart 0 --regionEnd 99999

4) Customizing figure
$ python 02_3C_profile_plotting.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNorm
Prob_byChr.hdf5 --anchorChromosome 1 --anchor 10000 --regionChromosome 2 --regionStart 0 --regionEnd 99999 --ymax 0.01

<Output>
Plot of 3C profile in PNG format.
'''

import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mirnylib import h5dict

# Main function
def main():
	CHROMOSOME_HDF5, ANCHOR_CHROMOSOME, ANCHOR, REGION_CHROMOSOME, REGION_START, REGION_END, YMAX= parseInput()
	matrix, anchorBin, regionStartBin, regionEndBin, resolution, FIGURE = get3CProfile(CHROMOSOME_HDF5, ANCHOR_CHROMOSOME, ANCHOR, REGION_CHROMOSOME, REGION_START, REGION_END)
	plot3CProfile(matrix, ANCHOR_CHROMOSOME, anchorBin, REGION_CHROMOSOME, regionStartBin, regionEndBin, resolution, FIGURE, YMAX)

# Parse input arguments
def parseInput():
	# Define arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--chromosomeHDF5', dest = 'chromosomeHDF5', required = True, help = 'Input HDF5 file of chromosome-wide contact probability matrices. [Required]')
	parser.add_argument('--anchorChromosome', dest = 'anchorChromosome', required = True, help = 'The chromosome of anchor bin (e.g. use "--anchorChromosome 1" to specify the anchor on chromosome 1). [Required]')
	parser.add_argument('--anchor', dest = 'anchor', required = True, type = int, help = 'The coordinate of anchor bin. [Required]')
	parser.add_argument('--regionChromosome', dest = 'regionChromosome', default = None, help = 'The chromosome of region for plotting (e.g. use "--regionChromosome 1" to specify the region on chromosome 1). By default, it plots the chromosome where anchor bin is located. [Optional]')
	parser.add_argument('--regionStart', dest = 'regionStart', default = None, type = int, help = 'The start coordinate of region for plotting. By default, it plots the entire chromosome as specified by --regionChromosome. [Optional]')
	parser.add_argument('--regionEnd', dest = 'regionEnd', default = None, type = int, help = 'The end coordinate of region for plotting. By default, it plots the entire chromosome as specified by --regionChromosome. [Optional]')
	parser.add_argument('--ymax', dest = 'ymax', type = float, default = None, help = 'maximum value of y-axis. [Optional; Default = None]')

	# Print help information if no argument is given
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	info = parser.parse_args()
	CHROMOSOME_HDF5 = info.chromosomeHDF5
	ANCHOR_CHROMOSOME = info.anchorChromosome
	ANCHOR = info.anchor
	REGION_CHROMOSOME = info.regionChromosome
	if REGION_CHROMOSOME == None:
		REGION_CHROMOSOME = ANCHOR_CHROMOSOME
	REGION_START = info.regionStart
	REGION_END = info.regionEnd
	YMAX = info.ymax

	# Return
	return CHROMOSOME_HDF5, ANCHOR_CHROMOSOME, ANCHOR, REGION_CHROMOSOME, REGION_START, REGION_END, YMAX

# Get matrix for plotting 3C profile
def get3CProfile(CHROMOSOME_HDF5, ANCHOR_CHROMOSOME, ANCHOR, REGION_CHROMOSOME, REGION_START, REGION_END):
	# Read hdf5
	f = h5dict.h5dict(CHROMOSOME_HDF5, mode='r')
	genomeIdxToLabel = f['genomeIdxToLabel']
	chromosomeStarts = f['chromosomeStarts']
        binNumber = f['binNumber']

	# Get chromosome information
	for i in range(len(genomeIdxToLabel)):
		chrmLabel = genomeIdxToLabel[i]
		if chrmLabel == ANCHOR_CHROMOSOME:
			anchor_chrmIdx = i
			anchor_chrmLen = _getChrmLen(i, chromosomeStarts, binNumber)
		if chrmLabel == REGION_CHROMOSOME:
			region_chrmIdx = i
			region_chrmLen = _getChrmLen(i, chromosomeStarts, binNumber)
	
	# Convert coordinates to bin numbers
	resolution = f['resolution']
	anchorBin = ANCHOR / resolution
	if REGION_START != None:
		regionStartBin = REGION_START / resolution
	else:
		regionStartBin = 0
	if REGION_END != None:
		regionEndBin = REGION_END / resolution + 1
	else:
		regionEndBin = region_chrmLen

	# Check bin numbers
	if anchorBin > anchor_chrmLen:
		print '[Error] Anchor coordinate (%s) exceeds chromosome length (%s).' % (ANCHOR, anchor_chrmLen * resolution)
		sys.exit(1)
	if regionEndBin < regionStartBin:
		print '[Error] Region start (%s) is larger than region end (%s).' % (regionStartBin * resolution, regionEndBin * resolution)
		sys.exit(1)
	if regionEndBin > region_chrmLen:
		print '[Error] Region (%s-%s) exceed chromosome length (%s).' % (regionStartBin * resolution, regionEndBin * resolution, region_chrmLen * resolution)
		esys.exit(1)

	# Get matrix
	key = str(anchor_chrmIdx) + ' ' + str(region_chrmIdx)
	matrix = f[key]
	matrix = matrix[anchorBin,regionStartBin:regionEndBin]

	# Name output figure
	FIGURE = CHROMOSOME_HDF5.split('/')[-1][:-5] + '_anchor_chr' + ANCHOR_CHROMOSOME + '_' + str(anchorBin * resolution) + '-' + str((anchorBin + 1) * resolution - 1) + '_region_chr' + REGION_CHROMOSOME + '_' + str(regionStartBin * resolution) + '-' + str(regionEndBin * resolution - 1)

	# Return
	return matrix, anchorBin, regionStartBin, regionEndBin, resolution, FIGURE

# Get chromosome length
def _getChrmLen(i, chromosomeStarts, binNumber):
	if i != len(chromosomeStarts) - 1:
		return chromosomeStarts[i+1] - chromosomeStarts[i]
	else:
		return binNumber - chromosomeStarts[i]

# Plot 3C profile
def plot3CProfile(matrix, ANCHOR_CHROMOSOME, anchorBin, REGION_CHROMOSOME, regionStartBin, regionEndBin, resolution, FIGURE, YMAX):
	# Plot anchor
	anchor_label = 'Anchor bin (chr' + ANCHOR_CHROMOSOME + ': ' + str(anchorBin * resolution) + '-' + str((anchorBin + 1) * resolution - 1) + ')'
	if ANCHOR_CHROMOSOME == REGION_CHROMOSOME:
		plt.axvspan(anchorBin * resolution, (anchorBin + 1) * resolution, facecolor = 'grey', alpha = 0.8, label = anchor_label)
	else:
		plt.axvspan(None, None, facecolor = 'none', label = anchor_label)
	plt.legend(bbox_to_anchor = (1, 1), loc = 'lower right', fontsize = 'small')
	
	# Plot 3C profile
	x = range(regionStartBin * resolution, regionEndBin * resolution, resolution)
	plt.scatter(x, matrix, s = 0.5)

	# Set axis scales
	plt.xlim((regionStartBin - 1) * resolution, regionEndBin * resolution)
	if YMAX != None:
		plt.ylim(0, YMAX)
		FIGURE += '_ymax' + str(YMAX)
	FIGURE += '.png'

	# Plot x-axis and y-axis titles
	plt.xlabel('Coordinates on chr%s' % REGION_CHROMOSOME)
	plt.ylabel('Contact probability')

	# Save figure
	plt.savefig(FIGURE, dpi = 500)
	plt.close()

#Execute main function
if __name__ == '__main__':
	main()
