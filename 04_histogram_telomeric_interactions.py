#! /usr/bin/python

'''
<Description>
This script generates histogram of non-zero interchromosomal or intrachromosomal long-range interactions. Mean contact probabilities of bulk chromatin as well as intercrhomosomal or intrachromosomal telomeric interactions are calculated and plotted on the histogram. Mann-Whitney U test is performed to compare telomeric interactions and the bulk chromatin.

<Input data>
The hdf5 file of chromosome-wide contact probability matrices.

<Usage>
python 04_04_histogram_telomeric_interactions.py [-h] --chromosomeHDF5
                                              CHROMOSOMEHDF5
                                              [--interChromosomal]
                                              [--intraChromosomal]
                                              [--longRange LONGRANGE]
                                              [--distanceTelomere DISTANCETELOMERE]
                                              [--xmax XMAX] [--ymax YMAX]

<Examples>
1) Plotting histogram of non-zero interchromosomal interactions (with mean value) together with mean contact probability of intercrhomosomal telomeric interactions.
$ python 04_histogram_telomeric_interactions.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --interChromosomal

2) Plotting histogram of non-zero intrachromosomal long-range interactions (with mean value) together with mean contact probability of intracrhomosomal telomeric interactions.
$ python 04_histogram_telomeric_interactions.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --intraChromosomal

3) Changing the definitions of long-range interactions and intrachromosomal telomeric interactions.
$ python 04_histogram_telomeric_interactions.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --intraChromosomal --longRange 200000 --distanceTelomere 20000

4) Customizing figure.
$ python 04_histogram_telomeric_interactions.py --chromosomeHDF5 ../test_mapped_reads_heatmap_10kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --interChromosomal --xmax 0.04 --ymax 70000

<Output>
Histogram of non-zero interchromosomal or intrachromosomal long-range interactions together with mean contact probabilities of bulk chromatin as well as intercrhomosomal or intrachromosomal telomeric interactions in PNG format.
'''

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mirnylib import h5dict
from scipy.stats import mannwhitneyu

# Main function
def main():
	CHROMOSOME_HDF5, INTER_CHROMOSOMAL, INTRA_CHROMOSOMAL, LONG_RANGE, DISTANCE_TELOMERE, XMAX, YMAX, FIGURE_INTER, FIGURE_INTRA = parseInput()
	f, resolution, genomeIdxToLabel = readDict(CHROMOSOME_HDF5)
	if INTER_CHROMOSOMAL:
		trans_tel, trans_all = getTransTelomeric(f, genomeIdxToLabel)
		stats(trans_tel, trans_all, 0)
		plotHistogram(trans_tel, trans_all, FIGURE_INTER, XMAX, YMAX, 0)
	if INTRA_CHROMOSOMAL:
		long_tel, long_all = getCisTelomeric(f, resolution, genomeIdxToLabel, LONG_RANGE, DISTANCE_TELOMERE)
		stats(long_tel, long_all, 0)
		plotHistogram(long_tel, long_all, FIGURE_INTRA, XMAX, YMAX, 1)

# Parse input arguments
def parseInput():
	# Define arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--chromosomeHDF5', dest = 'chromosomeHDF5', required = True, help = 'Input HDF5 file of chromosome-wide contact probability matrices. [Required]')
	parser.add_argument('--interChromosomal', dest = 'interChromosomal', action='store_true', help = 'Plot histogram of interchromosomal telomeric interactions as well as all interchromosomal interactions [Optional]')
	parser.add_argument('--intraChromosomal', dest = 'intraChromosomal', action='store_true', help = 'Plot histogram of intrachromosomal telomeric interactions as well as all intrachromosomal long-range interactions [Optional]')
	parser.add_argument('--longRange', dest = 'longRange', type = int, default = 100000, help = 'The distance between two intrachromosomal loci for defining long-range interactions (e.g. use "--longRange 100000" to define long-range interactions as interactions at distance >= 100 kb). [Optional; Please use when --intraChromosomal is specified; Default = 100000]')
	parser.add_argument('--distanceTelomere', dest = 'distanceTelomere', type = int, default = 10000, help = 'The distance threshold for defining intrachromosomal telomeric interactions (e.g. use "--distanceTelomeric" to define intrachromosomal telomeric interactions as interactions between loci whose distances to two telomeres have a sum of <= 10 kb. In this case, if bin i is 2 kb apart from 5\'-end telomere of one chromosome and bin j is located on the same chromosome but is 4 kb apart from the 3\'-end telomere, then the sum of distances is 2 kb + 4 kb = 6 kb. Since 6 kb <= 10 kb, the interactions between bin i and bin j is considered as intrachromosomal telomeric interaction). [Optional; Please use when --intraChromosomal is specified; Default = 10000]')
	parser.add_argument('--xmax', dest = 'xmax', type = float, default = None, help = 'Maximum value of x-axis. [Optional; Default = None]')
	parser.add_argument('--ymax', dest = 'ymax', type = int, default = None, help = 'Maximum value of y-axis. [Optional; Default = None]')

	# Print help information if no argument is given
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	info = parser.parse_args()
	CHROMOSOME_HDF5 = info.chromosomeHDF5
	INTER_CHROMOSOMAL = info.interChromosomal
	INTRA_CHROMOSOMAL = info.intraChromosomal
	LONG_RANGE = info.longRange
	DISTANCE_TELOMERE = info.distanceTelomere
	XMAX = info.xmax
	YMAX = info.ymax
	
	# Name output figure
	tmp = ''
	if XMAX != None:
		tmp += '_xmax' + str(XMAX)
	if YMAX != None:
		tmp += '_ymax' + str(YMAX)
	FIGURE_INTER = CHROMOSOME_HDF5.split('/')[-1][:-5] + '_histogram_trans_telomeric' + tmp + '.png'
	FIGURE_INTRA = CHROMOSOME_HDF5.split('/')[-1][:-5] + '_histogram_cis_telomeric' + str(DISTANCE_TELOMERE) + '_longRange' + str(LONG_RANGE) + tmp + '.png'

	# Return
	return CHROMOSOME_HDF5, INTER_CHROMOSOMAL, INTRA_CHROMOSOMAL, LONG_RANGE, DISTANCE_TELOMERE, XMAX, YMAX, FIGURE_INTER, FIGURE_INTRA

# Read hdf5 file
def readDict(CHROMOSOME_HDF5):
	# Read file
	f = h5dict.h5dict(CHROMOSOME_HDF5, mode='r')
	resolution = f['resolution']
	genomeIdxToLabel = f['genomeIdxToLabel']
	
	# Return
	return f, resolution, genomeIdxToLabel

# Get interchromosomal telomeric interactions and all interchromosomal interactions
def getTransTelomeric(f, genomeIdxToLabel):
	# Create lists
	trans_tel = []
	trans_all = []
	
	# Add non-zero values in each trans matrix to the lists
	for chrm1 in range(len(genomeIdxToLabel)):
		for chrm2 in range(chrm1 + 1, len(genomeIdxToLabel)):
			key = str(chrm1) + ' ' + str(chrm2)
			m = f[key]
			for i in range(len(m)):
				for j in range(len(m[i])):
					# Skip zero values
	                                if m[i][j] == 0:
        	                                continue
					# Add value to trans_tel list if i and j are from end bins
                        	        if i in [0, len(m) - 1] and j in [0, len(m[i]) - 1]:
                                	        trans_tel.append(m[i][j])
					# Add value to trans_all list
                                	trans_all.append(m[i][j])

	# Return lists
	return trans_tel, trans_all

# Get intrachromosomal telomeric interactions and all intrachromosomal long-range interaction
def getCisTelomeric(f, resolution, genomeIdxToLabel, LONG_RANGE, DISTANCE_TELOMERE):
	# Create lists
	long_tel = []
	long_all = []

	# Convert LONG_RANGE and DISTANCE_TELOMERE to bins
	longRangeBin = LONG_RANGE / resolution
	distanceTelBin = DISTANCE_TELOMERE / resolution

	# Add non-zero values in each cis matrix to the lists
	for chrm in range(len(genomeIdxToLabel)):
		key = str(chrm) + ' ' + str(chrm)
		m = f[key]
		for i in range(len(m)):
			# Take only long-range interactions
			if i >= len(m) - longRangeBin:
				break
			for j in range(i + longRangeBin, len(m)):
				# Skip zero values
                        	if m[i][j] == 0:
                                	continue
				# Add value to long_tel list if the distance between i and 5'-end telomere (d = i) + distance between j and 3'-end telomere (len(m) - 1 - j) <= distanceTelBin
                        	if i + len(m) - 1 - j <= distanceTelBin:
                                	long_tel.append(m[i][j])
				# Add value to long_all list
				long_all.append(m[i][j])
	
	# Return
	return long_tel, long_all

# Statistical analysis for comparing telomeric interactions and bulk chromatin
# Mode == 0: interchromosomal interactions
# Mode == 1: intrachromosomal interactions
def stats(list_tel, list_all, mode):
	if mode == 0:
		print 'Mann-Whitney U test for comparing interchromosomal telomeric interactions and all interchromosomal interactions'
	elif mode == 1:
		print 'Mann-Whitney U test for comparing intrachromosomal telomeric interactions and all intrachromosomal long-range interactions'
	print '[Bulk chromatin] sample size = %s, mean = %s' % (len(list_all), np.mean(list_all))
	print '[Telomeric interactions] sample size = %s, mean = %s' % (len(list_tel), np.mean(list_tel))
	print mannwhitneyu(list_all,list_tel)

# Plot histogram for all interchromosomal interactions or all intrachromosomal long-range interactions as well as the mean contact probabilities of bulk chromatin and telomeric interactions
# Mode == 0: interchromosomal interactions
# Mode == 1: intrachromosomal interactions
def plotHistogram(list_tel, list_all, FIGURE, XMAX, YMAX, mode):
	# Plot histogram of all interchromosomal interactions or all intrachromosomal long-range interactions
	if XMAX != None:
		plt.hist(list_all, bins = 1000, range = (0, XMAX), color = 'grey')
	else:
		plt.hist(list_all, bins = 1000, color = 'grey')

	# Add vertical lines for mean contact probabilities of bulk chromatin and telomeric interactions
	if mode == 0:
		label_all = 'All interchromosomal'
		label_tel = 'Interchromosomal telomeric'
	elif mode == 1:
		label_all = 'All long-range (>100 kb)'
		label_tel = 'Intrachromosomal telomeric'
	plt.axvline(np.mean(list_all), color = 'k', linestyle = '-', lw = 1.5, label = label_all)
	plt.axvline(np.mean(list_tel), color = 'b', linestyle = '-', lw = 1.5, label = label_tel)

	# Add legend
	plt.legend()

	# Add axis titles
	if mode == 0:
		ylabel = 'Counts of non-zero interchromosomal interactions'
	elif mode == 1:
		ylabel = 'Counts of non-zero long-range interactions'
	plt.ylabel(ylabel)
	plt.xlabel('Contact probability')

	# Set axis scales
	if XMAX != None:
		plt.xlim([0, XMAX])
	if YMAX != None:
		plt.ylim([0, YMAX])
	
	# Save figure
	plt.tight_layout()
	plt.savefig(FIGURE, dpi = 500)
	plt.close()

# Execute main function
if __name__ == '__main__':
	main()
