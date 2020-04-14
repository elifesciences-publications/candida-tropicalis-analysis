#! /usr/bin/python

'''
<Description>
This script performs aggregate signal analysis of contact probability matrix and generates the heatmap of mean contact probabilities surrounding all centromere-centromere or telomere-telomere interactions.

<Input data>
The hdf5 file of chromosome-wide contact probability matrices.

<Usage>
python aggregate_heatmap_at_centromere_telomere.py [-h] --chromosomeHDF5
                                                   CHROMOSOMEHDF5
                                                   [--plotTelomere PLOTTELOMERE]
                                                   [--plotCentromere PLOTCENTROMERE]
                                                   [--centromereBed CENTROMEREBED]
                                                   [--colorMap COLORMAP]
                                                   [--vmin VMIN] [--vmax VMAX]

<Examples>
1) Plotting heatmap of mean contact probabilities (bin size = 2 kb) in 50 kb upstream and 50 kb downsteam regions surrounding centromere-centromere interactions (50 kb upstream + 2 kb centromere mid bin + 50 kb downstream = total 102 kb in 51 bins)
$ python aggregate_heatmap_at_centromere_telomere.py --chromosomeHDF5 ../test_mapped_reads_heatmap_2kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotCentromere 50000 --centromereBed centromeres.bed

For Candida tropicalis, the centromeres.bed includes:
chr1	466066	475754
chr2	725185	735300
chr3	952377	962681
chr4	906302	916031
chr5	594006	604320
chr6	308570	331012
chrR	855697	865834

2) Plotting heatmap of mean contact probabilities (bin size = 2 kb) in 100 kb region adjacent to telomere-telomere interactions (100 kb upstream/downstream + 2 kb telomere end bin = 102 kb in 51 bins)
$ python aggregate_heatmap_at_centromere_telomere.py --chromosomeHDF5 ../test_mapped_reads_heatmap_2kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotTelomere 100000

3) Customizing output heatmap
$ python aggregate_heatmap_at_centromere_telomere.py --chromosomeHDF5 ../test_mapped_reads_heatmap_2kb_N0.5_noTruncTrans_IC_readNormProb_byChr.hdf5 --plotCentromere 50000 --centromereBed centromeres.bed --colorMap Blues --vmin -12 --vmax -7

<Output>
Heatmap of mean contact probabilities surrounding all centromere-centromere or telomere-telomere interactions in PNG format.
'''

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mirnylib import h5dict

# Main function
def main():
	CHROMOSOME_HDF5, PLOT_TELOMERE, PLOT_CENTROMERE, CENTROMERE_BED, COLORMAP, VMIN, VMAX, FIGURE_CENTROMERE, FIGURE_TELOMERE = parseInput()
	f, resolution, genomeIdxToLabel = readDict(CHROMOSOME_HDF5)
	if PLOT_CENTROMERE:
		agg = aggregateSignalCentromere(f, resolution, genomeIdxToLabel, PLOT_CENTROMERE, CENTROMERE_BED)
		plotHeatmap(agg, FIGURE_CENTROMERE, COLORMAP, VMIN, VMAX, 0)
	if PLOT_TELOMERE:
		agg = aggregateSignalTelomere(f, resolution, genomeIdxToLabel, PLOT_TELOMERE)
		plotHeatmap(agg, FIGURE_TELOMERE, COLORMAP, VMIN, VMAX, 1)

# Parse input arguments
def parseInput():
	# Define arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--chromosomeHDF5', dest = 'chromosomeHDF5', required = True, help = 'Input HDF5 file of chromosome-wide contact probability matrices. [Required]')
	parser.add_argument('--plotTelomere', dest = 'plotTelomere', type = int, default = 0, help = 'Number of base pairs plotted from the chromosome end (telomere). [Optional]')
	parser.add_argument('--plotCentromere', dest = 'plotCentromere', type = int, default = 0, help = 'Number of base pairs plotted from each side of centromere mid bin (e.g. use "--plotCentromere 50000" to plot a region containing 50000 bp upstream region, centromere mid bin and 50000 bp downstream region). [Optional]')
	parser.add_argument('--centromereBed', dest = 'centromereBed', default = None, help = 'BED file containing start and end coordinates of centromere for each chromosome. Format of columns (separated by tabs): Chromosome, Start, End. [Required if --plotCentromere is specified]')
	parser.add_argument('--colorMap', dest = 'colorMap', default = 'hot_r', help = 'Built-in colormaps in Matplotlib package (https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html) [Optional; Default = hot_r]')
	parser.add_argument('--vmin', dest = 'vmin', type = float, default = None, help = 'minimum value of color scale. [Optional; Default = None]')
	parser.add_argument('--vmax', dest = 'vmax', type = float, default = None, help = 'maximum value of color scale. [Optional; Default = None]')

	# Print help information if no argument is given
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	info = parser.parse_args()
	CHROMOSOME_HDF5 = info.chromosomeHDF5
	PLOT_TELOMERE = info.plotTelomere
	PLOT_CENTROMERE = info.plotCentromere
	CENTROMERE_BED = info.centromereBed
	COLORMAP = info.colorMap
	VMIN = info.vmin
	VMAX = info.vmax
	
	# Check arguments
	if CENTROMERE_BED == None and PLOT_CENTROMERE:
		print '[Error] --centromereBed is required when --plotCentromere is specified.'
		sys.exit(1)

	# Name output figure
	tmp = '_' + COLORMAP
	if VMIN != None:
		tmp += '_vmin' + str(VMIN)
	if VMAX != None:
		tmp += '_vmax' + str(VMAX)
	FIGURE_CENTROMERE = CHROMOSOME_HDF5.split('/')[-1][:-5] + '_aggregate_signal_centromere' + str(PLOT_CENTROMERE*2+1) + tmp + '.png'
	FIGURE_TELOMERE = CHROMOSOME_HDF5.split('/')[-1][:-5] + '_aggregate_signal_telomere' + str(PLOT_TELOMERE) + tmp + '.png'

	# Return
	return CHROMOSOME_HDF5, PLOT_TELOMERE, PLOT_CENTROMERE, CENTROMERE_BED, COLORMAP, VMIN, VMAX, FIGURE_CENTROMERE, FIGURE_TELOMERE

# Read hdf5 file
def readDict(CHROMOSOME_HDF5):
	# Read file
	f = h5dict.h5dict(CHROMOSOME_HDF5, mode='r')
	resolution = f['resolution']
	genomeIdxToLabel = f['genomeIdxToLabel']
	
	# Return
	return f, resolution, genomeIdxToLabel

# Aggregate signal at centromere
def aggregateSignalCentromere(f, resolution, genomeIdxToLabel, PLOT_CENTROMERE, CENTROMERE_BED):
	# Get centromere center locations
	cens = _getCentromereCenter(CENTROMERE_BED, resolution)

	# Aggregate signal analysis
	N = PLOT_CENTROMERE / resolution
	agg = np.zeros((N * 2 + 1, N * 2 + 1))
	count = 0
	for chrm1 in range(len(genomeIdxToLabel) - 1):
        	for chrm2 in range(chrm1 + 1, len(genomeIdxToLabel)):
                	cen1 = cens[genomeIdxToLabel[chrm1]]
                	cen2 = cens[genomeIdxToLabel[chrm2]]
                	key = str(chrm1) + ' ' + str(chrm2)
                	agg += f[key][cen1 - N : cen1 + 1 + N, cen2 - N : cen2 + 1 + N]
                	count += 1
	agg /= count

	# Return aggregated signal
	return agg

# Get centromere center locations
def _getCentromereCenter(CENTROMERE_BED, resolution):
	# Create dictionary
	cens = {}

	# Read BED file of centromere locations
	with open(CENTROMERE_BED, 'r') as f:
	        for l in f:
        	        if not l.strip().startswith('#'):
                	        tmp = l.strip().split('\t')
                        	chrmLabel = tmp[0][3:]
                        	center = (int(tmp[1]) + int(tmp[2])) / 2
                        	cens[chrmLabel] = center / resolution
	
	# Return dictionary
	return cens

def aggregateSignalTelomere(f, resolution, genomeIdxToLabel, PLOT_TELOMERE):
	# Aggregated signal analysis
	N = PLOT_TELOMERE / resolution + 1
	agg = np.zeros((N, N))
	count = 0
	for l in genomeIdxToLabel:
        	for r in genomeIdxToLabel:
        	        key = str(l) + ' ' + str(r)
                	agg += f[key][:N, -N:]
                	count += 1
	agg /= count

	# Return aggregated signal
	return agg

# Plot aggregated signal
# Mode == 0: plot for centromere-centromere interactions
# Mode == 1: plot for telomere-telomere interactions
def plotHeatmap(agg, FIGURE, COLORMAP, VMIN, VMAX, mode):
	# Plot heatmap
	im = plt.imshow(np.log2(agg), cmap = COLORMAP, vmin = VMIN, vmax = VMAX)
	plt.colorbar(im, label = 'Contact probability (log2)')

	# Set tick labels
	if mode == 0:
		ticks = [len(agg) / 2]
		ticklabels = ['Centromere']
		plt.xticks(ticks, ticklabels, ha = 'center')
		plt.yticks(ticks, ticklabels, rotation = 90, va = 'center')
	elif mode == 1:
		ticklabels = ['Telomere']
		plt.xticks([len(agg) - 1], ticklabels, ha = 'right')
		plt.yticks([0], ticklabels, rotation = 90, va = 'top')

	# Set label location
	plt.gca().xaxis.set_tick_params(label1On = False, label2On = True, tick1On = False, tick2On = True)
        plt.gca().yaxis.set_tick_params(label1On = True, label2On = False, tick1On = True, tick2On = False)

	# Save figure
	plt.savefig(FIGURE,dpi=500)
	plt.close()

# Execute main function
if __name__ == '__main__':
	main()
