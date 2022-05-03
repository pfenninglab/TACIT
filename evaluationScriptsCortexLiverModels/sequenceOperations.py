import sys
import os
import subprocess
import numpy as np
import pybedtools as bt
from itertools import izip
from Bio import SeqIO
import pyximport
importers = pyximport.install()
pyximport.uninstall(*importers)

class narrowPeak(object):
	# narrowPeak class
	def __init__(self, peakLine):
		peakLineElements = peakLine.strip().split("\t")
		self.chrom = peakLineElements[0]
		self.start = int(peakLineElements[1])
		self.end = int(peakLineElements[2])
		if len(peakLineElements) > 3:
			# The line is more than a 3-column bed
			self.name = peakLineElements[3]
			self.browserScore = int(peakLineElements[4])
			self.strand = peakLineElements[5]
			if len(peakLineElements) > 6:
				# The line is more than a 6-column bed
				self.signal = float(peakLineElements[6])
				self.pVal = float(peakLineElements[7])
				self.qVal = float(peakLineElements[8])
				self.summitOffset = int(peakLineElements[9])
			else:
				# The line is only a 6-column bed
				self.signal = -1
				self.pVal = -1
				self.qVal = -1
				self.summitOffset = -1
		else:
			# The line is only a 3-column bed
			self.name = "."
			self.browserScore = 0
			self.strand = "."
			self.signal = -1
			self.pVal = -1
			self.qVal = -1
			self.summitOffset = -1
	def getPeakWidth(self):
		# Get the width of the peak
		return self.end - self.start
	def beyondChromosome(self, chromSizesFileName):
		# Determine if the peak is beyond the edge of a chromosome
		if self.start < 0:
			# The peak starts at a negative coordinate
			return True
		chromSizesFile = open(chromSizesFileName)
		for line in chromSizesFile:
			# Iterate through the files with the chromosome sizes and check if the peak ends after the end of a chromosome
			lineElements = line.strip().split("\t")
			if lineElements[0] == self.chrom:
				# At the chromosome of the narrowPeak
				if self.end >= int(lineElements[1]):
					# The current peak goes beyond the end of the chromosome, so do not use it
					return True
				else:
					return False
		print("Problem: Chromosome " + self.chrom + " not in list!")
		return False
	def getPeakBins(self, binSize, sequenceLength, binStride, chromSizesFileName):
		# Get the bins for the peak and create a new peak for each bin
		peaksFromBins = []
		for binStart in xrange(self.start, self.end - (binSize/2), binStride):
			# Iterate through the bins and create a new peak for each bin
			binEnd = binStart + binSize
			peakFromBin =\
				narrowPeak("\t".join([self.chrom, str(binStart - ((sequenceLength - binSize) / 2)), \
					str(binEnd + ((sequenceLength - binSize) / 2)), self.name, str(self.browserScore), self.strand, str(self.signal), \
					str(self.pVal), str(self.qVal), str(sequenceLength / 2)]))
			if peakFromBin.beyondChromosome(chromSizesFileName):
				# The peak that has been created from the current bin is too long
				continue
			peaksFromBins.append(peakFromBin)
		return peaksFromBins
	def overlappingPeak(self, peak):
		# Check if 2 narrow peaks have overlapping coordinates
		if self.chrom == peak.chrom:
			# The peaks are on the same chromosome
			if self.start <= peak.start:
				# The narrowPeak starts before peak
				if peak.start < self.end:
					# peak starts within narrowPeak
					return True
				else:
					# peak is after narrowPeak
					return False
			elif self.start < peak.end:
				# narrowPeak starts within peak
				return True
			else:
				# narrowPeak is after peak
				return False
		else:
			# The peaks are on different chromosomes
			return False
	def mergePeak(self, peak):
		# Merge 2 narrow peaks
		assert(self.chrom == peak.chrom)
		mergedPeak = None
		if self.start <= peak.start:
			# The start of the merged peak will be the start of narrowPeak
			if self.end >= peak.end:
				# The end of the merged peak will be the end of narrowPeak
				mergedPeak =\
					narrowPeak("\t".join([self.chrom, str(self.start), str(self.end), ".", "0", ".", str((self.signal + peak.signal)/2), str((self.pVal + peak.pVal)/2), \
						str((self.qVal + peak.qVal)/2), str((self.summitOffset + peak.summitOffset)/2)]))
			else:
				# The end of the merged peak will the end of peak
				mergedPeak =\
					narrowPeak("\t".join([self.chrom, str(self.start), str(peak.end), ".", "0", ".", str((self.signal + peak.signal)/2), str((self.pVal + peak.pVal)/2), \
						str((self.qVal + peak.qVal)/2), str((self.summitOffset + peak.summitOffset)/2)]))
		else:
			# The start of the merged peak will be the start of peak
			if self.end >= peak.end:
				# The end of the merged peak will be the end of narrowPeak
				mergedPeak =\
					narrowPeak("\t".join([self.chrom, str(peak.start), str(self.end), ".", "0", ".", str((self.signal + peak.signal)/2), str((self.pVal + peak.pVal)/2), \
						str((self.qVal + peak.qVal)/2), str((self.summitOffset + peak.summitOffset)/2)]))
			else:
				# The end of the merged peak will the end of peak
				mergedPeak =\
					narrowPeak("\t".join([self.chrom, str(peak.start), str(peak.end), ".", "0", ".", str((self.signal + peak.signal)/2), str((self.pVal + peak.pVal)/2), \
						str((self.qVal + peak.qVal)/2), str((self.summitOffset + peak.summitOffset)/2)]))
		return mergedPeak, peak.signal
		
class peaksInfo(object):
	def __init__(self, narrowPeakList, labels, signals):
		self.peaks = narrowPeakList
		self.labels = labels
		self.scores = signals
		if len(labels.shape) == 1:
			# Make the labels and scores 2D
			self.labels = labels[:,np.newaxis]
			self.scores = signals[:,np.newaxis]
	def removeFullyAmbiguousEntries(self):
		# Remove all entries for which all of the labels are -1
		labelsSum = np.sum(self.labels, axis = 1)
		fullyAmbiguousSum = -1 * (self.labels.shape[1])
		notFullyAmbiguousIndices = np.where(labelsSum > fullyAmbiguousSum)[0]
		peaksNotFullyAmbiguous = []
		for i in range(len(self.peaks)):
			# Iterate through the peaks and save those that do not have exclusively ambiguous labels
			if i in notFullyAmbiguousIndices:
				peaksNotFullyAmbiguous.append(self.peaks[i])
		labelsNotFullyAmbiguous = self.labels[notFullyAmbiguousIndices]
		scoresNotFullyAmbiguous = self.scores[notFullyAmbiguousIndices]
		peaksInfoNotFullyAmbiguous = peaksInfo(peaksNotFullyAmbiguous, labelsNotFullyAmbiguous, scoresNotFullyAmbiguous)
		return peaksInfoNotFullyAmbiguous
	def subsetDataByChrom(self, chromList):
		# Select the entries whose chromosomes are in the list of chromosomes
		subsetPeaksList = []
		subsetIndexList = []
		for i in range(len(self.peaks)):
			# Iterate through the peaks and get those that are on the list of chromosomes
			if self.peaks[i].chrom in chromList:
				# Add the current peak to the list
				subsetPeaksList.append(self.peaks[i])
				subsetIndexList.append(i)
		subsetIndices = np.array(subsetIndexList)
		subsetLabels = self.labels[subsetIndices]
		subsetScores = self.scores[subsetIndices]
		subsetPeaksInfo = peaksInfo(subsetPeaksList, subsetLabels, subsetScores)
		return subsetPeaksInfo

def show_value(s):
	"""
	Convert unicode to str under Python 2;
	all other values pass through unchanged
	"""
	if sys.version_info.major == 2:
		if isinstance(s, unicode):
			return str(s)
	return s

def oneHotEncode(sequence):
	# Also in sequenceOperationsCore.py
	encodedSequence = np.zeros((4, len(sequence)), dtype=np.int8)
	sequenceDict = {}
	sequenceDict["A"] = np.array([1, 0, 0, 0])
	sequenceDict["a"] = np.array([1, 0, 0, 0])
	sequenceDict["C"] = np.array([0, 1, 0, 0])
	sequenceDict["c"] = np.array([0, 1, 0, 0])
	sequenceDict["G"] = np.array([0, 0, 1, 0])
	sequenceDict["g"] = np.array([0, 0, 1, 0])
	sequenceDict["T"] = np.array([0, 0, 0, 1])
	sequenceDict["t"] = np.array([0, 0, 0, 1])
	sequenceDict["N"] = np.array([0, 0, 0, 0])
	sequenceDict["n"] = np.array([0, 0, 0, 0])
	# These are all 0's even though they should ideally have 2 indices with 0.5's because storing ints requires less space than storing floats
	sequenceDict["R"] = np.array([0, 0, 0, 0])
	sequenceDict["r"] = np.array([0, 0, 0, 0])
	sequenceDict["Y"] = np.array([0, 0, 0, 0])
	sequenceDict["y"] = np.array([0, 0, 0, 0])
	sequenceDict["M"] = np.array([0, 0, 0, 0])
	sequenceDict["m"] = np.array([0, 0, 0, 0])
	sequenceDict["K"] = np.array([0, 0, 0, 0])
	sequenceDict["k"] = np.array([0, 0, 0, 0])
	sequenceDict["W"] = np.array([0, 0, 0, 0])
	sequenceDict["w"] = np.array([0, 0, 0, 0])
	sequenceDict["S"] = np.array([0, 0, 0, 0])
	sequenceDict["s"] = np.array([0, 0, 0, 0])
	for i in range(len(sequence)):
		# Iterate through the bases in the sequence and record each
		encodedSequence[:,i] = sequenceDict[sequence[i]]
	numNs = len(sequence) - np.sum(encodedSequence)
	return encodedSequence, numNs
		
def reverse_complement(encoded_sequences):
	# Because the encoding is A, C, G, T in that order, can just reverse each sequence along both axes.
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	return encoded_sequences[..., ::-1, ::-1]

def loadPerBaseTracks(perBaseTrackFileNames):
	# Load the per base tracks
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	perBaseTracks = []
	if perBaseTrackFileNames:
		# Load the per base track files
		perBaseTracks = [np.loadtxt(pbtfn) for pbtfn in perBaseTrackFileNames]
	return perBaseTracks
	
def createPerBaseTracksMat(perBaseTracks, width, sampleCount, divisor):
	# Create a matrix with the per base tracks for the current sample
	# ASSUMES THAT sampleCount IS A MULTIPLE OF divisor
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	perBaseTracksIndex = sampleCount / divisor
	perBaseTracksMat = np.empty((0, width))
	for pbt in perBaseTracks:
		# Iterate through the per base 
		perBaseTracksMat = np.vstack((perBaseTracksMat, pbt[perBaseTracksIndex, :]))
	return perBaseTracksMat
	
def makeMultiModedData(allData, dataShape, numPerBaseTracks):
	# Convert data into the format for multi-moding
	# ASSUMES per-base tracks are 1 high
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	assert(numPerBaseTracks > 0)
	allDataList = []
	allDataList.append(allData[:,:,0:dataShape[1],:])
	for i in range(numPerBaseTracks):
		# Iterate through the per-base tracks and add the data for that track to the list
		allDataList.append(allData[:,:,dataShape[1] + i - 1:dataShape[1] + i,:])
	return allDataList

def convertFastaFileToSequencesFile(fastaFileName, padUpstream=0, padDownstream=0):
	# Convert a fasta file to a sequences file
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	numSequences = 0
	sequencesFileName = ".".join(fastaFileName.split(".")[-0:-1]) + "_sequences.txt"
	sequencesFile = open(sequencesFileName, 'w+')
	sequenceIDs = []
	for record in SeqIO.parse(fastaFileName, "fasta"):
		sequence = str(record.seq)
		for i in range(padUpstream):
			# Add an N to the sequence upstream
			sequence = "N" + sequence
		for i in range(padDownstream):
			# Add an N to the sequence downstream
			sequence = sequence + "N"
		sequencesFile.write(sequence + "\n")
		numSequences = numSequences + 1
		sequenceIDs.append(record.id.strip())
	sequencesFile.close()
	return sequencesFileName, numSequences, sequenceIDs
		
def makeSequenceInputArraysNoLabels(sequenceFileName, dataShape, numSequences, perBaseTrackFileNames=[], multiMode=False, maxFracNs = 1.0):
	# Convert each sequence into a numpy array, but do not load any labels/signals files
	# ASSUMES THAT THE SEQUENCES ARE LISTS AND NOT IN FASTA FORMAT
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	sequenceFile = open(sequenceFileName)
	perBaseTracks = loadPerBaseTracks(perBaseTrackFileNames)
	channel1 = dataShape[0];
	channel2 = dataShape[1] + len(perBaseTracks);
	channel3 = dataShape[2];
	allData = np.zeros((numSequences*2, channel1, channel2, channel3), dtype=np.int8);
	if perBaseTracks:
		# There are additional per-base tracks that might not be ints
		allData = np.zeros((numSequences*2, channel1, channel2, channel3), dtype=np.float16);
	sampleCount = 0
	skippedIndices = []
	totalNs = 0
	index = 0
	for sequence in sequenceFile:
		# Iterate through the fasta sequences and create the alphabets for the sequence and the reverse complement of each
		perBaseTracksMat = createPerBaseTracksMat(perBaseTracks, channel3, sampleCount, 2)
		sequenceArray, numNs = oneHotEncode(sequence.strip())
		sequenceFracNs = float(numNs)/float(dataShape[2])
		index = index + 1
		if sequenceFracNs > maxFracNs:
			# The percentage of N's in the current sequence is too high
			numSequences = numSequences - 1
			skippedIndices.append(index - 1)
			continue
		if sequenceArray.shape[1] != dataShape[2]:
			# The current sequences is the wrong length, so skip it
			print("This sequence is the wrong length: " + sequence)
			skippedIndices.append(index - 1)
			numSequences = numSequences - 1
			continue
		totalNs = totalNs + numNs
		sequenceArrayReshape = np.reshape(np.vstack((sequenceArray, perBaseTracksMat)), (channel1, channel2, channel3))
		allData[sampleCount,:,:,:] = sequenceArrayReshape
		sampleCount = sampleCount + 1
		# Repeat for the reverse complement
		sequenceArrayRC = reverse_complement(sequenceArray)
		sequenceArrayReshapeRC = np.reshape(np.vstack((sequenceArrayRC, perBaseTracksMat)), (channel1, channel2, channel3))
		allData[sampleCount,:,:,:] = sequenceArrayReshapeRC
		sampleCount = sampleCount + 1
	assert (sampleCount == numSequences*2)
        allData = allData[0:sampleCount,:,:,:] # Remove empty entries for skipped indices
	fracNs = float(totalNs)/float(dataShape[2] * numSequences)
	print("The fraction of Ns is: " + str(fracNs))
	sequenceFile.close()
        if multiMode:
		# Re-format the data for multi-moding
		return makeMultiModedData(allData, dataShape, len(perBaseTrackFileNames))
	return allData, skippedIndices
	
def createBedToolForFilteredList(regionList, createBedFilt, chroms, bedFiltFileName=None):
	# Get the BedTool for a filtered set of peaks
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	if createBedFilt:
		# Create a bed file for a filtered set of peaks
		regionListFilt = bt.BedTool([region for region in regionList if region[0] in chroms])
		regionListFilt.saveas(bedFiltFileName)
	if bedFiltFileName:
		# Save the region list to a file
		regionListFilt = bt.BedTool(bedFiltFileName)
	return regionListFilt
	
def defineInterval(r, halfWindowSize, summitPresent, windowSizeOdd=False):
	# Create an interval that the CNN can take from a region from a bed file
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	chrom = show_value(r[0])
	if summitPresent:
		# Convert the region to a summit-centered interval
		start = int(show_value(r[1])) + int(show_value(r[9])) - halfWindowSize
		if windowSizeOdd:
			# Subtract 1 from the start
			start = start - 1
		end = int(show_value(r[1])) + int(show_value(r[9])) + halfWindowSize
		return [chrom, start, end]
	else:
		# Use the centers of the peaks instead of summits
		start = int(show_value(r[1])) + int(round((float(show_value(r[2])) - float(show_value(r[1])))/2.0)) - halfWindowSize
		if windowSizeOdd:
			# Subtract 1 from the start
			start = start - 1 
		end = int(show_value(r[1])) + int(round((float(show_value(r[2])) - float(show_value(r[1])))/2.0)) + halfWindowSize
		return [chrom, start, end]

def createSetForDeepLearning(genomeFileName, regionList, peakFileNamePrefix, halfWindowSize, summitPresent=True, maxPeakLength=None, \
	chromSizesFileName=None, windowSizeOdd=False, chromEdgeDistLimit=0):
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	chromSizesDict = None
	if chromSizesFileName != None:
		# Create a dictionary mapping chromosomes to their sizes
		chromSizesFile = open(chromSizesFileName)
		chromSizesDict = {}
		for line in chromSizesFile:
			# Iterate through the chromosome sizes and make an entry in the dictionary for each
			lineElements = line.strip().split("\t")
			chromSizesDict[lineElements[0]] = int(lineElements[1])
	intervalList = []
	regionListFiltList = []
	for r in regionList:
		# Convert the list of regions into intervals
		[chrom, start, end] = defineInterval(r, halfWindowSize, summitPresent, windowSizeOdd)
		if start < chromEdgeDistLimit:
			# Do not use the current region because it is too close to the start of the chromosome
			print ("Start < chromEdgeDistLimit for region: " + str(r))
			continue
		if chromSizesDict != None:
			# Check if the current region is too close to the end of the chromosome
			if chrom not in chromSizesDict:
				# The current chromosome is not in the dictionary, so skip it
				print ("Chromosome " + chrom + " is not in the list of chromosomes")
				continue
			if end > chromSizesDict[chrom] - chromEdgeDistLimit:
				# Do not use the current region because it is too close to the end of the chromosome
				print ("End greater than chromosome length - chromEdgeDistLimit for region: " + str(r))
				continue
		if (maxPeakLength != None) and (int(round(float(show_value(r[2])) - float(show_value(r[1])))) > maxPeakLength):
			# The current region is too log, so skip it
			continue
		regionListFiltList.append(r)
		intervalList.append(bt.Interval(chrom, start, end, show_value(r[4])))
	regionListFilt = bt.BedTool(regionListFiltList)
	summitPlusMinus = bt.BedTool(intervalList)
	fastaFileName = None
	if not windowSizeOdd:
		# Do not add 1 to the half window size in the name of the fasta file
		fastaFileName = ".".join([peakFileNamePrefix, "plusMinus" + str(halfWindowSize) + "bp", "fa"])
	else:
		# Add 1 to the half window size in the name of the fasta file
		fastaFileName = ".".join([peakFileNamePrefix, "plusMinus" + str(halfWindowSize + 1) + "bp", "fa"])
	fasta = summitPlusMinus.sequence(fi = genomeFileName, fo = fastaFileName)
	return summitPlusMinus, fastaFileName, regionListFilt
	
def addSequenceToArray(channel1, channel2, channel3, sequenceRecord, perBaseTracks, allData, sampleCount, perBaseTracksDivisor=2):
	# Add a sequence and its reverse complement to a numpy array
	# Also in sequenceOperationsForZnfs.py
	perBaseTracksMat = createPerBaseTracksMat(perBaseTracks, channel3, sampleCount, perBaseTracksDivisor)
	sequenceArray, numNs = oneHotEncode(str(sequenceRecord.seq).strip())
	sequenceArrayReshape = np.reshape(np.vstack((sequenceArray, perBaseTracksMat)), (channel1, channel2, channel3))
	allData[sampleCount,:,:,:] = sequenceArrayReshape
	sampleCount = sampleCount + 1
	# Repeat for the reverse complement
	sequenceArrayRC = reverse_complement(sequenceArray)
	sequenceArrayReshapeRC = np.reshape(np.vstack((sequenceArrayRC, perBaseTracksMat)), (channel1, channel2, channel3))
	allData[sampleCount,:,:,:] = sequenceArrayReshapeRC
	sampleCount = sampleCount + 1
	return [allData, sampleCount]

def addSequenceToArrayNoRC(channel1, channel2, channel3, sequenceRecord, perBaseTracks, allData, sampleCount):
	# Add a sequence and its reverse complement to a numpy array
	# Also in sequenceOperationsForZnfs.py
	perBaseTracksMat = createPerBaseTracksMat(perBaseTracks, channel3, sampleCount, 2)
	sequenceArray, numNs = oneHotEncode(str(sequenceRecord.seq).strip())
	sequenceArrayReshape = np.reshape(np.vstack((sequenceArray, perBaseTracksMat)), (channel1, channel2, channel3))
	allData[sampleCount,:,:,:] = sequenceArrayReshape
	sampleCount = sampleCount + 1
	return [allData, sampleCount]
	
def createPositiveSetFromNarrowPeaks(optimalPeakFileName, genomeFileName, dataShape, createOptimalBed=False, createOptimalBedFilt=True, \
	maxPeakLength=None, chroms=None, chromSizesFileName=None, chromEdgeDistLimit=0):
	# Create the positive set for the deep learning model
	# Also in sequenceOperationsForZnfs.py and sequenceOperationsCore.py
	optimalPeakFileNameElements = optimalPeakFileName.split(".")
	optimalPeakFileNamePrefix = ".".join(optimalPeakFileNameElements[0:-2])
	optimalBedFileName = optimalPeakFileNamePrefix + "_optimal.bed"
	if createOptimalBed:
			# Create a bed file for the optimal peaks
			os.system(" ".join(["zcat", optimalPeakFileName, "| grep -v chrUn | grep -v random | grep chr | sort -k1,1 -k2,2n -k3,3n -k10,10n >", \
					optimalBedFileName]))
	else:
			os.system(" ".join(["zcat", optimalPeakFileName, "| sort -k1,1 -k2,2n -k3,3n -k10,10n >", optimalBedFileName]))
	optimalRegionList = bt.BedTool(optimalBedFileName)
	if chroms != None:
			# Filter for specific chromosomes
			optimalBedFiltFileName = optimalPeakFileNamePrefix + ".train.bed"
			optimalRegionListFilt = createBedToolForFilteredList(optimalRegionList, createOptimalBedFilt, chroms, optimalBedFiltFileName)
	else:
			# Include all of the chromosomes
			optimalRegionListFilt = optimalRegionList
	halfWindowSize = int(dataShape[2]/2)
	windowSizeOdd = False
	if dataShape[2] % 2 > 0:
		# The window size is odd, so put an extra base on the upstream end
		windowSizeOdd = True
	summitPlusMinus, positiveFastaFileName, optimalRegionListFiltPlus =\
			createSetForDeepLearning(genomeFileName, optimalRegionListFilt, optimalPeakFileNamePrefix, halfWindowSize, \
				maxPeakLength=maxPeakLength, chromSizesFileName=chromSizesFileName, windowSizeOdd=windowSizeOdd, \
				chromEdgeDistLimit=chromEdgeDistLimit)
	return optimalPeakFileNamePrefix, optimalRegionList, optimalRegionListFiltPlus, halfWindowSize, summitPlusMinus, positiveFastaFileName

def makePositiveSequenceInputArraysFromFasta(positiveFastaFileName, dataShape=(1,4,1000), labels=np.array([]), RC=False):
	# Convert each peak into a numpy array, where each peak is inputted in fasta format
	# Also in sequenceOperationsForZnfs.py
	if labels.size == 0:
                # Need to get the labels
                cmd = "grep -c '>' $1 " + positiveFastaFileName
		if RC:
			# Count the reverse complements in the labels
                	numPos = int(subprocess.check_output(cmd, shell=True).strip()) * 2
		else:
			# Use only 1 label for each example
			numPos = int(subprocess.check_output(cmd, shell=True).strip())
		labels = np.ones(numPos, dtype=np.int8)
	channel1 = dataShape[0];
        channel2 = dataShape[1];
        channel3 = dataShape[2];
        allData = np.zeros((len(labels), channel1, channel2, channel3), dtype=np.int8)
        sampleCount = 0
        for positiveRecord in SeqIO.parse(positiveFastaFileName, "fasta"):
                # Iterate through the positive fastas and make a numpy array for each
                if RC:
                        # Include the reverse complement
                        [allData, sampleCount] = addSequenceToArray(channel1, channel2, channel3, positiveRecord, [], allData, sampleCount)
                else:
                        # Do not include the reverse complement
                        [allData, sampleCount] = addSequenceToArrayNoRC(channel1, channel2, channel3, positiveRecord, [], allData, sampleCount)
	return allData, labels

def encodeShuffledData(fastaFileName, kmer, copies, shuffledFastaFileName, data, dataIndicesToSkip=[]):
	# Make one-hot-encoded versions of shuffled sequences
	os.system(" ".join(["fasta-shuffle-letters -kmer", str(kmer), "-copies", str(copies), "-dna", fastaFileName, shuffledFastaFileName]))
        dataIndex = 0
        for seqRecord in SeqIO.parse(shuffledFastaFileName, "fasta"):
        	# Iterate through the shuffled fasta file and create a one-hot-encode of each sequence
                if dataIndex in dataIndicesToSkip:
                        # Skip the current index because the corresponding sequence was skipped
                        continue
                sequence = str(seqRecord.seq)
                encodedSequence, numNs = oneHotEncode(sequence)
                data[dataIndex, 0, :, :] = np.transpose(encodedSequence)
                dataIndex = dataIndex + 1
	return data

def makePositiveAndNegativeSequenceInputArraysFromFastas(positiveFastaFileNameList, negativeFastaFileName, dataShape, labels=np.array([]), \
	removeFastas=False, allDataFileName=None, labelsFileName=None):
	# Convert each peak sequence from a fasta file into a numpy array, assumes that tasks do not have overlapping positive sets
	channel1 = dataShape[0];
        channel2 = dataShape[1];
        channel3 = dataShape[2];
	if labels.size == 0:
		# Need to get the labels
		numPos = np.zeros(len(positiveFastaFileNameList), dtype=int)
		for i in range(len(positiveFastaFileNameList)):
			# Get the number of positives in each fasta file
			cmd = "grep -c '>' $1 " + positiveFastaFileNameList[i]
			numPos[i] = int(subprocess.check_output(cmd, shell=True).strip()) * 2
		cmd = "grep -c '>' $1 " +  negativeFastaFileName
		numNeg = int(subprocess.check_output(cmd, shell=True).strip()) * 2
		print("The number of negative examples is: " + str(numNeg))
		labels = np.zeros((np.sum(numPos) + numNeg, len(positiveFastaFileNameList)), dtype=np.int8)
		for i in range(len(positiveFastaFileNameList)):
			# Set the labels for each task
			if i == 0:
				# At the first task, so first labels should be set to positives
				labels[0:numPos[0], 0] = 1
			else:
				# Set the next set of labels to positives
				labels[np.sum(numPos[0:i]): np.sum(numPos[0:i+1]), i] = 1
        allData = np.zeros((labels.shape[0], channel1, channel2, channel3), dtype=np.int8)
        sampleCount = 0
        for positiveFastaFileName in positiveFastaFileNameList:
                # Iterate through the fasta files for the positives and add each to the dataset
                for positiveRecord in SeqIO.parse(positiveFastaFileName, "fasta"):
                        # Iterate through the positive fastas and make a numpy array for each
                        [allData, sampleCount] = addSequenceToArray(channel1, channel2, channel3, positiveRecord, [], allData, sampleCount)
		if removeFastas:
			# Delete the current fasta file
                	os.remove(positiveFastaFileName)
	print("The number of positive examples is: " + str(sampleCount))
	if negativeFastaFileName != None:
		# Negatives have been added
        	for negativeRecord in SeqIO.parse(negativeFastaFileName, "fasta"):
                	# Iterate through the positive fastas and make a numpy array for each
                	[allData, sampleCount] = addSequenceToArray(channel1, channel2, channel3, negativeRecord, [], allData, sampleCount)
		print("The total number of examples is: " + str(sampleCount))
		if removeFastas:
			# Delete the current fasta file
        		os.remove(negativeFastaFileName)
        if not allDataFileName:
                # Do not save the data and the labels
                assert (not labelsFileName)
                return allData, labels
        np.save(allDataFileName, allData)
        np.save(labelsFileName, labels)

def makePositiveAndNegativeSequenceInputArraysFromNarrowPeaks(positivePeakFileNameList, genomeFileName, negativePeakFileName=None, \
	allDataFileName=None, labelsFileName=None, createOptimalBed=True, createOptimalBedFilt=True, dataShape=(1,4,1000), maxPeakLength=None, \
	regression=False, logLabels=False, 
			chroms=["chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", \
				"chr3", "chr4", "chr5", "chr6", "chr7", "chrX"], \
	addFlanksToNegatives=False, chromSizesFileName=None, separationDistance=500, filterRegionsToLength=True, ambiguousPeakFileName=None, \
	addShufflesToNegatives=False, kmer=3, copies=2, datasetName="train", outputFileNamePrefix="", \
	shuffledFastaFileNameSuffix="positivesShuffled.fa", chromEdgeDistLimit=0):
	# Convert each peak into a numpy array, assumes that tasks do not have overlapping positive sets
	positiveFastaFileNameList = []
	halfWindowSize = dataShape[2]/2
	numPositives = np.zeros(len(positivePeakFileNameList))
	currentTaskNum = 0
	positiveRegionListFiltList = []
	for positivePeakFileName in positivePeakFileNameList:
		# Iterate through the positives for each task and get the fasta file for each
		_, _, positiveRegionListFilt, _, _, positiveFastaFileName =\
			createPositiveSetFromNarrowPeaks(positivePeakFileName, genomeFileName, dataShape, createOptimalBed=createOptimalBed, \
				createOptimalBedFilt=createOptimalBedFilt, maxPeakLength=maxPeakLength, chroms=chroms, \
				chromEdgeDistLimit=chromEdgeDistLimit)
		positiveFastaFileNameList.append(positiveFastaFileName)
		numPositives[currentTaskNum] = int(positiveRegionListFilt.count())
		currentTaskNum = currentTaskNum + 1
		positiveRegionListFiltList.append(positiveRegionListFilt)
	posLabels = np.zeros((int(2 * np.sum(numPositives)), len(positivePeakFileNameList)), dtype=np.int8)
	print(posLabels.shape)
	totalPositives = int(0)
	for i in range(len(positivePeakFileNameList)):
		# Iterate through the positive examples and set the labels appropriately
		print("i: " + str(i))
		print("totalPositives: " + str(totalPositives))
		print("numPositives[i]: " + str(numPositives[i]))
		if not regression:
			# Use binary labels
			posLabels[totalPositives:int(totalPositives + (2 * numPositives[i])), i] = 1
		else:
			# Use the signal values
			assert(len(positiveRegionListFiltList[i][0]) >= 7)
			signalsList = []
			for r in positiveRegionListFiltList[i]:
				# Iterate through the positive sets and add the signals for each
				if not logLabels:
					# Use the raw signals
					signalsList.append(float(show_value(r[7])))
					signalsList.append(float(show_value(r[7])))
				else:
					# Use the log of the signals
					signalsList.append(np.log2(float(show_value(r[7]))))
					signalsList.append(np.log2(float(show_value(r[7]))))
			signals = np.array(signalsList)
			posLabels[totalPositives:int(totalPositives + (2 * numPositives[i])), i] = signals
		totalPositives = int(totalPositives + (2 *numPositives[i]))
	negativeRegionListFilt = bt.BedTool()
	negativeFastaFileName = None
	if negativePeakFileName != None:
		negativePeakFileNamePrefix, negativeRegionList, negativeRegionListFilt, halfWindowSize, summitPlusMinus, negativeFastaFileName =\
			createPositiveSetFromNarrowPeaks(negativePeakFileName, genomeFileName, dataShape, createOptimalBed=createOptimalBed, \
				createOptimalBedFilt=createOptimalBedFilt, maxPeakLength=maxPeakLength, chroms=chroms, \
				chromEdgeDistLimit=chromEdgeDistLimit)
		print("The total number of negative regions without flanks is: " + str(negativeRegionListFilt.count()))
	if addFlanksToNegatives:
		# Add flanking regions to the negative set, implemented for only single-task models
		assert(len(positivePeakFileNameList) == 1)
                print("Adding flanking regions")
                flankingRegions =\
                        makeFlankingRegions(positiveRegionListFiltList[0], chromSizesFileName, flankDistance=dataShape[2], \
				separationDistance=separationDistance, filterRegionsToLength=filterRegionsToLength)
        	flankingRegionsFilt =\
                	flankingRegions.subtract(positiveRegionListFiltList[0], A = True)
        	if ambiguousPeakFileName != None:
                	# Remove the ambiguous regions from the flanking regions
			ambigRegions = bt.BedTool(ambiguousPeakFileName[0]).sort()
                	flankingRegionsFilt = flankingRegionsFilt.subtract(ambigRegions, A = True)
		if negativePeakFileName != None:
                        # Remove the existing negatives from the flanking regions
                        flankingRegionsFilt = flankingRegionsFilt.subtract(negativeRegionListFilt, A = True)
        		negativeRegionListFilt = negativeRegionListFilt.cat(flankingRegionsFilt).sort() # Note that cat merges after concatenating
		else:
			# Sort the negative regions, but there are no existing negative regions to combine with the flanking regions
			negativeRegionListFilt = flankingRegionsFilt.sort()
		if not negativeFastaFileName:
			# Create the file name for the negative set
			positivePeakFileNameElements = positivePeakFileName.split(".")
        		positivePeakFileNamePrefix = ".".join(positivePeakFileNameElements[0:-2])
			negativeFastaFileName = positivePeakFileNamePrefix + "_flanks.bed"
        	fasta =\
			negativeRegionListFilt.filter(lambda region: region.end - region.start == dataShape[2]).sequence(fi = genomeFileName, \
				fo = negativeFastaFileName)
	labels = np.empty(1)
	if negativeFastaFileName != None:
		# The negatives have been generated
		negLabels = np.zeros((2 * negativeRegionListFilt.count(), len(positivePeakFileNameList)), dtype=np.int8)
		labels = np.concatenate((posLabels, negLabels), axis=0)
	else:
		# The negatives have not been generated
		assert(addShufflesToNegatives)
		labels = posLabels
	print ("The total number of examples without shuffles is: " + str(labels.shape[0]))
	if addShufflesToNegatives:
                # Add the shuffled data to the negative set, implemented for only single-task models
		assert(len(positivePeakFileNameList) == 1)
                shuffledData = np.zeros((totalPositives * copies, 1, dataShape[2], 4))
		shuffledFastaFileName = outputFileNamePrefix + "_" + datasetName + "_" + shuffledFastaFileNameSuffix
		shuffledData = encodeShuffledData(positiveFastaFileNameList[0], kmer, copies, shuffledFastaFileName, shuffledData)
		negShuffledLabels = np.zeros((totalPositives * copies, len(positivePeakFileNameList)), dtype=np.int8)
		labels = np.concatenate((labels, negShuffledLabels), axis=0)
		print ("The total number of examples with shuffles is: " + str(labels.shape[0]))
		if not negativeFastaFileName:
			# Make the negative fasta file be the shuffled data
			negativeFastaFileName = shuffledFastaFileName
		else:
			# Make the negative fasta file a concatenation of the negative fasta file and the shuffled fasta file
			concatenatedFastaFileName = ".".join(negativeFastaFileName.split(".")[0:-1]) + "_catShuf.fa"
			with open(concatenatedFastaFileName, 'wb') as concatenatedFastaFile:
				# Concatenate the negative fasta files
				for fileName in [negativeFastaFileName, shuffledFastaFileName]:
					# Iterate through the negative fasta files
					with open(fileName, 'rb') as f:
						# Concatenate the current negative fasta file
						shutil.copyfileobj(f, concatenatedFastaFile)
			negativeFastaFileName = concatenatedFastaFileName
	makePositiveAndNegativeSequenceInputArraysFromFastas(positiveFastaFileNameList, negativeFastaFileName, dataShape, labels=labels, removeFastas=True, \
        	allDataFileName=allDataFileName, labelsFileName=labelsFileName)
	
def makeFlankingRegions(optimalRegionList, chromSizesFileName, flankDistance=1000, separationDistance=0, filterRegionsToLength=True):
	# Make flanking regions
	# Also in sequenceOperationsForZnfs.py
	slopRegionList = optimalRegionList.slop(g=chromSizesFileName, b=separationDistance)
	flankingRegionListPre = slopRegionList.flank(g=chromSizesFileName, b=flankDistance)
	flankingRegionList = flankingRegionListPre.subtract(slopRegionList)
	if filterRegionsToLength:
		# Remove the flanks that have been partially removed
		flankingRegionList = flankingRegionList.filter(lambda x: len(x) == flankDistance)
	flankingRegionList = flankingRegionList.sort()
	return flankingRegionList

def makeSequenceInputArraysFromDifferentialPeaks(DESeq2OutputFileName, genomeFileName, backgroundFileName, dataShape, \
                        createOptimalBed=False, backgroundSummitPresent=False, backgroundSummitOnly=True, createModelDir=False, \
                        chroms=["chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr22", \
                                "chr3", "chr4", "chr5", "chr6", "chr7", "chr9", "chrX"], \
                        bigWigFileNames=[], multiMode=False, streamData=False, dataFileName="", RC=True, removeFastas=True, strictNegativeSet=False, \
                        fcCutoff=-1, swapLabels=True, useDESeq2OutputFileNameForFastaFileName=False, chromEdgeDistLimit=0):
        # Convert each peak into a numpy array, where optimal peaks are positives and background peaks that do not overlap relaxed peaks are negatives
        # ASSUMES THAT THIS IS FOR BINARY CLASSIFICATION
        backgroundRegionList = None
        if not backgroundSummitOnly:
                # The background dataset has peaks, not just summits
                backgroundRegionList = bt.BedTool(backgroundFileName).filter(lambda x: x[0] in chroms)
        else:
                # Create regions from the summits
                backgroundRegionList =\
                        makeSummitPlusMinus(backgroundFileName, createOptimalBed=False, dataShape=dataShape, summitCol=1, startCol=None)
        DESeq2Output = np.loadtxt(DESeq2OutputFileName)
        positiveRegionList = []
        negativeRegionList = []
        i = 0
        positiveIndicesList = []
        negativeIndicesList = []
        for backgroundRegion in backgroundRegionList:
                # Iterate through the regions and add the up-regulated ones to the positive set and the down-regulated ones to the negative set
                if backgroundRegion[0] not in chroms:
                        # The current region is not on the list of chromosomes to include
                        i = i + 1
                        continue
                if np.isnan(DESeq2Output[i,1]) or np.isnan(DESeq2Output[i,5]):
                        # The log2 fold change or the padj is nan, so skip this row
                        i = i + 1
                        continue
                if (DESeq2Output[i,1] < fcCutoff) and (DESeq2Output[i,5] < 0.05):
                        # The current row has a large fold-change in the wild-type direction and a small p-value, so add the current region to the negative set
                        negativeRegionList.append(backgroundRegion)
                        negativeIndicesList.append(i)
                elif (not strictNegativeSet) and (DESeq2Output[i,1] >= 0):
                        # The wild-type is not up-regulated, so add the current region to the positive set
                        positiveRegionList.append(backgroundRegion)
                        positiveIndicesList.append(i)
                elif strictNegativeSet and ((DESeq2Output[i,1] > fcCutoff) and (DESeq2Output[i,5] < 0.05)):
                        # The wild-type is down-regulated, so add the current region to the positive set
                        positiveRegionList.append(backgroundRegion)
                        positiveIndicesList.append(i)
                i = i + 1
        print ("Number of positives: " + str(len(positiveRegionList)))
        print ("Number of negatives: " + str(len(negativeRegionList)))
        backgroundFileNameElements = backgroundFileName.split(".")
	halfWindowSize = dataShape[2]/2
        modelDir = ".".join(backgroundFileNameElements[0:-2]) + ".KerasModels"
        if createModelDir:
                # Create a new directory where the model parameters will be written
                os.mkdir(modelDir)
        backgroundPosFileNamePrefix = ".".join(backgroundFileNameElements[0:-2]) + ".Pos"
        backgroundNegFileNamePrefix = ".".join(backgroundFileNameElements[0:-2]) + ".Neg"
        if useDESeq2OutputFileNameForFastaFileName:
                # Use the DESeq2 file name instead of the background file name to create the names of the fasta file
                DESeq2OutputFileNameElements = DESeq2OutputFileName.split(".")
                backgroundPosFileNamePrefix = ".".join(DESeq2OutputFileNameElements[0:-1]) + ".Pos"
                backgroundNegFileNamePrefix = ".".join(DESeq2OutputFileNameElements[0:-1]) + ".Neg"
        summitPlusMinus, positiveFastaFileName, _ =\
                createSetForDeepLearning(genomeFileName, bt.BedTool(positiveRegionList), backgroundPosFileNamePrefix, halfWindowSize,
                        summitPresent=backgroundSummitPresent, chromEdgeDistLimit=chromEdgeDistLimit)
        summitPlusMinus, negativeFastaFileName, _ =\
                createSetForDeepLearning(genomeFileName, bt.BedTool(negativeRegionList), backgroundNegFileNamePrefix, halfWindowSize,
                        summitPresent=backgroundSummitPresent, chromEdgeDistLimit=chromEdgeDistLimit)
        allData, labels, positiveFastaFileName, negativeFastaFileName =\
                makeSequenceInputArraysFromNarrowPeaks(bt.BedTool(positiveRegionList), bt.BedTool(negativeRegionList),
                        positiveFastaFileName, negativeFastaFileName, dataShape, bigWigFileNames=bigWigFileNames, multiMode=multiMode, streamData=streamData,
                        dataFileName=dataFileName, RC=RC, removeFastas=removeFastas)
        positiveIndicesList.extend(negativeIndicesList)
        if swapLabels:
                # Swap the positive and negative labels
                labels = 1 - labels
        bt.helpers.cleanup()
        return allData, labels, positiveFastaFileName, negativeFastaFileName, modelDir, np.array(positiveIndicesList)

def makeSequenceInputArrays(sequenceFileName, labelsFileName, dataShape, logLabels=False, standardizeLabels=False, perBaseTrackFileNames=[], multiMode=False, RC=True, maxFracNs = 1.0):
        # Convert each sequence into a numpy array
        # ASSUMES THAT EACH LABEL CORRESPONDS TO THE SEQUENCE ENTRY WITH THE SAME INDEX
        # ASSUMES THAT THE SEQUENCES ARE LISTS AND NOT IN FASTA FORMAT
        sequenceFile = open(sequenceFileName)
        perBaseTracks = loadPerBaseTracks(perBaseTrackFileNames)
        labelsNoRC = np.loadtxt(labelsFileName)
        labels = np.repeat(labelsNoRC, 2, axis = 0)
        if not RC:
                # Do not use the reverse complements
                labels = labelsNoRC
        if logLabels:
                # Log the values
                labels = labels + 0.0001 # Making sure no labels are 0 before log2ing
                assert(labels.all() > 0)
                labels = np.log2(labels)
        if standardizeLabels:
                # Standardize the values
                labels = (labels - self.meanTrain)/self.stdTrain
        channel1 = dataShape[0];
        channel2 = dataShape[1] + len(perBaseTracks);
        channel3 = dataShape[2];
        allData = np.zeros((len(labels), channel1, channel2, channel3), dtype=np.int8)
        print ("The dimensions of the data are: " + str(allData.shape))
        if perBaseTracks:
                # There are additional per-base tracks that might not be ints
                allData = np.zeros((len(labels), channel1, channel2, channel3), dtype=np.float16)
        sampleCount = 0
        totalNs = 0
        for sequence in sequenceFile:
                # Iterate through the fasta sequences and create the alphabets for the sequence and the reverse complement of each
                perBaseTracksMat = createPerBaseTracksMat(perBaseTracks, channel3, sampleCount, 2)
                sequenceArray, numNs = one_hot_encode(sequence.strip())
                fracNsSequence = float(numNs)/float(dataShape[2])
                if fracNsSequence > maxFracNs:
                        # The current sequence has too high a percentage of N's, so do not include it
                        labels = np.vstack((labels[0:sampleCount], labels[sampleCount + 2:len(labels)]))
                        allData = allData[0:-2,:]
                        continue
                totalNs = totalNs + numNs
                sequenceArrayReshape = np.reshape(np.vstack((sequenceArray, perBaseTracksMat)), (channel1, channel2, channel3))
                allData[sampleCount,:,:,:] = sequenceArrayReshape
                sampleCount = sampleCount + 1
                if RC:
                        # Repeat for the reverse complement
                        sequenceArrayRC = reverse_complement(sequenceArray)
                        sequenceArrayReshapeRC = np.reshape(np.vstack((sequenceArrayRC, perBaseTracksMat)), (channel1, channel2, channel3))
			allData[sampleCount,:,:,:] = sequenceArrayReshapeRC
                        sampleCount = sampleCount + 1
        print("The number of examples is: " + str(sampleCount))
        print("The number of labels is: " + str(labels.shape[0]))
        fracNs = float(totalNs)/float(dataShape[2] * sampleCount)
        if RC:
                # Multiply the fraction of N's by 2 because it does not include the reverse complements
                fracNs = fracNs * 2.0
        print("The fraction of Ns is: " + str(fracNs))
        assert (sampleCount == labels.shape[0])
        sequenceFile.close()
        if multiMode:
                # Re-format the data for multi-moding
                return makeMultiModedData(allData, dataShape, len(perBaseTrackFileNames)), labels
        return allData, labels
	
def updatePeakList(regionsFile, mergePeaks, narrowPeakList, labelList, signalsList, label):
	# Add peaks to the list of peaks
	# ASSUMES SINGLE-TASK MODEL
	atFirstRegion = True
	for line in regionsFile:
		# Iterate through the negative regions and make a narrowPeak for each
		peak = narrowPeak(line)
		if atFirstRegion:
			# At the first negative region, so change the flag
			atFirstRegion = False
		elif mergePeaks and narrowPeakList[-1].overlappingPeak(peak):
			# The current peak overlaps with the previous peak, so merge them
			mergedPeak, peakSignal = narrowPeakList[-1].mergePeak(peak)
			narrowPeakList[-1] = mergedPeak
			signalsList[-1].append(peakSignal)
			continue
		narrowPeakList.append(peak)
		labelList.append(label)
		signalsList.append([peak.signal])
	return narrowPeakList, labelList, signalsList

def loadPeaksAndLabelsFromBed(negativeRegionsFile=None, positiveRegionsFileList=[], binSize=200, sequenceLength=500, binStride=50, maxPeakLength=1000, \
	scoreIndex=6, mergePeaks=True, chromSizesFileName=None):
	# Create positive and negative sets by binning regions
	# negativeRegionsFile: File with negative regions; should be sorted by chromosome, start, end, summit offset; should not include any positive or ambiguous peaks
	# positiveRegionsFileList: List of files with positive regions for each task; should be sorted by chromosome, start, end, summit offset; no regions should overlap the corresponding ambiguous peaks
	# binSize: Size of bins
	# sequenceLength: Length of sequence that will be used for training; flank size = (sequenceLength - binSize) / 2
	# binStride: Stride that will be used for binning
	# maxPeakLength: Maximum length of peaks to use (peaks longer than maxPeakLength will not be used)
	# scoreIndex: Column with peak signals, -1 if the peak files are in 3-column bed format
	# mergePeaks: Merge peaks before binning
	# chromSizesFileName: File with list of chromosome sizes
	# Output: peaksAndLabels -- peaksInfo data structure with one-hot-encoded binned peaks and the corresponding labels and signals
	assert(len(positiveRegionsFileList) > 0)
	assert(len(positiveRegionsFileList) == 1) # Multi-tasking has not yet been implemented
	narrowPeakList = []
	positiveRegionsFile = positiveRegionsFileList[0]
	labelList = []
	signalsList = []
	narrowPeakList, labelList, signalsList = updatePeakList(positiveRegionsFile, mergePeaks, narrowPeakList, labelList, signalsList, 1)
	narrowPeakList, labelList, signalsList = updatePeakList(negativeRegionsFile, mergePeaks, narrowPeakList, labelList, signalsList, 0)
	peakSignalList = [np.mean(np.array(s)) for s in signalsList]
	peakBinList = []
	labelBinList = []
	signalBinList = []
	for peak, label, signal in izip(narrowPeakList, labelList, peakSignalList):
		# Iterate through the peaks and bin each of them
		peakBins = peak.getPeakBins(binSize, sequenceLength, binStride, chromSizesFileName)
		peakBinList.extend(peakBins)
		for i in range(len(peakBins)):
			# Iterate through the bins and add a label and signal for each
			labelBinList.append(label)
			signalBinList.append(signal)
	print("The number of bins is " + str(len(peakBinList)))
	peaksAndLabels = peaksInfo(peakBinList, np.array(labelBinList), np.array(signalBinList))
	return peaksAndLabels
	
def getOneHotEncodedSequences(peaks, genomeFastaFile, RC=True, shuffle=False, fastaFileName=None, shuffledFastaFileName=None, kmer=3, copies=2, removeFastaFiles=True):
	# Get one-hot-encoded sequences for a list of peaks and optionally shuffle them
	# peaks: narrowPeaks with peaks whose sequences need to be encoded
	# genomeFastaFile: pysam FastaFile of file with genomeFastaFile
	# RC: True if reverse complements should be included
	# shuffle: True if sequences should be shuffledFastaFileName
	# fastaFileName: Name of fasta file where sequences will be stored
	# shuffledFastaFileName: Name of fasta file where shuffled sequences will be stored
	# kmer: Nucleotide frequence that will be preserved when shuffling
	# copies: Number of shuffles to create
	# removeFastaFiles: Remove fasta files that get created when shuffling sequences
	sequenceLength = peaks[0].getPeakWidth()
	assert(p.getPeakWidth == sequenceLength for p in peaks)
	data = np.zeros((len(peaks), 1, sequenceLength, 4))
	if not shuffle:
		# Do not shuffle the sequences
		for i, p in enumerate(peaks):
			# Iterate through the peaks and one-hot-encode each
			sequence = genomeFastaFile.fetch(p.chrom, p.start, p.end)
			encodedSequence, numNs = oneHotEncode(sequence)
			data[i, 0, :, :] = np.transpose(encodedSequence)
	else:
		# Shuffle all of the sequences and then encode them
		fastaFile = open(fastaFileName, 'w+')
		for i, p in enumerate(peaks):
			sequence = genomeFastaFile.fetch(p.chrom, p.start, p.end)
			# skip sequences overrunning the contig boundary
			fastaFile.write(">peak" + str(i) + "\n")
			fastaFile.write(sequence + "\n")
		fastaFile.close()
		data = encodeShuffledData(fastaFileName, kmer, copies, shuffledFastaFileName, data)
		if removeFastaFiles:
			# Remove the fasta files
			os.remove(fastaFileName)
			os.remove(shuffledFastaFileName)
	return data
