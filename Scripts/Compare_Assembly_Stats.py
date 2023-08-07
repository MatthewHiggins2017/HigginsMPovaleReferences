import sys
import numpy
import numpy as np
import pandas as pd
import regex as re
import matplotlib.pyplot as plt

def FastaToDict(InputFile):
	# Converts Fasta File to Dict
	fastadict = {}
	with open(InputFile) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				active_sequence_name = line[1:].split()[0]
				if active_sequence_name not in fastadict:
					fastadict[active_sequence_name] = ''
				continue
			sequence = line
			# If Uracil Present, Substitute for T (I.e. stick with 4 DNA bases)
			fastadict[active_sequence_name] += re.sub('U','T',sequence.upper())
	return fastadict


def SplitAssembly(self):
	# Split assembly at each gap point where
	# 3 N nucleotides are present.
	ContigsList = []
	for Id, Seq in self.AssemblyDict.items():
		ContigsList += Seq.upper().split('NNN')
	ContigsList = [C for C in ContigsList if (C.count('N') != len(C) and len(C)>0)]
	return ContigsList


def AssemblySize(self):
	AS = np.sum([len(x) for x in list(self.AssemblyDict.values())])
	return AS



def GapCount(self):
	GC = len(self.Contigs) - len(self.AssemblyDict)
	return GC


def NCount(self):
	NC = 0
	for Id, Seq in self.AssemblyDict.items():
		NC += Seq.upper().count('N')
	return NC

def GCContent(self):
	PGC = 0
	TC = 0
	for Id, Seq in self.AssemblyDict.items():
		PGC += Seq.upper().count('C')
		PGC += Seq.upper().count('G')
		TC += len(Seq)
	GC = (PGC/TC) * 100
	return GC


def ContigLen(self):
	ContigLen = sorted([len(x) for x in self.Contigs], reverse=True)
	return ContigLen


def N50(self):
	csum=numpy.cumsum(self.ContigLengths)
	n2=int(sum(self.ContigLengths)/2)
	# get index for cumsum >= N/2
	csumn2=min(csum[csum >= n2])
	ind=numpy.where(csum == csumn2)
	n50 = self.ContigLengths[int(ind[0])]
	return n50




################################################################################
################################################################################


class Assembly:
	def __init__(self,ID,fastapath):
		self.ID = ID
		self.AssemblyPath = fastapath
		self.AssemblyDict = FastaToDict(fastapath)
		self.AssemblySize = AssemblySize(self)
		self.ScaffoldCount = len(self.AssemblyDict)
		self.Contigs = SplitAssembly(self)
		self.GapCount = GapCount(self)
		self.NCount = NCount(self)
		self.GCcontent = GCContent(self)
		self.ContigLengths = ContigLen(self)
		self.N50 = N50(self)





ID = sys.argv[1]
FastaPath = sys.argv[2]
Assembly(ID,FastaPath)

################################################################################
################################################################################
################################################################################
################################################################################
