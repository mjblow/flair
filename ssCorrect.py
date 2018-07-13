#!/usr/bin/env python3


########################################################################
# File: 
#  executable: 
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 05/01/2018 Created
#
########################################################################


########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
import re
from subprocess import Popen, PIPE
from multiprocessing import Pool
from intervaltree import Interval, IntervalTree
import random
from pybedtools import BedTool
import peakutils
########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
	'''
	Handle the command line, usage and help requests.
	CommandLine uses argparse, now standard in 2.7 and beyond. 
	it implements a standard command line argument parser with various argument options,
	and a standard usage and help,
	attributes:
	myCommandLine.args is a dictionary which includes each of the available command line arguments as
	myCommandLine.args['option'] 
	
	methods:
	
	'''
	
	def __init__(self, inOpts=None) :
		'''
		CommandLine constructor.
		Implements a parser to interpret the command line argv string using argparse.
		'''
		import argparse
		self.parser = argparse.ArgumentParser(description = ' TBD - TBD.',
											 epilog = 'Please feel free to forward any questions/concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s -i reads.bed -g annotations.gtf -j other_junctions.bed')
		# Add args
		self.parser.add_argument('-i', "--input_bed", action = 'store', required=True, help='Input reads in bed12 format.')
		self.parser.add_argument('-g', "--annotated_junctions", action = 'store', required=False, help='Gencode annotation file.')
		self.parser.add_argument('-j', "--other_junctions", action = 'store', required=False, help='Junction bed file.')
		self.parser.add_argument('-w', '--wiggle_window', action = 'store', required=False, default = 15, help='Splice site correction window.')
		self.parser.add_argument('-p', "--threads", action = 'store', required=False, default = 2, help='Num threads.')
		
		if inOpts is None :
			self.args = vars(self.parser.parse_args())
		else :
			self.args = vars(self.parser.parse_args(inOpts))


########################################################################
# BED File
########################################################################

class BED12(object):
	'''
	Handles BED format file input and output.
	'''
	def __init__(self, fname=None):
		self.fname = fname
		
		if not os.path.isfile(fname):
			print("%s does not exist. Exiting.", file=sys.stderr)
			sys.exit(1)

	def getLine(self):

		with open(self.fname,'r') as entries:
			for entry in entries:
				cols = entry.rstrip().split()
				self.chrom, self.start, self.end, self.read = cols[0], int(cols[1]), int(cols[2]), cols[3]
				self.score, self.strand, self.c1, self.c2 = int(cols[4]), cols[5], int(cols[6]), int(cols[7])
				self.color, self.exons, self.sizes, self.starts = cols[8], int(cols[9]), [int(x) for x in cols[10].split(",")[:-1]], [int(x) for x in cols[11].split(",")[:-1]] 
				yield cols

	def bed12toJuncs(self):

		junctions = list()
		for num, st in enumerate(self.starts,0):

			if num+1 >= len(self.starts):
				break
			ss1 = self.start + st + self.sizes[num]
			ss2 = self.start + self.starts[num+1]
			junctions.append((ss1,ss2))

		return junctions

	def bed12toExons(self):

		exons = list()
		for num, st in enumerate(self.starts,0):
			c1 = self.start + st
			c2 = c1 + self.sizes[num]
			exons.append((c1,c2))
		return exons

	 

########################################################################
# Functions
########################################################################

def juncsToBed12(start, end, coords):
	sizes, starts = [],[]
	# initial start is 0
	if len(coords) > 0:
		for num,junc in enumerate(coords,0):
			ss1, ss2 = junc
			if num == 0:
				st = 0
				size = abs(start-ss1)
			else:
				st = coords[num-1][1] - start
				size =  ss1 - (st + start)
			starts.append(st)
			sizes.append(size)
		st = coords[-1][1] - start
		size =  end - (st + start)
		starts.append(st)
		sizes.append(size)
		return len(starts), sizes, starts
	else:
		return 1, [end-start], [0] 


def buildOtherDB(ssDB, spliceSites, bedJuncs, wiggle):

	ss = ssDB
		
	with open(bedJuncs,'r') as bedLines:
		for line in bedLines:
			cols = line.rstrip().split()
			chrom, c1, c2, name, score, strand = cols[0], int(cols[1]), int(cols[2]), cols[3], cols[4], cols[5]
			
		if chrom not in ss:
			ss[chrom] = IntervalTree()

		#if c1 not in spliceSites:

		ss[chrom][c1-wiggle:c1+wiggle] = ('other',c1, strand)
	
		#if c2 not in spliceSites:
		ss[chrom][c2-wiggle:c2+wiggle] = ('other',c2, strand)

	return ss

def buildGTFDB(file, wiggle):

	ss = dict()
	exons = dict()
	junctionSet = set()

	with open(file,'r') as lines:
		for l in lines:
			if l[0] == "#":
				continue

			cols = l.split()

			if "exon" == cols[2]:
				chrom, c1, c2, strand =  cols[0], int(cols[3]), int(cols[4]), cols[6]
				txn = cols[11]
				key = (chrom, txn, strand)
				if key not in exons:
					exons[key] = list()
				exons[key].append((c1,c2))
				
	for exonInfo, coords in exons.items():
		chrom, txn, strand = exonInfo
		if strand == "-":
			coords = coords[::-1]

		if chrom not in ss:
			ss[chrom] = IntervalTree()

		for pos,x in enumerate(coords,0):
			if pos+1<len(coords):
				junctionSet.add(x[1])
				junctionSet.add(coords[pos+1][0])
				ss[chrom][x[1]-wiggle:x[1]+wiggle] = ('gtf',x[1],strand)
				ss[chrom][coords[pos+1][0]-wiggle:coords[pos+1][0]+wiggle] = ('gtf',coords[pos+1][0],strand)

	return ss, junctionSet

def resolveJunc(juncs, known):
	for junc in juncs:

		left, right = junc
		print(ss[ch][junctionCoords[0][0]],ss[ch][junctionCoords[0][1]])




def main():
	'''
	TDB
	'''
	myCommandLine = CommandLine()
	
	bed = myCommandLine.args['input_bed']
	gtf = myCommandLine.args['annotated_junctions']
	otherJuncs = myCommandLine.args['other_junctions']
	wiggle = myCommandLine.args['wiggle_window']
	threads = myCommandLine.args['threads']
	
	ssDB, spliceSites = buildGTFDB(gtf, wiggle)
	ssDB = buildOtherDB(ssDB, spliceSites, otherJuncs, wiggle)
	
	data = BED12(bed)

	for i in ssDB["chr1"]:
		print(i) 
	#print(ssDB["chrY"])
	sys.exit(1)

	for line in data.getLine():
		junctionCoords = data.bed12toJuncs()
		print(junctionCoords)
		print(data.read)
		ch, st, end = data.chrom, data.start, data.end
		
		fixed = resolveJunc(ss[ch], junctionCoords)
		blocks, sizes, lens = juncsToBed12(st,end,fixed)
		print(line)
		print(blocks)
		print(sizes)

if __name__ == "__main__":
	main()
