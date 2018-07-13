#!/usr/bin/env python3


########################################################################
# File: samJuncs.py
#  executable: samJuncs.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 02/12/2018 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
import pysam

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
        self.parser = argparse.ArgumentParser(description = 'samJuncs.py - lorem ipsium.',
                                             epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -i input_sam ')
        # Add args
        self.parser.add_argument('-i', '--sam_file', type=str, action = 'store', required=True, help='Input SAM/BAM file.')
        self.parser.add_argument('--isSam', action = 'store_true',  default=False, help='File is sam')
        self.parser.add_argument('--is_hisat', action = 'store_true', required=False, default=False,  help='Infer junction strand using minimap tag nomenclature.')
        
                
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Sequence Alignment File
########################################################################

class SAM(object):
    '''
    Handles sequence alignment format file input and output.
    '''

    def __init__(self, inFile=None, isSam=False, isHISAT=False):
        # Attributes
        self.isSam = isSam
        self.inFile = inFile
        self.juncCounts = dict()

        # Start pysam object either as bam or sam object
        if self.isSam == False:
            readSam = "r"
        else:
            readSam = "rb"

        try:
            self.reader = pysam.AlignmentFile(self.inFile, readSam)
        except:
            #File does not exist.
            print("ERROR: Cannot find file %s. Exiting!" % self.inFile, file=sys.stderr)
            sys.exit(1)
    
        self.strandInfo = {0:'+', 16:'-'}
        #print(self.reader.find_introns((read for read in self.reader.fetch() if read.is_reverse)))

        #sys.exit(1)

        if isHISAT:
            self.inferJuncStrand = self.inferHISATJuncStrand
        else:
            self.inferJuncStrand = self.inferMiniJuncStrand

    def inferMiniJuncStrand(self, read):
                    # Next will be junctions
        junctions = list()
        orientation = read.flag

        if orientation == 0:
            orientation = "+"
        elif orientation == 16:
            orientation = "-"
        
        tags = read.get_tags()

        
        juncDir = [x[-1] for x in tags if x[0] == 'ts']

        if not juncDir:
            juncDir = 'nan'

        else:
            juncDir = juncDir[0]
            if orientation == "+" and juncDir == "+":
                juncDir = "+"
            elif orientation == "+" and juncDir == "-":
                juncDir = "-"
            elif orientation == "-" and juncDir == "+":
                juncDir = "-"
            elif orientation == "-" and juncDir == "-":
                juncDir = "+"

        return juncDir


    def inferHISATJuncStrand(self, read):
                    # Next will be junctions
        junctions = list()
        orientation = read.flag

        if orientation == 0:
            orientation = "+"
        elif orientation == 16:
            orientation = "-"
        
        tags = read.get_tags()


        juncDir = [x[-1] for x in tags if x[0] == 'XS']

        return juncDir


    def countJuncs(self):
        '''
        Reads self.reader and reports junction counts for each junction.
        Returns dictionary with junction counts per junction per chromosome.
        '''

        strandInfo = {0:'+', 16:'-'}


        for read in self.reader.fetch():

            try:
                strand = strandInfo[read.flag]
            except:
                continue



            qName = read.query_name
            chromosome = read.reference_name
            
            refPos = read.pos
            refEnd = read.pos
            

            startPos = read.pos
            cigar = read.cigar
    
            juncStrand = self.inferJuncStrand(read)

            if len(juncStrand)<1:
                continue
            else:
                juncStrand = juncStrand[0]

            for num, flagTuple in enumerate(cigar,1):
                flag, length = flagTuple 
                if flag not in [0,2,3,7,8]:
                    continue
                    
                if flag == 3:
                    if chromosome not in self.juncCounts:
                        self.juncCounts[chromosome] = dict()
                    if (refEnd+1, refEnd+length, juncStrand) not in self.juncCounts[chromosome]:
                        self.juncCounts[chromosome][(refEnd+1, refEnd+length, juncStrand)] = int()

                    self.juncCounts[chromosome][(refEnd+1, refEnd+length, juncStrand)] += 1

                refPos = refEnd+length
                refEnd = refPos



        return self.juncCounts


    def readJuncs(self):
        '''
        Returns start, end and junctions from a single read.
        '''

        for read in self.reader.fetch():

            try:
                # Skip unmapped or multimapped reads.
                strand = self.strandInfo[read.flag]
            except:
                continue

            qName = read.query_name
            chromosome = read.reference_name
            
            refPos = read.pos
            refEnd = read.pos
            

            startPos = read.pos
            cigar = read.cigar
            
            # Here is the read starts
            rstart = int(read.pos)

            # Next will be junctions
            junctions = list()
            orientation = read.flag

            # minimap gives junction strand denoted as 'ts'
            # the sign corresponds to the alignment orientation, where + agrees and - disagrees
            try:
                juncDir = read.get_tag('ts')
            except:
                juncDir = None

            # Try to resolve strand by looking for polyA
            if not juncDir:
                left, right = cigar[0], cigar[-1]
                s1, s2 = read.seq[:50], read.seq[-50:]
                #pa = str()
                if ("T"*10 in s1 and left[0] == 4 and left[1] >= 10) and ("A"*10 in s2 and right[0] == 4 and right[1] >= 10):
                    # probably internal priming
                    juncDir = "ambig"

                elif ("T"*10 in s1 and left[0] == 4 and left[1] >= 10):
                    # maps to positive strand but has a rev comp polyA
                    juncDir = "-" if orientation == 16 else "+"
                    #print("anti")
                    #pa = "ppa"
                elif ("A"*10 in s2 and right[0] == 4 and right[1] >= 10):
                    # maps to positive strand but has a sense polyA
                    juncDir = "+" if orientation == 16 else "-"
                    #print("sense")
                    #pa = "ppa"
                else:
                    # no polyA or polyT. Fragment?
                    juncDir = "ambig"
                    #pa = "nan"

            else:
                if orientation == 0 and juncDir == "+":
                    juncDir = "+"
                elif orientation == 0 and juncDir == "-":
                    juncDir = "-"
                elif orientation == 16 and juncDir == "+":
                    juncDir = "-"
                elif orientation == 16 and juncDir == "-":
                    juncDir = "+"

            for num, flagTuple in enumerate(cigar,1):
                flag, length = flagTuple 
                if flag not in [0,2,3,7,8]:
                    continue
                    
                if flag == 3:
                    junctions.append((refEnd, refEnd+length))

                refPos = refEnd+length
                refEnd = refPos

            # Last is the end
            rend = refEnd

            yield (qName, chromosome, rstart, junctions, rend, orientation, juncDir, read.mapq)#, pa)

def main():
    '''
    TDB
    '''
    myCommandLine = CommandLine()
    
    alignmentFile = myCommandLine.args['sam_file']
    isSam = myCommandLine.args['isSam']
    isHISAT = myCommandLine.args['is_hisat']

    sObj = SAM(alignmentFile, isSam, isHISAT)


    d = sObj.countJuncs()

    for c,j in d.items():
        for i in j:
            print(c,i[0]-1, i[1], ".", d[c][i], i[-1],  sep="\t")
    

########################################################################
# Main
# Here is the main program
# 
########################################################################

if __name__ == "__main__":
    main()      
