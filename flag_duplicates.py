#!/usr/local/bin/python

"""Simple script that reads in a genome start site sorted sam file and filters for nonduplicated and
duplicated records. Duplicated records are further split into PCR duplicates and optical duplicates
(e.g. barcoding errors due to sequence proximity on Illumina chip). Optical duplicates are defined
as records with the same genome start site on the same reference with the same tile and cigar string
that fall within 100 pixels in both the y and x coordinates as defined by the read header. Input file
must have sample ID as first field ie: sample1:911:HJNNYADXX:2:2215:16157:7489
Usage: python flag_duplicates.py sorted.sam
Output files include a list of nonduplicated records, PCR duplicate records (only first instance), 
optical duplicate records, and a crosstabulation of duplicate record groups per sample. Most common
sample per record group (eg most likely origin of duplicate) may be printed to sys.out by uncommenting
region in script.
"""

import pandas as pd
from pandas import DataFrame
import sys
import numpy as np

def find_pcr_opt_dups(dups):
	"""From a dataframe generated using read_sam_get_nondups splits putative PCR duplicates from
	putative optical duplicates
	"""
	#initalize empty lists
	optDups = []
	pcrDups = []
	possibleOptNames = []
	
	#add bin groups and group dataframe
	dupsBin = dups.assign(bin_group = dups['tileCig'] + dups['ref'])
	dupGroup = dupsBin.groupby(['bin_group', 'start'], as_index=False)

	#if there are duplicates but all from same sample == PCR duplicate, else add group name to list
	for name, group in dupGroup:
		grouped_samples = dupGroup.get_group(name)['sampleID']
		if len(grouped_samples.unique()) == 1:
			sameSamp = dupGroup.get_group(name)['record'].reset_index()
			pcrDups.append(sameSamp['record'][0]) #only keeps first instance of the duplicate
		else:
			#Uncomment print statements to get most common sample per duplicate group
			#print "---------------------"
			#print "Duplicates in group:"
			#print name, "\n"
			#print "Most likely from this sample:"
			#print grouped_samples.value_counts().idxmax() #most common id in group (will not work if counts are equal across samples)
			#print "---------------------"
			possibleOptNames.append(name)

	#get crosstab table by sample per bin group
	dupsGroup_counts = dups.assign(duplicates = dups['tileCig'] + "_" + dups['ref'] + "_" + dups['start'])
	count_table = pd.crosstab([dupsGroup_counts.ref, dupsGroup_counts.start, dupsGroup_counts.tileCig], dupsGroup_counts.sampleID, margins=True)
	with open("count_table.txt", "w") as countTable:
		count_table.to_csv(countTable)
	countTable.close()
	
	#process duplicate records
	nit = 0
	print "Processing possible optical duplicates"
	for i in range(0, len(possibleOptNames)):
		xvals = dupGroup.get_group(possibleOptNames[nit]).sort('x')['x'].reset_index()
		yvals = dupGroup.get_group(possibleOptNames[nit]).sort('x')['y'].reset_index()
		records = dupGroup.get_group(possibleOptNames[nit])['record'].reset_index()
		nit += 1
		fin = 0
		sin = 1

		for k in range(0, len(xvals)):
			try:
				t = xvals['x'][sin] #if you have reached last record break the loop
			except KeyError:
				break

			xi = xvals['x'][fin]
			xj = xvals['x'][sin]
			yi = yvals['y'][fin]
			yj = yvals['y'][sin]
			if abs(int(xi)-int(xj)) > 100 or abs(int(yi)-int(yj)) > 100: #if either x or y coordinates more than 100 pixels == PCR duplicate
				if records['record'][fin] not in pcrDups:
					pcrDups.append(records['record'][fin]) #add first record
				if records['record'][sin] not in pcrDups:
					pcrDups.append(records['record'][sin]) #add second record
				fin += 1
				sin += 1
			elif abs(int(yi)-int(yj)) <= 100 and abs(int(xi)-int(xj)) <= 100: #if y and x coordinates under 100 pixels == optical duplicate
				if records['record'][fin] not in optDups:			
					optDups.append(records['record'][fin])
				if records['record'][sin] not in optDups:			
					optDups.append(records['record'][sin])
				fin += 1
				sin += 1
				break		
			else:
				print "Records do not match critiera:"
				print records['record'][fin]
				print records['record'][sin]
	print "Complete!"
	print "Found %i PCR duplicates" % len(pcrDups)
	print "Found %i optical duplicates" % len(optDups)

	with open("pcr_duplicates.txt", "w") as pcrOut:
		for record in pcrDups:
			print >> pcrOut, record
	pcrOut.close()
	with open("optical_duplicates.txt", "w") as optOut:
		for record in optDups:
			print >> optOut, record
	optOut.close()


def read_sam_get_nondups(inputfile):
	"""Loads and extracts data from sorted sam file
	"""
	data = []
	with open(inputfile, "r") as f:
		print "Reading in %s..." % inputfile
		for line in f:
			try:
				record = line.split()
				t=record[2]
			except IndexError:
				break
	
			if line.startswith("@"):
				pass
			else:
				record = line.split()
				ref = record[2]
				start = record[3]
				cigar = record[5]
				read = record[0].split(":")
				sampleID = read[0]
				tile = read[4]
				x = read[5]
				y = read[6]
				tileCig = tile + "_" + cigar
				data.append([tileCig, ref, start, x, y, sampleID, line])
		dfSam = DataFrame(data)
		print "Found %i records in %s" % (len(dfSam), inputfile)
	
	print "Finding nondupicated records..."
	nondupsOut = "nonduplicates.txt"

	dfSam.columns = ['tileCig', 'ref', 'start', 'x', 'y', 'sampleID', 'record']
	dfSam['count'] = dfSam.groupby('start')['start'].transform('count')

	nondups = dfSam[dfSam['count'] == 1]
	print "Found %i nonduplicated samples" % len(nondups)
	nondups.record.to_csv(nondupsOut, index=False, header=False)
	dups = dfSam[dfSam['count'] > 1]
	find_pcr_opt_dups(dups)


def main():
	inputfile = sys.argv[1]
	read_sam_get_nondups(inputfile)

main()

	
		

	
		




		
	


		

