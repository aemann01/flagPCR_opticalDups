#!/usr/local/bin/python

"""Reads in a sam file and filters for nonduplicated and duplicated records. Duplicated records are further split into PCR duplicates and optical duplicates. Usage: python flag_duplicates.py sorted.sam
"""

import pandas as pd
from pandas import DataFrame
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from matplotlib.font_manager import FontProperties
import random

#print colored text to screen (should work linux, OS X, and windows -- only tested on linux)
class bcolors:
	WARNING = '\033[91m'
	COMPLETE = '\033[92m'
	
def find_pcr_opt_dups(dups, pixDist):
	#From a dataframe generated using read_sam_get_nondups splits putative PCR duplicates from putative optical duplicates
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
			possibleOptNames.append(name)

	#get crosstab table by sample per bin group
	dupsGroup_counts = dups.assign(duplicates = dups['tileCig'] + "_" + dups['ref'] + "_" + dups['start'])
	count_table = pd.crosstab([dupsGroup_counts.ref, dupsGroup_counts.start, dupsGroup_counts.tileCig], dupsGroup_counts.sampleID, margins=True)
	with open("count_table.txt", "w") as countTable:
		count_table.to_csv(countTable)
	countTable.close()
	
	#process duplicate records
	nit = 0
	print "Processing possible optical duplicates..."
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
			else:
				print "Records do not match critiera:"
				print records['record'][fin]
				print records['record'][sin]
	print "Found %i PCR duplicates" % len(pcrDups)
	print "Found %i optical duplicates" % len(optDups)

	with open("pcr_duplicates.txt", "w") as pcrOut:
		for record in pcrDups:
			pcrOut.write(record)
	pcrOut.close()
	with open("optical_duplicates.txt", "w") as optOut:
		for record in optDups:
			optOut.write(record)
	optOut.close()
	
	#if optical duplicates detected, make figures?
	if len(optDups) >= 1:
		print "-----------------------"
		next = raw_input("Generate tile figures?: ")
		if 'y' in next:
			plotOpticalDups()
		else:
			print bcolors.COMPLETE + "Complete, no figures generated"

def plotOpticalDups():
	#read in optical duplicates file, extract info
	data = []
	with open("optical_duplicates.txt", "r") as f:
		for line in f:
			if not line.strip():
				continue
			else:
				try:
					optRecord = line.split()
					t = optRecord[2]
				except IndexError:
					break
				optRecord = line.split()
				read = optRecord[0].split(":")
				tile = read[4]
				x = read[5]
				y = read[6]
				sampleID = read[0]
				data.append([tile, sampleID, x, y])
			dfOptDups = DataFrame(data)
			dfOptDups.columns = ['tile', 'sampleID', 'x', 'y']
	tileGroups = dfOptDups.groupby('tile')
	tiles = []
	for name, group in tileGroups:
		tiles.append(name)
	
	#initalize values
	git = 0 #group iterator
	fontP = FontProperties()
	fontP.set_size('small')

	#get colors, x y values for group
	for i in range(0, len(tiles)):
		colorDict = {}
		currentGroup = tileGroups.get_group(tiles[git])
		currentName = currentGroup.tile[0]
		x = map(int, currentGroup.x.tolist())
		y = map(int, currentGroup.y.tolist())
		xnorm = [float(i)/sum(x) for i in x]
		ynorm = [float(i)/sum(y) for i in y]
		
		for i in range(0, len(currentGroup.sampleID.unique())):
			colorDict[currentGroup.sampleID.unique()[i]] = matplotlib.colors.cnames.values()[random.randrange(1, 150)] #link random color to each sample
		#loop to create figures	
		fig = plt.figure()
		backgroundColor = '#B9B9B9'
		ax = fig.add_subplot(111, axisbg=backgroundColor)
		ax.scatter(xnorm, ynorm, color=colorDict.values(), marker=',', s=5)
		ax.spines['top'].set_visible(False) #remove plot borders
		ax.spines['right'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		plt.xlim(0, 1.0)
		plt.ylim(0, 1.0)
		plt.tick_params(bottom='off', top='off', right='off', left='off') #no tick marks
		plt.grid(color="white")
		markers = [plt.Line2D([0,0], [0,0], color=color, marker='s', linestyle='') for color in colorDict.values()]
		lgd = plt.legend(markers, colorDict.keys(), numpoints=1, prop=fontP, bbox_to_anchor=(0.5, -0.1), ncol=5, fancybox=True)
		fig.savefig(currentName + '.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight', format='pdf')
		git += 1	
	print bcolors.COMPLETE + "Complete! Figures successfully generated"

def read_sam_get_nondups(inputfile):
	#Loads and extracts data from sam file
	data = []
	header = []
	with open(inputfile, "r") as f:
		print "Reading in %s..." % inputfile
		for line in f:
			try:
				record = line.split()
				t=record[2]
			except IndexError:
				break
			if line.startswith("@"):
				header.append(line)
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

	dfSam.columns = ['tileCig', 'ref', 'start', 'x', 'y', 'sampleID', 'record']
	dfSam['count'] = dfSam.groupby('start')['start'].transform('count')
	
	nondups = dfSam[dfSam['count'] == 1]['record']
	print "Found %i nonduplicated samples" % len(nondups)

	with open("nonduplicates.txt", "w") as f:
		f.write(''.join(header)) #add header to nonduplicates file
		for line in nondups:
			f.write(line)
	dups = dfSam[dfSam['count'] > 1]
	
	if len(nondups) == len(dfSam): #break if no duplicates are found
		print "-----------------------"
		print bcolors.WARNING + "Exit: No duplicated records found"
		sys.exit(1)

	#set pixel distance
	print "-----------------------"
	next = raw_input("Set pixel distance to value other than default (100)? ")
	if 'y' in next:
		pixDist = raw_input("New pixel distance: ")
	else:
		pixDist = 100
	find_pcr_opt_dups(dups, pixDist)

def main():
	print """
  ___ _             ___            _ _         _          
 | __| |__ _ __ _  |   \ _  _ _ __| (_)__ __ _| |_ ___ ___
 | _|| / _` / _` | | |) | || | '_ \ | / _/ _` |  _/ -_|_-<
 |_| |_\__,_\__, | |___/ \_,_| .__/_|_\__\__,_|\__\___/__/
            |___/            |_|                          
	""" 
	print "Author: Allison E Mann (allison.e.mann@ou.edu)"
	print "Usage: python flag_duplicates.py input.sam \n"
	if len(sys.argv) != 2:
		print bcolors.WARNING + "Error! No sam file specified"
		sys.exit(1)
	inputfile = sys.argv[1]
	assert os.path.exists(inputfile), bcolors.WARNING + 'Error! File does not exist: %s. Is the path correct?' %inputfile
	read_sam_get_nondups(inputfile)

main()

	
		

	
		




		
	


		

