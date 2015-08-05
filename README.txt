###############################
Flag Optical and PCR duplicates
###############################
Author: Allison E Mann
Email: allison.e.mann@ou.edu

This is a simple script that takes in a sam file sorted by genome location and returns
multiple output files: nonduplicates.txt, pcr_duplicates.txt, optical_duplicates.txt, and 
count_table.txt (crosstabulation of duplicate record groups per sample).

Optical duplicates are defined as duplicated records that map to the same start site on the same
reference with the same tile and cigar string that fall within 100 pixels of another record in both 
the x and y coordinates as defined by the read header. Input file must have the sample ID as the first
field (sep ":") eg: sample1:911:HJNNYADXX:2:2215:16157:7489. A short sam example file is included.

Usage: python flag_duplicates.py sorted.sam

The most common sample per record group (or the most likely origin of the duplicate record) may be
printed to standard output by uncommenting a region in the script. Note: if counts per record is 
equivalent across sample the last sample in the group will be flagged as the likely origin. 



