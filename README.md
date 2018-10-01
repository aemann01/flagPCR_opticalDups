# Flag PCR and optical duplicates

Script that calculates the number of pcr and optical duplicates (duplicated records that map to the same start size on the same reference with the same tile and cigar string that fall within 100 pixels of another record) in a sam file. Returns multiple output files -- nonduplicates.txt, pcr_duplicates.txt, optical_duplicates, and count_table.txt (crosstabulation of optical duplicate record groups per sample). 

Input sam file must have the sample ID as the first field separated by ":" i.e., sample1:911:HJNNYADXX:2:2215:16157:7489

```
python flag_duplicates.py sorted.sam
```
