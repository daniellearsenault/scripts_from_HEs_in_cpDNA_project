import re
import os

from datetime import date
today = date.today()
date = today.strftime("%b-%d-%Y")

#original script purpose:
#input == MSFASTA of blast hits organized by ascending e-value (best at top, worst at bottom), potentially multiple entries per genome
#goal == we only want to keep the single best hit from each genome
#output == copy of MSFASTA where each genome is only represented ONCE, and is represented by its best hit according to evalues

#NOTES:
#!!!script was written based on the format of the FASTA complete seq download option
#       from the web version of BLAST Suite tools example below - be sure to compare
#       your output to this format first to assess any differences.
#       >lcl|Query_12345 Name_NC_123.1
#       SEQUENCE
#save this script in the directory that contains all fasta files (.fst, .fasta, or .fas all accepted) you want remove
#       dupes/keep top hits from.
#unless you edit the script yourself to do differently, it will run on !ALL! fasta
#       files in the current directory. unrelated fasta files in the directory may
#       cause errors.
#output: copy of ea. fasta file with all dupes removed (only 1st entry from ea. genome kept),
#       with today's date at the beginning (ex. 'test1.fst' -> 'Oct-10-2021_test1.fst')

def remove_dupes(filename):
    raw_file = open(filename,'r')
    file = raw_file.readlines()
    output = open(str(date+'_'+filename),'w')

    seen = []

    i = 0
    while i < len(file):
        line = file[i]
        if line[0]=='>':
            if line in seen:
                print('removing duplicate in ',line,' from ',filename)
                i+=1
                while i < len(file) and '>' not in file[i]:
                    i+=1
            if line not in seen:
                seen.append(line)
                output.write(line)
                i+=1
        else:
            output.write(line)
            i+=1
    output.close()

for filename in os.listdir('.'):
    if filename.endswith('.fst') or filename.endswith('.fasta'):
        remove_dupes(filename)
