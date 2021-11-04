import re
import os

from datetime import date
today = date.today()
date = today.strftime("%b-%d-%Y")

#original script purpose:
#input == MSFASTA of hits organized by best to worst evalue, potentially multiple entries per genome
#goal == we only want to keep the single best hit from each genome
#output == copy of MSFASTA where each genome is only represented ONCE, and is represented by its best hit according to evalues

#NOTES:
#save this script in the directory that contains all fasta files you want remove
#       dupes/keep top hits from.
#unless you edit the script yourself to do differently, it will run on ALL fasta
#       files in the current directory. unrelated fasta files in the directory may
#       cause errors.
#script was written based on the format of the FASTA complete seq download option
#       from the web version of BLAST Suite tools example below - be sure to compare
#       your output to this format first to assess any differences.
#       >lcl|Query_12345 Name_NC_123.1
#       SEQUENCE
#output: copy of ea. fasta file with all dupes removed (only 1st entry from ea. genome kept),
#       with today's date at the beginning (ex. 'test1.fst' -> 'Oct-10-2021_test1.fst')

def get_ann_line(L):
    ann_end = L.find('.1')
    name = L[:ann_end+2]
    return name

def remove_dupes(filename):
    if filename.endswith('.fasta') or filename.endswith('.fst'):
        #establish input and output files
        with open(filename, 'r') as raw_fst:
            fst = raw_fst.read()
            fst = fst.replace("\n", "")
        output = open((date+'_'+filename),'w')
        #print('fst,',fst)
        out=''
        names = {}

        name_begins = []
        name_ends = []
        #
        i=0
        for i in range(len(fst)):
            if fst[i] == '>':
                name_begins.append(i)
                print('name_begins',name_begins)
            if i!=0:
                if fst[i-1]+fst[i]=='.1':
                    name_ends.append(i)
                    print('name_ends',name_ends)

        for i in range(len(name_begins)):
            cur_start=name_begins[i]
            cur_stop=name_ends[i]+1
            cur_name=fst[cur_start+17:cur_stop]
            cur_seq=[]

            if i == len(name_begins)-1:
                cur_seq = fst[cur_stop:]
            if i < len(name_begins)-1:
                cur_next_start = name_begins[i+1]
                print('cur_seq = fst[cur_stop:cur_next_start]',fst[cur_stop:cur_next_start])
                cur_seq = fst[cur_stop:cur_next_start]

            if cur_name not in names:
                names[cur_name]=cur_seq

        for key in names:
            print(str('>'+key+'\n'+names[key]+'\n'))
            out+=('>'+key+'\n'+names[key]+'\n')

        output.write(out)
        output.close()


#again make sure you run this script IN the directory containing
#all the .fasta or .fst files you would like to run the script on
for filename in os.listdir('.'):
    remove_dupes(filename)
