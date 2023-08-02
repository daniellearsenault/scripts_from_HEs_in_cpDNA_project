#designed for use with python 3 or later (just due to dif. in input functions)

#input: MFannot .new output file
#output: gff3 annotation file

#get input filename from user (MFannot .new output file)
infile = input('Enter mfannot .new output filename:\n')
with open(infile, 'r') as raw_mfa:
    mfa = raw_mfa.readlines()

#get name for .gff3 output from user
outfile_name = input('Enter name for gff3 output file:\n')
#provide extension if user did not
if outfile_name[-4:]!='.gff':
    outfile_name += '.gff'
outfile = open(outfile_name,'w')

#write gff version at top of output file
outfile.write('##gff-version 3\n')

#establish seqid
seqid_i = 0
while mfa[seqid_i][0]!='>':
    seqid_i += 1
seqid_line = mfa[seqid_i]
seqid_stop = seqid_line.find('gc=')
seqid = seqid_line[1:seqid_stop]

#establish transl_table for attributes
transl_table = ''
if seqid_line[-1]=='\n':
    transl_table = seqid_line[(seqid_stop+3):-1]
if seqid_line[-1]!='\n':
    transl_table = seqid_line[(seqid_stop+3):]

#establish source
source = 'mfannot'

#entries will be collected as a dictionary, with the keys being the element names, and their values being the gff3 format vlaues as a list
#entries[elem_name]=[seqid, source, type, start, end, score, strand, phase, attributes]
entries = {}

#return 'start' or 'end' depending on which bound given line corresponds to
def start_or_end(line):
    bound = 'start'
    if ' end' in line:
        bound = 'end'
    return bound

#return the element name for a given line starting with ';' (ex. ';     G-rnl_2 ==> end' returns 'rnl_2')
def get_name(line):
    n_l = 1
    while line[n_l] == ' ':
        n_l += 1
    n_r = n_l
    while line[n_r]!= ' ':
        n_r += 1  
    to_ret = line[n_l:n_r]
    if to_ret[:2] == 'G-':
        to_ret = to_ret[2:]
    return to_ret

#marked_dupes_<input filename>.new, copy of the input file with all duplicates tagged
# created for easy parsing when establishing entry pairs
raw_uniq_mfa_filename = str('marked_dupes_'+infile)
raw_uniq_mfa = open(raw_uniq_mfa_filename,'w')
for i in range(0, seqid_i):
    raw_uniq_mfa.write(mfa[i]) 
seen = {}
#for entries that start with ;;mfannot
misc_counter = 1
used = False
end_round = False
for i in range(seqid_i, len(mfa)):
    line = mfa[i]
    #if normal comment, not ;; mfannot: line
    if (line[0]==';') and (line[:10] != ';; mfannot'):
        #if begin w ;;, delete first ; for easier name retrieval
        if line[:2]==';;':
            line = line[1:]
        name = get_name(line)
        # if start or stop for elem has already been passed
        if name in seen:
            name_r = line.find(name) + len(name)
            if start_or_end(line) == 'start':
                #if on start and start has already been seen at least once
                if (seen[name][0]>0):
                    print('Duplicate element found (start):')
                    seen[name][0]+=1
                    raw_uniq_mfa.write(str(line[:name_r] + '_' + (str(seen[name][0])) + line[name_r+1:]))
                    print(name,' start, copy #',seen[name][0],'\n')
                #if on start but only end has been seen
                if (seen[name][0]==0):
                    seen[name][0]+=1
                    raw_uniq_mfa.write(line)
            if start_or_end(line) == 'end':
                #if on end and end has already been seen at least once
                if (seen[name][1]>0):
                    print('Duplicate element found (end):')
                    seen[name][1]+=1
                    raw_uniq_mfa.write(str(line[:name_r] + '_' + (str(seen[name][1])) + line[name_r+1:]))
                    print(name,' end, copy #',seen[name][1],'\n')
                #if on end but only start has been seen
                if(seen[name][1]==0):
                    seen[name][1]+=1
                    raw_uniq_mfa.write(line)       
                    
        #if elem has not been seen yet
        if name not in seen:
            #if at start position, add elem entry in seen and set start counter to 1
            if start_or_end(line) == 'start':
                seen[name] = [1,0]
            #if at end position, add elem entry in seen and set end counter to 1
            if start_or_end(line) == 'end':
                seen[name] = [0,1]
            #write current line as is, since it is first instance
            raw_uniq_mfa.write(line)

    #pairs ;;mfannot annotations that border random intronic hits, and ;; --> ;
    #gives user message abt potential alternate start positions ("Def by similarity")
    if line[:10] == ';; mfannot':
        if ('Def by similarity' not in line):
            if used and not end_round:
                raw_uniq_mfa.write(str(line[1:10])+'_'+str(misc_counter)+' ==> end '+str(line[11:]))
                misc_counter += 1
                used = False
                end_round = True
                
            if not used and not end_round:
                raw_uniq_mfa.write(str(line[1:10])+'_'+str(misc_counter)+' ==> start '+str(line[11:]))
                used = True
                end_round = True        
        #do not use def by sim as additional gff3 coord, but print message for user
        if('Def by similarity' in line):
            hold = line[(line.find(':')+1):]
            tru_name = get_name(hold)
            print('**NOTE** potential alternate coordinate for '+tru_name+':')
            print(line[:-1])
            print(mfa[i+1])
    end_round = False
    #if non-comment line, write to raw_uniq_mfa
    if line[0]!=';':
        raw_uniq_mfa.write(line)        

#convert raw_uniq_mfa to accessible lines
raw_uniq_mfa.close()
with open(raw_uniq_mfa_filename,'r') as raw_raw_mfa:
    uniq_mfa = raw_raw_mfa.readlines()

#return strand a given element is on (ex. ';     G-rnl_2 ==> end' returns '+')
def get_strand(line):
    if '<==' in line:
        return '-'
    if '==>' in line:
        return '+'  

#return bound position of given line
def get_pos(line):
    case_1 = ((get_strand(line) == '-') and (start_or_end(line) == 'end')) or ((get_strand(line) == '+') and (start_or_end(line) == 'start'))
    case_2 = ((get_strand(line) == '-') and (start_or_end(line) == 'start')) or ((get_strand(line) == '+') and (start_or_end(line) == 'end'))

    to_ret = 0

    #find line below containing position (sometimes not immediately beneath)
    indx = uniq_mfa.index(line)
    to_add = 0
    while uniq_mfa[indx+to_add][0]==';':
        to_add+=1

    #if space is first char on position line
    if uniq_mfa[indx+to_add][0]==' ':
        to_ret = int(get_name(uniq_mfa[indx+to_add]))

    #if space is not first char on position line
    if uniq_mfa[indx+to_add][0]!=' ':
        r_i = 0
        while uniq_mfa[indx+to_add][r_i]!=' ':
            r_i+=1
        to_ret = int(uniq_mfa[indx+to_add][:r_i])
    
    if case_1:
        return to_ret
    if case_2:
        return to_ret - 1

#return element type for a given line
def get_type(line):
    name = get_name(line)
    elem_type = 'gene'
    if 'trn' in name:
        elem_type = 'tRNA'
    if ('rns' in name) or ('rnl' in name):
        elem_type = 'rRNA'
    if '-E' in name:
        elem_type = 'exon'
    if '-I' in name:
        elem_type = 'intron'
    if 'mfannot' in name:
        elem_type = 'misc_feature'
    if 'orf' in name:
        elem_type = 'CDS'
    return elem_type

#return codon frame for phase if applicable (dep on get_type)
def get_phase(line):
    phase = '.'
    if get_type(line)=='CDS':
        phase = '0'
    return phase

#return attributs for a given line (ID, Name, and transl_table separate, rest as 'note')
def get_attributes(line):
    name = get_name(line)
    i = 0
    if 'start' in line:
        i = line.find('start') + 5
    if 'end' in line:
        i = line.find('end') + 3
    att_chunk = line[i:]
    att_elems = str('ID='+name+';Name='+name+';transl_table='+transl_table)
    skip = ['\n', '(', ')']
    space = [';', ':', '=', ' ', '/']
    if len(att_chunk) > 1:
        att_elems += ';note='
        for char in att_chunk:
            if char not in skip:
                if char in space:
                    att_elems += '_'
                if char not in space:
                    att_elems += char
    att_elems += '\n'
    return att_elems

#optional function copies_report:
#reports how many bounds are present for each element (a dupe-free element will return 2: 1 start, 1 end)
#uses *original* input file
def copies_report():
    names = {}
    for i in range (seqid_i, len(uniq_mfa)):
        line = uniq_mfa[i]
        if line[0]==';':
            name = get_name(line)
            if name in names:
                names[name] += 1
            if name not in names:
                names[name] = 1
    for name in names:
        print(str(names[name])+' bounds for '+name)
#copies_report()




#create entries
for i in range(seqid_i, len(uniq_mfa)):
    line = uniq_mfa[i]
    if line[0]==';':
        #get name of current element
        name = get_name(line)
        #if there is already an entry for this element
        if name in entries:
            #get attributes for this line
            cur_att = get_attributes(line)
            #if this attribute contains additional info, add it
            if (len(cur_att) > len(entries[name][8])) and ((entries[name][8]) in cur_att):
                entries[name][8] = cur_att
            #add current bound (start or end)
            if (start_or_end(line)=='start') and (entries[name][3] == '.'):
                entries[name][3] = get_pos(line)
            if (start_or_end(line)=='end') and (entries[name][4] == '.'):
                entries[name][4] = get_pos(line) 
        #if there is not already an entry for this element
        if name not in entries:
            entries[name] = [seqid, source, get_type(line), '.', '.', '.', get_strand(line), get_phase(line), get_attributes(line)]
            if (start_or_end(line)=='start') and (entries[name][3] == '.'):
                entries[name][3] = get_pos(line)
            if (start_or_end(line)=='end') and (entries[name][4] == '.'):
                entries[name][4] = get_pos(line)
                
#write entries in gff3 format to outfile
for name in entries:
    for i in range(8):
        outfile.write(str(str(entries[name][i])+'\t'))
    outfile.write(str(entries[name][8]))
outfile.close()
