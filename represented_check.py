import os
import re     

#original purpose: 
#ea. gene msfasta (ex. atpA.fst) in a directory containing blastp results limited to a certain set of genome accessions
#using acc_list of genome accessions:
#report on:
#1. which genomes were underrepresented for each gene (since we expect at least one from each genome) (genome_rep_results.txt)
#and 2. which genes were not found in each genome (commonly_unrep_genomes.txt)


output = open('genome_rep_results.txt','w')
output2 = open('commonly_unrep_genomes.txt','w')

#acc_list: list, hardcoded list of genome accessions. you could add a bit that reads a txt file but my list was constant and small enough that I just hard-coded.
acc_list = ['NC_008097.1','NC_053616.1','NC_053615.1','NC_053614.1','NC_054192.1','NC_054284.1','NC_053613.1','NC_053612.1','NC_053611.1','NC_050739.1','NC_047441.1','NC_047440.1','NC_047431.1','NC_047430.1','NC_045914.1','NC_045893.1','NC_045059.1','NC_044704.1','NC_043860.1','NC_043776.1','NC_042489.1','NC_042488.1','NC_042487.1','NC_042486.1','NC_042485.1','NC_042484.1','NC_042483.1','NC_042255.1','NC_042251.1','NC_042250.1','NC_042241.1','NC_042182.1','NC_042181.1','NC_032043.1','NC_039920.1','NC_039766.1','NC_039755.1','NC_039754.1','NC_039528.1','NC_039527.1','NC_039526.1','NC_039525.1','NC_039524.1','NC_039523.1','NC_039522.1','NC_039521.1','NC_039520.1','NC_039377.1','NC_037366.1','NC_037365.1','NC_037363.1','NC_037923.1','NC_037922.1','NC_037921.1','NC_037920.1','NC_037919.1','NC_037918.1','NC_034714.1','NC_037480.1','NC_037479.1','NC_037450.1','NC_037367.1','NC_037364.1','NC_037007.1','NC_036806.1','NC_036805.1','NC_036668.1','NC_036137.1','NC_035823.1','NC_034712.1','NC_034711.1','NC_034713.1','NC_034710.1','NC_034709.1','NC_034655.1','NC_034654.1','NC_025524.1','NC_033387.1','NC_031510.1','NC_028586.1','NC_028582.1','NC_028581.1','NC_028578.1','NC_026796.1','NC_026795.1','NC_025544.1','NC_025546.1','NC_025541.1','NC_025538.1','NC_025533.1','NC_025526.1','NC_025534.1','NC_025547.1','NC_025549.1','NC_025536.1','NC_025535.1','NC_032284.1','NC_032109.1','NC_032042.1','NC_031368.1','NC_031367.1','NC_030312.1','NC_030169.1','NC_029807.1','NC_029735.1','NC_029673.1','NC_029676.1','NC_029675.1','NC_029674.1','NC_029672.1','NC_029671.1','NC_029670.1','NC_029669.1','NC_029040.1','NC_024811.1','NC_028587.1','NC_028585.1','NC_028584.1','NC_028583.1','NC_028580.1','NC_028579.1','NW_014013626.1','NC_024829.1','NC_024828.1','NC_023835.1','NC_025531.1','NC_025548.1','NC_008101.1','NC_025530.1','NC_025528.1','NC_025537.1','NC_025529.1','NC_025532.1','NC_025525.1','NC_025527.1','NC_025540.1','NC_025543.1','NC_025545.1','NC_025539.1','NC_023775.1','NC_018569.1','NC_021109.1','NC_020438.1','NC_016732.1','NC_016733.1','NC_015645.1','NC_015359.1','NC_015084.1','NC_014346.1','NC_005353.1','NC_001865.1','NC_013359.1','NC_012978.1','NC_012575.1','NC_012101.1','NC_012097.1','NC_012099.1','NC_011031.1','NC_009681.1','NC_008372.1','NC_008289.1','NC_008114.1','NC_008100.1','NC_008099.1','NC_000927.1']

#unreps: dictionary, format unreps[acc]=<list of genes not found in genome at this accession number>
unreps = {}

#for loop: in dict unreps-- making a key for each accession number, with its value initialized to an empty list
for acc in acc_list:
    unreps[acc]=[]

#genome_rep_check: function, return string result; result: string, count how many times each acc appears in the curent gene msfasta and report as a string
def genome_rep_check(file,name):
    unrep = []
    result = ''
    for acc in acc_list:
        result+=str(acc)+' is represented '+str(file.count(acc))+' times'+'\n'
        if file.count(acc) == 0:
            unrep.append(acc)
    if unrep != []:
        result+=('\n'+'The following genomes are !!NOT!! represented for '+name+': '+str(unrep)+'\n'+str(len(unrep))+' genomes unrepresented in total.')
    return result

#for loop: use genome_rep_check function on all gene msfastas in directory with extension fmt .fst, creates genome_rep_results output
for gene_filename in os.listdir('.'):
    if '.fst' in gene_filename:
        with open(gene_filename,'r') as raw:
            gene_file = raw.read()
        output.write('\n')
        output.write(gene_filename)
        output.write('\n')
        output.write(genome_rep_check(gene_file,gene_filename))
        output.write('\n')

        for acc in acc_list:
            if gene_file.count(acc)==0:
                unreps[acc].append(gene_filename)

#for loop: uses unreps dict to write report of genes that were not found for each genome, creates commonly_unrep_genomes output
for acc in unreps:
    output2.write('*********************** \n'+acc+' is unrepresented in the following gene sets: '+str(unreps[acc])+'\n *********************** \n')
output2.close()
output.close()
