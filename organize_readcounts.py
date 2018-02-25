# Script by Muzi Li, Liam Spurr
import os 
import glob
import sys

#---------------------------------
path = str(sys.argv[1])
out_name = str(sys.argv[2])

columns = ['Chrom', 'Pos', 'Ref', 'AltA', 'AltB']
output = dict()
    
for folder in sorted(glob.glob(os.path.join(path, '*'))): 
    print folder
    sample = folder.split('/')[-1]
    for filename in glob.glob(os.path.join(folder, 'readCounts.tsv')):
        print filename
        table = open(filename,'r')
        table.readline()
        for line in table:
            value = line.strip().split('\t')
            dtype = value[4].split('_')[1]
            R = value[13]
            if ',' in value[3]:
                bases = value[3].strip().split(',')
                for base in bases:
                    index = value[0] + ',' + value[1] + ',' + value[2] + ',' + base
                    if index not in output:
                        output[index] = {sample:{dtype:R}}
                    else:
                        if sample not in output[index]:
                            output[index][sample] = {dtype:R}
                        else:
                            output[index][sample][dtype] = R
            if ',' not in value[3]:
                index = value[0] + ',' + value[1] + ',' + value[2] + ',' + value[3]
                if index not in output:
                    output[index] = {sample:{dtype:R}}
                else:
                    if sample not in output[index]:
                        output[index][sample] = {dtype:R}
                    else:
                        output[index][sample][dtype] = R
        table.close()
    
                
fout = open(out_name+'_readcounts_all.txt','w')
columns.extend(['Sample','Nex', 'Ntr', 'Tex', 'Ttr','Samplenum'])
head = '\t'.join(columns)
fout.write(head + '\n')
               
               
outv = []
for item in output.items():
    for j in sorted(item[1].keys()):
        outv.append(str(item[1][j].get('Nex', 'NA')))
        outv.append(str(item[1][j].get('Ntr', 'NA')))
        outv.append(str(item[1][j].get('Tex', 'NA')))
        outv.append(str(item[1][j].get('Ttr', 'NA')))
        fout.write('\t'.join(item[0].split(',')) + '\t' + item[0].split(',')[-1] + '\t' + j + '\t' + '\t'.join(map(str,outv)) + '\t' + str(len(item[1])) + '\n')
        outv = []
    
    
    
fout.close()
