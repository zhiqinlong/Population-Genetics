# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 21:33:15 2021

@author: Administrator
"""
import os
import csv
import sys
import re
import pandas as pd
if len(sys.argv)!=3:
    sys.stderr.write('Usage:python vcf2sweepfinder2.py outgroup_dir vcf_dir rec_rates_dir ')
    sys.exit()
    
sweep_path=sys.argv[2]+"/sweep_input/"
if not os.path.exists(sweep_path):
    os.makedirs(sweep_path)

def outgroup_gt(outgroup_file_path):
    outgroup_dict={}
        outgroup=open(outgroup_file_path,'r')
        outgroup_lines=csv.reader(outgroup,delimiter='\t')
        for line in outgroup_lines:
                if any('#' in strings for strings in line ):
                    continue
                else:            
                    gts=[haplotype.split(":")[0] for haplotype in line[9:]]
                    if [True for gt in gts if gt in ['./.', '0/1' ,'1/0','0|1','1|0']]:
                        continue
                    else:
                        chrom,pos,id,ref,alt,qual,filter,info,format=line[0:9]
                        gts=[ ref+'/'+ref if gt=='0/0' else alt+'/'+alt for gt in gts]
                        id=chrom+":"+pos
                        outgroup_dict[id]=gts
        outgroup.close()
        return(outgroup_dict)
  
files = [f for f in os.listdir(sys.argv[2]) if 'vcf' in f]
for file in files:
    chr_list=re.split('\.',file)
    outgroup_file=sys.argv[1]+chr_list[0]+'.DP10.GQ10.highmappability.biallelic.recode.vcf'
    outgroup_dict=outgroup_gt(outgroup_file)
    rec_path=sys.argv[3]+chr_list[0]+"_rec_rates"
    rec=pd.read_csv(rec_file_path,sep='\t')
    vcf_path=sys.argv[2]+file
    outfile_path=sweep_path+chr_list[0]
    vcf=open(vcf_path,'r')
    vcf_data=csv.reader(vcf,delimiter='\t')
    sweep_out = open(outfile_path, 'w')
    sweep_out.write("position\tx\tn\tfolded\n") # write header
    ids=[]
    for line in vcf_data:
        if any('#' in strings for strings in line ):
            continue
        chrom,pos,id,ref,alt,qual,filter,info,format=line[0:9]
        id=chrom+":"+pos
        ids.append(id)
        gts=[haplotype.split(":")[0] for haplotype in line[9:]]
        n=len(gts)*2-gts.count('./.')*2
        x=gts.count('1/1')*2+gts.count('0/1')+gts.count('1/0')+gts.count('1|0')+gts.count('0|1')
        if id not in outgroup_dict.keys():
            if x!=0 and x!=n:
                sweep_out.write(str(pos)+"\t"+str(n-x)+"\t"+str(n)+"\t"+'1'+"\n")
        elif ref+'/'+ref in outgroup_dict[id]:
                ###alt is derived 
                sweep_out.write(str(pos)+"\t"+str(x)+"\t"+str(n)+"\t"+'0'+"\n")
        elif alt+'/'+alt in outgroup_dict[id]:
                ###ref is derived 
                sweep_out.write(str(pos)+"\t"+str(n-x)+"\t"+str(n)+"\t"+'0'+"\n")
        elif x!=0 and x!=n:
                ###can't polarize
                sweep_out.write(str(pos)+"\t"+str(n-x)+"\t"+str(n)+"\t"+'1'+"\n")
        out_rec_path=sweep_path+chr_list[0]+"_rec_rates"
        rec[rec['position'].isin(ids)].to_csv(out_rec_path,index=False)   
    sweep_out.close()
    vcf.close()