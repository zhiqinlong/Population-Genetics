#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 22:16:42 2024

@author: zhiqinlong
"""

import sys
import csv
import pandas as pd
if len(sys.argv)!=3:
    sys.stderr.write('Usage:python vcf2fastphase.py vcf fast ')
    sys.exit()
vcf_path=sys.argv[1]
vcf_data=pd.read_csv(vcf_path,sep='\t',compression='gzip',skiprows=25150)
sites=vcf_data.shape[0]
ind_n=vcf_data.shape[1] - 9
f1=open(sys.argv[2],'w')
f1.write(str(ind_n)+'\n')
f1.write(str(sites)+'\n')
n=1
print(vcf_data.head(2))
for col in vcf_data.columns[9:]:
    f1.write('# '+col+' '+str(n)+'\n')
    n+=1
    temp=pd.DataFrame({col:vcf_data[col]})
    temp[['gt','inf']]=temp[col].str.split(':', expand=True,n=1)
    temp[['gt1','gt2']]=temp['gt'].str.split('/|\|', expand=True,n=1)
    gt1_l="".join([str(x) for x in temp['gt1']])
    gt1_l=gt1_l.replace(".", "?")
    gt1_l=gt1_l.replace("None", "?")
    f1.write(gt1_l+'\n')
    gt2_l="".join([str(x) for x in temp['gt2']])
    gt2_l=gt2_l.replace(".", "?")
    gt2_l=gt2_l.replace("None", "?")
    f1.write(gt2_l+'\n')
f1.close()
    
