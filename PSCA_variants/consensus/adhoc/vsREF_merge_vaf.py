# To run this script:
# python3 merge_vaf.py

import pandas as pd
import numpy as np
import os
import glob
import re
import sys


# Get the Bcell_list (/nfs/data05/otgk/rDNA/PSCA_variants/consensus/Bcell_IDlist.txt)
Bcell_list = pd.read_csv('/nfs/data05/otgk/rDNA/Bcell_IDlist.txt', sep='\t', header=None).values.reshape(-1).tolist()
# print(Bcell_list)
# Get the list of files
df_merged = pd.DataFrame()
for Bcell in Bcell_list:
    vaf_file = f'/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/{Bcell}/{Bcell}toConsPvalue.multiple_test.bonferroni.accepted.csv'
    # vaf file
    # POS,log_P_value_adj,Variant_type,Variant_Frequency,DP,variant_count
    # 48,-43.02886131838572,SNV,0.07513416815742398,559,42
    vaf_df = pd.read_csv(vaf_file, sep=',', header=0, usecols=['POS', 'Variant_type', 'Variant_Frequency'])
    vaf_df['POS_Variant_type'] = vaf_df['POS'].astype(str) + '_' + vaf_df['Variant_type']
    vaf_df = vaf_df[['POS_Variant_type', 'Variant_Frequency']]
    # set index "POS_Variant_type"
    vaf_df = vaf_df.set_index('POS_Variant_type')
    vaf_df = vaf_df.rename(columns={'Variant_Frequency': f"{Bcell}_VAF"})
    # print(len(vaf_df))
    df_merged = pd.concat([df_merged, vaf_df], axis=1)


# NULL -> zero
df_merged = df_merged.fillna(0)
# sort by position: df_merged['POS'] = df_merged['POS_Variant_type'].str.split('_').str[0].astype(int)

df_merged['POS'] = df_merged.index.str.split('_').str[0].astype(int)
df_merged = df_merged.sort_values(by=['POS'])
df_merged = df_merged.drop(columns=['POS'])



print(df_merged)
    
# Save the merged VAF file
df_merged.to_csv('/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/merged_VAF.csv', sep=',', index=True)