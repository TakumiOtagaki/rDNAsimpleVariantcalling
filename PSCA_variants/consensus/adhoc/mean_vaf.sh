# calculate the mean VAF.
# input:
    # POS,VARIANT_TYPE,B003_VAF,B005_VAF,...
        # 307 Bcells.
# output:
    # POS,mean_VAF,std_VAF
        # 1,0.5,0.1
        # 2,0.3,0.2

echo -n "" > /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.tmp.csv

pos_type_len=$(cat /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/merged_VAF.csv | tail -n +2 | wc -l)

# pos_type
cat /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/merged_VAF.csv | tail -n +2 |\
 cut -d , -f 1 | datamash -t , transpose >> /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.tmp.csv

# mean
cat /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/merged_VAF.csv | tail -n +2 |\
 datamash -t,  transpose | tail -n +2 | \
 datamash -t , mean 1-$pos_type_len >> /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.tmp.csv

# std 
cat /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/merged_VAF.csv | tail -n +2 |\
 datamash -t,  transpose | tail -n +2 | \
 datamash -t , sstdev 1-$pos_type_len >> /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.tmp.csv


# transpose
# header
echo "POS,TYPE,mean_VAF,std_VAF" > /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.csv
cat /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.tmp.csv | \
  datamash -t , transpose | sed -e "s/_/,/g" >> /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.csv

# rm tmp
rm /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/Bcellgroupby_VAF.tmp.csv