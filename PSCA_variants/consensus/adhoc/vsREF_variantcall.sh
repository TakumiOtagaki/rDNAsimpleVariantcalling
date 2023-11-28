# To run this job script, 
# qsub -t 1-307 -tc 20 vsREF_variantcall.sh

#!/bin/bash



# ----------------------------------qsub setting---------------------------------
#$ -S /bin/bash
#$ -N Bxxx_pvalue
#$ -cwd 
#$ -e ./log
#$ -o ./out
#$ -l hostname=a03|a02
#$ -q all.q
#$ -pe smp 1


# ---------------------------------module load-----------------------------------
. ~/.bashrc
module purge
ml samtools
conda activate vsREF_rDNAVariantCall




# ----------------------------------変数定義-----------------------------------------
# $SGE_TASK_ID = 1, 2, ..., 307
# Bxxx = Bcell_IDlist[SGE_TASK_ID]
export Bxxx=$(cat /nfs/data05/otgk/rDNA/PSCA_variants/consensus/Bcell_IDlist.txt | sed -n "${SGE_TASK_ID},${SGE_TASK_ID}p")
# echo "${SGE_TASK_ID}"
# echo $Bxxx
mkdir -p /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx}

# python rDNAadhoc_variantcaller.py --reflen 13332 \
# -I /nfs/data05/otgk/rDNA/PSCA_mapping/consensus/mapping/B003toCons.paf\
# -O /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/B003/B003toCons.csv\
# --graph_dir ./result/B003 \
# --rlen_min 9000

# ひとまず一回目の解析フローは以下の parameter で行う。
python vsREF_rDNAadhoc_variantcaller.py --reflen 13332 \
    -I /nfs/data05/otgk/rDNA/PSCA_mapping/consensus/mapping/${Bxxx}toCons.paf \
    -O /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx}/vsREF_${Bxxx}toConsPvalue \
    --graph_dir /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx} \
    --rlen_min 9000 \
    --graph_out \
    --sample_name ${Bxxx} \
    --multiple_test_correction bonferroni \
    --seq_error_rate 0.01

# python rDNAadhoc_variantcaller.py --reflen 13332 \
#     -I /nfs/data05/otgk/rDNA/PSCA_mapping/consensus/mapping/${Bxxx}toCons.paf\
#     -O /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx}/${Bxxx}toConsPvalue.csv\
#     --graph_dir /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx} \
#     --rlen_min 9000 \
#     --graph_out \
#     --sample_name ${Bxxx} \
#     --multiple_test_correction fdr_bh \
#     --seq_error_rate 0.01

echo "${Bxxx} pvalue calculation done."

