# To run this job script, 
# bash merge_vaf.sh

#!/bin/bash

# To do:
# joint /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/$Bcell_ID/$Bcell_IDtoConsPvalue.multiple_test.bonferroni.accepted.csv
# into /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/merged_vaf.csv
# merged_vaf.csv:
    # 1. rDNA position
    # 2. [Variant_Frequency($Bcell_ID) for each $Bcell_ID]

cat /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/B003/B003toConsPvalue.multiple_test.bonferroni.accepted.csv | tail -n +2 | awk -F, '{print $1"_"$3,$4}' | sort -k 1 > B003.tmp
cat /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/B005/B005toConsPvalue.multiple_test.bonferroni.accepted.csv | tail -n +2 | awk -F, '{print $1"_"$3,$4}' | sort -k 1 > B005.tmp
# join in the first columns
join -e 0  -j 1 -a 1 -a 2  B003.tmp B005.tmp -o auto  | sort -k 1n
    # > /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/merged_vaf.csv

# このままだと、B003 に存在するものだけが join されているっぽいが、そうでないものも join したい
# そして、その上で、Null の時は 0 で fill したい。
# 