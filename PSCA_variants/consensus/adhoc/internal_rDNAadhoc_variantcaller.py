# This program calculates the Variant Frequency of each variant in the Bcell (given ID).
# Then, for each variant, it calculates the mean Variant Frequency of all Bcells.
# Finally, for each variant, p-values is calculated by comparing the mean Variant Frequency with the Variant Frequency of the Bcell.
    # Variant Frequency of variant j in Bcell i: VAF_ij
    # VAF_ij ~ poisson(lambda_j)
    # lambda_j = mean(VAF_ji) = mean(VAF_j1, VAF_j2, ..., VAF_jn)
    # p_value_ij = poisson(VAF_ij, lambda_j)
    # poisson(VAF_ij, lambda_j) = exp(-lambda_j) * lambda_j^VAF_ij / (VAF_ij)!
        # VAF_ij: float
        # So, (VAF_ij)! should be replaced with gamma(VAF_ij + 1)
# input: 
    # paf file
        # /nfs/data05/otgk/rDNA/mapping/PSCA_mapping/ref=consensusmapping/mapping_result/{Bcell_ID}toCons.paf
    # output_file_prefix
        # /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx}/internal_${Bxxx}toConsPvalue
    # --graph_dir:
        # /nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/result/${Bxxx}
    # the minimum length of matching region (min_rlen):
        # 9000
    # --samplename
        # Bxxx
    



# To run this script:
    # conda activate rDNAVariantCall
    # python3 merge_vaf.py


import os
import glob
import pandas as pd
import cstag

min_rlen = 9000





# load the Bcell ID list
Bcell_IDlist_file = "/nfs/data05/otgk/rDNA/Bcell_IDlist.txt"
Bcell_IDlist = []
with open(Bcell_IDlist_file, "r") as f:
    for line in f:
        Bcell_IDlist.append(line.strip())
print(Bcell_IDlist)

# Define the Variant counter
class VariantCounter:
    # for each site(0, 1, ..., ref_len-1), count the number of variants
    # variant type: SNV, INS, DEL
    # and the read depth: DP
    def __init__(self, ref_len):
        self.ref_len = ref_len
        self.variant_count = [0] * ref_len
        self.variant_count_SNV = [0] * ref_len
        self.variant_count_INS = [0] * ref_len
        self.variant_count_DEL = [0] * ref_len
        # self.DP = [0] * ref_len # read depth <-- not used because of slow speed
        # instead of the above, use the segment tree:
        self.DP_segtree = SegmentTree(ref_len)
    def add_variant(self, variant_type, start, end):
        # variant_type: SNV, INS, DEL
        # start, end: 0-based
        if variant_type == "SNV":
            for i in range(start, end+1):
                self.variant_count_SNV[i] += 1
        elif variant_type == "INS":
            for i in range(start, end+1):
                self.variant_count_INS[i] += 1
        elif variant_type == "DEL":
            for i in range(start, end+1):
                self.variant_count_DEL[i] += 1
        else:
            print("Error: variant_type should be SNV, INS or DEL")
            sys.exit(1)
        for i in range(start, end+1):
            self.variant_count[i] += 1
    def add_DP(self, start, end, DP):
        # start, end: 0-based
        # for i in range(start, end+1):
        #     self.DP[i] += DP




# load the list of paf files
for Bcell_ID in Bcell_IDlist:
    paf_file = f"/nfs/data05/otgk/rDNA/mapping/PSCA_mapping/ref=consensusmapping/mapping_result/{Bcell_ID}toCons.paf"
    output_file_prefix = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/{Bcell_ID}/internal_{Bcell_ID}toConsPvalue"
    graph_dir = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/{Bcell_ID}"

    # ------------------------------------------- input -------------------------------------------
    # load the paf file and use only the columns of "cstag"(), "Target start on original strand (0-based)" and "Target end on original strand (0-based)", "Query start on original strand (0-based)", "Query end on original strand (0-based)", "Number of matching bases in the mapping", "Number bases in the mapping"
        # https://github.com/lh3/miniasm/blob/master/PAF.md
    paf_cs_df_ = pd.read_csv(paf_file, sep='\t', header=None, usecols=[7, 8, 22], names=[ "ref_start", "ref_end", "long_cstag"])
    paf_cs_df_["rlen"] = paf_cs_df_["ref_end"] - paf_cs_df_["ref_start"]
    # paf_cs_df_["short_cstag"]: for each long_cstag, apply "cstag.shorten(long_cstag)"
    paf_cs_df_["short_cstag"] = paf_cs_df_["long_cstag"].apply(cstag.shorten)
    print(paf_cs_df_.head())

    # extract the entries with rlen >= min_rlen
    paf_cs_df = paf_cs_df_[paf_cs_df_["rlen"] >= min_rlen]
    paf_cs_df.drop(columns=["long_cstag", "rlen"], inplace=True)
    print(paf_cs_df.head())


    # ------------------------------------------- caluculate Variant count for each site and 
    break
