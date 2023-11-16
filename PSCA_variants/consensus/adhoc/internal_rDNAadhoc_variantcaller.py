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
        # /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/
        # rather than 
        # /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/${Bxxx}/internal_${Bxxx}toConsPvalue
    # --graph_dir:
        # /nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/${Bxxx}/
    # the minimum length of matching region (min_rlen):
        # 9000
    # --samplename
        # Bxxx
    



# To run this script:
    # conda activate rDNAVariantCall
    # python3 merge_vaf.py


import os
import glob
import sys
import pandas as pd
import cstag
from VariantCounter import VariantCounter

min_rlen = 9000
reflen = 13332




# load the Bcell ID list
Bcell_IDlist_file = "/nfs/data05/otgk/rDNA/Bcell_IDlist.txt"
Bcell_IDlist = []
with open(Bcell_IDlist_file, "r") as f:
    for line in f:
        Bcell_IDlist.append(line.strip())
# print(Bcell_IDlist)


# for test:
Bcell_IDlist = Bcell_IDlist[:5]

# # Define the Variant counter
# class VariantCounter:
#     # for each site(0, 1, ..., ref_len-1), count the number of variants
#     # variant type: SNV, INS, DEL
#     # and the read depth: DP
#     def __init__(self, ref_len):
#         self.ref_len = ref_len
#         self.variant_count = [0] * ref_len
#         self.variant_count_SNV = [0] * ref_len
#         self.variant_count_INS = [0] * ref_len
#         self.variant_count_DEL = [0] * ref_len
#         self.DP = [0] * (ref_len + 1) # read depth 
        
#         # variant frequency
#         self.SNV_freq = [0] * ref_len
#         self.INS_freq = [0] * ref_len
#         self.DEL_freq = [0] * ref_len
#     def add_variant(self, variant_type, start, end):
#         # variant_type: SNV, INS, DEL
#         # start, end: 0-based
#         # if variant_type == "SNV":
#         if variant_type == "*":
#             for i in range(start, end+1):
#                 self.variant_count_SNV[i] += 1
#         # elif variant_type == "INS":
#         elif variant_type[0] == "+":
#             # start should be equal to end
#             for i in range(start, end+1):
#                 self.variant_count_INS[i] += 1
#         # elif variant_type == "DEL":
#         elif variant_type == "-":
#             for i in range(start, end+1):
#                 self.variant_count_DEL[i] += 1
#         elif variant_type == ":": # match
#             pass
#         else:
#             print("Error: variant_type should be SNV, INS or DEL")
#             sys.exit(1)
#         for i in range(start, end+1):
#             self.variant_count[i] += 1

#     def add_DP(self, start, end):
#         # start, end: 0-based
#         for i in range(start, end+1):
#         # for i in range(start, max(end+1, self.ref_len)):
#             # ここでずるいことをしている。PAF では len = 13332 なのに
#             # ref_start = 0 and ref_end = 13332 というマッピングがある。
#             # そのリードの長さは 13333 になるはず。どういうこと？
#             self.DP[i] += 1
#     def print(self):
#         print(self.variant_count)
#         print(self.variant_count_SNV)
#         print(self.variant_count_INS)
#         print(self.variant_count_DEL)
#         print(self.DP)

#     def calc_freq(self):
#         self.SNV_freq = [x / y for x, y in zip(self.variant_count_SNV, self.DP)]
#         self.INS_freq = [x / y for x, y in zip(self.variant_count_INS, self.DP)]
#         self.DEL_freq = [x / y for x, y in zip(self.variant_count_DEL, self.DP)]
#         # print(self.SNV_freq)

total_variant_counter = VariantCounter(reflen)

df_posfreq_vs_Bcell = pd.DataFrame(columns=["POS"])
df_posfreq_vs_Bcell["POS"] = list(range(reflen))
# df_posfreq_vs_Bcell will be like
    # position, SNV(B003), DEL(B003), INS(B003), SNV(B005), DEL(B005), INS(B005), ..., SNV(B_end), DEL(B_end), INS(B_end)
# load the list of paf files

# this process takes 1 sec for each step. 
    # this is slow but acceptable, so if affords to be slow, necessary to parallelize this process.
for Bcell_ID in Bcell_IDlist:
    print("Bcell_ID: ", Bcell_ID)
    paf_file = f"/nfs/data05/otgk/rDNA/mapping/PSCA_mapping/ref=consensusmapping/mapping_result/{Bcell_ID}toCons.paf"
    output_file_prefix = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/{Bcell_ID}/internal_{Bcell_ID}toCons"
    graph_dir = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/{Bcell_ID}"

    # ------------------------------------------- input -------------------------------------------
    # load the paf file and use only the columns of "cstag"(), "Target start on original strand (0-based)" and "Target end on original strand (0-based)", "Query start on original strand (0-based)", "Query end on original strand (0-based)", "Number of matching bases in the mapping", "Number bases in the mapping"
        # https://github.com/lh3/miniasm/blob/master/PAF.md
    paf_cs_df_ = pd.read_csv(paf_file, sep='\t', header=None, usecols=[7, 8, 22], names=[ "ref_start", "ref_end", "long_cstag"])
    paf_cs_df_["rlen"] = paf_cs_df_["ref_end"] - paf_cs_df_["ref_start"]
    # paf_cs_df_["short_cstag"]: for each long_cstag, apply "cstag.shorten(long_cstag)"
    paf_cs_df_["short_cstag"] = paf_cs_df_["long_cstag"].apply(cstag.shorten)
    # print(paf_cs_df_.head())

    # extract the entries with rlen >= min_rlen
    paf_cs_df = paf_cs_df_[paf_cs_df_["rlen"] >= min_rlen]
    paf_cs_df = paf_cs_df.drop(columns=["long_cstag", "rlen"])
    # print(paf_cs_df.head())


    # ------------------------------------------- caluculate Variant count for each site and DP -------------------------------------------
    Bcell_variant_counter = VariantCounter(reflen)
    # for each site, judge if there is a variant and accumulate the read depth
    for ref_start, ref_end, short_cstag in paf_cs_df.values:
        # ref_start, ref_end: 0-based integer
        # short_cstag: string.
            # cstag is:
                # :n means n bases are matched
                # -seq means seq bases are deleted
                    # seq = [ATCG]+
                # +seq means seq bases are inserted
                    # seq = [ATCG]+
                # *BeforeAfter means Before bases changed After bases (SNV)
                    # Before, After = [ATCG]
        # ~~~~~~~~~~~~~ 1. add DP: ~~~~~~~~~~~~~
        Bcell_variant_counter.add_DP(ref_start, ref_end)


        # ~~~~~~~~~~~~~ 2. manipulate the short_cstag ~~~~~~~~~~~~~
        # parse the short_cstag
        parsed_cstag = cstag.split(short_cstag)
        # print(parsed_cstag)

        # read the cstag
        position = ref_start
        for cstag_unit in parsed_cstag:
            cstag_type = cstag_unit[0]
            # print(f"position: {position}, cstag_unit: {cstag_unit}")
            
            if cstag_type == ":":
                # match
                match_len = int(cstag_unit[1:])
                position += match_len
            elif cstag_type == "-":
                # deletion
                del_seq = cstag_unit[1:]
                del_len = len(del_seq)
                Bcell_variant_counter.add_variant("-", position, position + del_len - 1)
                position += del_len

            elif cstag_type == "+":
                # insertion
                ins_seq = cstag_unit[1:]
                ins_len = len(ins_seq)
                Bcell_variant_counter.add_variant("+", position, position)
                # position += 1
            elif cstag_type == "*":
                # SNV
                snv_unit = cstag_unit[1:]
                # snv_before = snv_unit[0]
                # snv_after = snv_unit[1]
                Bcell_variant_counter.add_variant("*", position, position)
                position += 1
            else:
                print("Error: cstag_type should be :+-*")
                sys.exit(1)
        # print("END with the position: ", position)
            
    # ~~~~~~~~~ 3. add the Bcell_variant_counter to the total_variant_counter ~~~~~~~~~~~
    total_variant_counter.variant_count = [x + y for x, y in zip(total_variant_counter.variant_count, Bcell_variant_counter.variant_count)]
    total_variant_counter.variant_count_SNV = [x + y for x, y in zip(total_variant_counter.variant_count_SNV, Bcell_variant_counter.variant_count_SNV)]
    total_variant_counter.variant_count_INS = [x + y for x, y in zip(total_variant_counter.variant_count_INS, Bcell_variant_counter.variant_count_INS)]
    total_variant_counter.variant_count_DEL = [x + y for x, y in zip(total_variant_counter.variant_count_DEL, Bcell_variant_counter.variant_count_DEL)]
    total_variant_counter.DP = [x + y for x, y in zip(total_variant_counter.DP, Bcell_variant_counter.DP)]
    # total_variant_counter.print()


    # ------------------------------------------- calculate the Variant Frequency and p-value for each site -------------------------------------------
    # calculate the Variant Frequenc (for each ) in the Bcell we are looking at
    Bcell_variant_counter.calc_freq()
    # print(Bcell_variant_counter.SNV_freq)

    # ------------------------------------------- concat dataframe -------------------------------------------
    df_posfreq_vs_Bcell[f"SNV({Bcell_ID})"] = Bcell_variant_counter.SNV_freq
    df_posfreq_vs_Bcell[f"INS({Bcell_ID})"] = Bcell_variant_counter.INS_freq
    df_posfreq_vs_Bcell[f"DEL({Bcell_ID})"] = Bcell_variant_counter.DEL_freq
    # print(df_posfreq_vs_Bcell.head())
    # print(df_posfreq_vs_Bcell.tail())
    # 3 (num of variant) * 307 (num of Bcells) * 13332 (reflen) * 8 byte (float size) = 93 MB


# ------------------------------------------- export -------------------------------------------
# export the dataframe
print(f"writing df_posfreq_vs_Bcell to {output_file_prefix}.vaf.csv")
df_posfreq_vs_Bcell.to_csv(f"{output_file_prefix}.vaf.csv", sep=',', index=False)


# ------------------------------------------- calculate the mean Variant Frequency for each site -------------------------------------------
# calculate the mean Variant Frequency for each site
# mean_variant_freq_snv_j (j : position) = mean(VAF_SNV_ji) (i : Bcell)



mean_vaf_SNV = [sum(df_posfreq_vs_Bcell[f"SNV({Bcell_ID})"][j] for Bcell_ID in Bcell_IDlist) / len(Bcell_IDlist) for j in range(reflen)]
mean_vaf_INS = [sum(df_posfreq_vs_Bcell[f"INS({Bcell_ID})"][j] for Bcell_ID in Bcell_IDlist) / len(Bcell_IDlist) for j in range(reflen)]
mean_vaf_DEL = [sum(df_posfreq_vs_Bcell[f"DEL({Bcell_ID})"][j] for Bcell_ID in Bcell_IDlist) / len(Bcell_IDlist) for j in range(reflen)]


# ------------------------------------------- calculate the log p-value for each site for each Bcell -------------------------------------------

# output_df = pd.DataFrame(columns=["POS", "log_P_SNV(Bcell_ID)", "log_P_INS(Bcell_ID)", "log_P_DEL(Bcell_ID)"