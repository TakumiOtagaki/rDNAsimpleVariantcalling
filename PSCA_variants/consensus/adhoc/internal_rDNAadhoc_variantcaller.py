# This program calculates the Variant Frequency of each variant in the Bcell (given ID).
# Then, for each variant, it calculates the mean Variant Frequency of all Bcells.
# Finally, for each variant, p-values is calculated by comparing the mean Variant Frequency with the Variant Frequency of the Bcell.
    # Variant Frequency of variant j in Bcell i: VAF_ij (0 ~ 1)
    # VAF_ij ~ poisson(lambda_j)
    # lambda_j = mean(VAF_ji) = mean(VAF_j1, VAF_j2, ..., VAF_jn) <- this is false;
    # mu_j = mean(VAF_ji) = mean(VAF_j1, VAF_j2, ..., VAF_jn) <- Maximum Likelihood Estimation
    # And then,
    # lambda_j = mu_j * DP_ij

    # p_value_ij = poisson(VariantCount_ij, lambda_j)
    # poisson(VariantCount, lambda_j) = exp(-lambda_j) * lambda_j^VariantCount / (VariantCount)!
        # VariantCount: float
        # So, (VariantCount)! should be replaced with gamma(VariantCount + 1)

# input: 
    # paf file
        # /nfs/data05/otgk/rDNA/mapping/PSCA_mapping/ref=consensusmapping/mapping_result/{Bcell_ID}toCons.paf
        # with the cstag:
            # short_cstag: string.
            # cstag is:
                # :n means n bases are matched
                # -seq means seq bases are deleted
                    # seq = [ATCG]+
                # +seq means seq bases are inserted
                    # seq = [ATCG]+
                # *BeforeAfter means Before bases changed After bases (SNV)
                    # Before, After = [ATCG]
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
# from mystats import log_continuous_poisson
from mystats import log_continuous_cumulative_poisson
from math import log
from VariantCounter import VariantCounter, load_cstag
import subprocess
import gc

min_rlen = 9000
reflen = 13332
output_file_prefix = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/"
graph_dir = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/"

skip_counting = True

print("PROTECTION: Do you really want to run this script? again?(y/n) (It will takes 510 min to run this script))")
if input() != "y":
    print("exit")
    sys.exit()




# load the Bcell ID list
Bcell_IDlist_file = "/nfs/data05/otgk/rDNA/Bcell_IDlist.txt"
Bcell_IDlist = []
with open(Bcell_IDlist_file, "r") as f:
    for line in f:
        Bcell_IDlist.append(line.strip())
# print(Bcell_IDlist)


# for test:
# Bcell_IDlist = Bcell_IDlist



# this process takes 1 sec for each step. 
    # this is slow but acceptable, so if affords to be slow, necessary to parallelize this process.
# Bcell_IDlist = ["B357"]
if not skip_counting:
    for Bcell_ID in Bcell_IDlist:
        print("Bcell_ID: ", Bcell_ID)
        paf_file = f"/nfs/data05/otgk/rDNA/mapping/PSCA_mapping/ref=consensusmapping/mapping_result/{Bcell_ID}toCons.paf"
        # output_file_prefix = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/result/{Bcell_ID}/internal_{Bcell_ID}toCons"


        # ------------------------------------------- input processing -------------------------------------------
        # load the paf file and use only the columns of "cstag"(), "Target start on original strand (0-based)" and "Target end on original strand (0-based)", "Query start on original strand (0-based)", "Query end on original strand (0-based)", "Number of matching bases in the mapping", "Number bases in the mapping"
            # https://github.com/lh3/miniasm/blob/master/PAF.md
        # print(pd.read_csv(paf_file, sep='\s+', header=None, usecols=[7,8,22,16] ).tail(50))
        # paf_cs_df_ = pd.read_csv(paf_file, sep='\t', header=None, usecols=[7, 8, 22], names=[ "ref_start", "ref_end", "long_cstag"])
        paf_cs_df_ = pd.read_csv(paf_file, sep='\s+', header=None, usecols=[7, 8, 16, 22], names=[ "ref_start", "ref_end", "tp", "long_cstag"])
        paf_cs_df_ = paf_cs_df_[paf_cs_df_["tp"] == "tp:A:P"] # primay aligment
        # print(paf_cs_df_.describe())
        # print(paf_cs_df_[index = 3353])



        paf_cs_df_["rlen"] = paf_cs_df_["ref_end"] - paf_cs_df_["ref_start"]
        # paf_cs_df_["short_cstag"]: for each long_cstag, apply "cstag.shorten(long_cstag)"

        paf_cs_df_["short_cstag"] = paf_cs_df_["long_cstag"].apply(cstag.shorten)

        # extract the entries with rlen >= min_rlen
        paf_cs_df = paf_cs_df_[paf_cs_df_["rlen"] >= min_rlen]
        paf_cs_df = paf_cs_df.drop(columns=["long_cstag", "rlen", "tp"])
        del paf_cs_df_
        gc.collect()


        # ------------------------------------------- caluculate Variant count for each site and DP -------------------------------------------
        Bcell_variant_counter = VariantCounter(reflen)
        # for each site, judge if there is a variant and accumulate the read depth
        for ref_start, ref_end, short_cstag in paf_cs_df.values:
            # ref_start, ref_end: 0-based integer

            Bcell_variant_counter.add_DP(ref_start, ref_end)

            # read the parsed_cstag(=splitted_cstag) from the short_cstag
            Bcell_variant_counter = load_cstag(cstag.split(short_cstag), ref_start, Bcell_variant_counter)
        
        # export the Bcell_variant_counter.
        Bcell_variant_counter.to_csv(f"{output_file_prefix}result/{Bcell_ID}/{Bcell_ID}2Cons.variant.csv", Bcell_ID)


        # ------------------------------------------- calculate the Variant Frequency and p-value for each site -------------------------------------------
        # calculate the Variant Frequenc (for each ) in the Bcell we are looking at
        Bcell_variant_counter.calc_freq()


    # memory free:
    del paf_cs_df, Bcell_variant_counter
    gc.collect()



samples_vs_SNV_count = pd.DataFrame(columns=["POS"])
samples_vs_INS_count = pd.DataFrame(columns=["POS"])
samples_vs_DEL_count = pd.DataFrame(columns=["POS"])
samples_vs_DP = pd.DataFrame(columns=["POS"])

samples_vs_SNV_count["POS"] = list(range(reflen))
samples_vs_INS_count["POS"] = list(range(reflen))
samples_vs_DEL_count["POS"] = list(range(reflen))
samples_vs_DP["POS"] = list(range(reflen))

for Bcell_ID in Bcell_IDlist:
    samples_vs_SNV_count[f"{Bcell_ID}"] = pd.read_csv(f"{output_file_prefix}result/{Bcell_ID}/{Bcell_ID}2Cons.variant.csv", sep=',', header=0, )["SNV_COUNT"]
    samples_vs_INS_count[f"{Bcell_ID}"] = pd.read_csv(f"{output_file_prefix}result/{Bcell_ID}/{Bcell_ID}2Cons.variant.csv", sep=',', header=0,)["INS_COUNT"]
    samples_vs_DEL_count[f"{Bcell_ID}"] = pd.read_csv(f"{output_file_prefix}result/{Bcell_ID}/{Bcell_ID}2Cons.variant.csv", sep=',', header=0,)["DEL_COUNT"]
    samples_vs_DP[f"{Bcell_ID}"] = pd.read_csv(f"{output_file_prefix}result/{Bcell_ID}/{Bcell_ID}2Cons.variant.csv", sep=',', header=0,)["DP"]

samples_vs_SNV_freq = samples_vs_SNV_count.drop(columns=["POS"]) / samples_vs_DP.drop(columns=["POS"])
samples_vs_INS_freq = samples_vs_INS_count.drop(columns=["POS"]) / samples_vs_DP.drop(columns=["POS"])
samples_vs_DEL_freq = samples_vs_DEL_count.drop(columns=["POS"]) / samples_vs_DP.drop(columns=["POS"])

samples_vs_SNV_freq["POS"] = list(range(reflen))
samples_vs_INS_freq["POS"] = list(range(reflen))
samples_vs_DEL_freq["POS"] = list(range(reflen))

print(samples_vs_SNV_freq.head())
print("shape of samples_vs_SNV_freq: ", samples_vs_SNV_freq.shape)


# ------------------------------------------- export -------------------------------------------
# export the dataframe
# print(f"writing df_posfreq_vs_Bcell to {output_file_prefix}vaf.csv")
# df_posfreq_vs_Bcell.to_csv(f"{output_file_prefix}vaf.csv", sep=',', index=False) #--> this is memory consuming, have to be improved.

samples_vs_SNV_count.to_csv(f"{output_file_prefix}POSvsBcellSNVcount.csv", sep=',', index=False)
samples_vs_INS_count.to_csv(f"{output_file_prefix}POSvsBcellINScount.csv", sep=',', index=False)
samples_vs_DEL_count.to_csv(f"{output_file_prefix}POSvsBcellDELcount.csv", sep=',', index=False)
samples_vs_DP.to_csv(f"{output_file_prefix}POSvsBcellDP.csv", sep=',', index=False)

samples_vs_SNV_freq.to_csv(f"{output_file_prefix}POSvsBcellSNVfreq.csv", sep=',', index=False)
samples_vs_INS_freq.to_csv(f"{output_file_prefix}POSvsBcellINSfreq.csv", sep=',', index=False)
samples_vs_DEL_freq.to_csv(f"{output_file_prefix}POSvsBcellDELfreq.csv", sep=',', index=False)



# ------------------------------------------- calculate the mean Variant Frequency for each site -------------------------------------------
# calculate the mean Variant Frequency for each site
# mean_variant_freq_snv_j (j : position) = mean(VAF_SNV_ji) (i : Bcell)



# mean_vaf_SNV = [sum(df_posfreq_vs_Bcell[f"SNV({Bcell_ID})"][j] for Bcell_ID in Bcell_IDlist) / len(Bcell_IDlist) for j in range(reflen)]
# mean_vaf_INS = [sum(df_posfreq_vs_Bcell[f"INS({Bcell_ID})"][j] for Bcell_ID in Bcell_IDlist) / len(Bcell_IDlist) for j in range(reflen)]
# mean_vaf_DEL = [sum(df_posfreq_vs_Bcell[f"DEL({Bcell_ID})"][j] for Bcell_ID in Bcell_IDlist) / len(Bcell_IDlist) for j in range(reflen)]
# using numpy, instead of above code:
# mean vaf: (sum_i (VAF_ij)) / len(Bcell_IDlist)
import numpy as np

mean_vaf_SNV = np.mean(samples_vs_SNV_freq.drop(columns=["POS"]).values, axis=1)
mean_vaf_INS = np.mean(samples_vs_INS_freq.drop(columns=["POS"]).values, axis=1)
mean_vaf_DEL = np.mean(samples_vs_DEL_freq.drop(columns=["POS"]).values, axis=1)

del samples_vs_SNV_freq, samples_vs_INS_freq, samples_vs_DEL_freq
gc.collect()




# ------------------------------------------- calculate the log prob for each site for each Bcell -------------------------------------------

# prob_df = pd.DataFrame(columns=["POS", "log_P_SNV(Bcell_ID)", "log_P_INS(Bcell_ID)", "log_P_DEL(Bcell_ID)" for Bcell_ID in Bcell_IDlist])
prob_df = pd.DataFrame(columns=["POS"])
prob_df["POS"] = list(range(reflen))

# calculate the log p-value for each site for each Bcell
for Bcell in Bcell_IDlist:
    # poisson.log probability density function not probability mass function
    # prob_df[f"log_P_{variant_type}({Bcell})"] = [log_continuous_poisson(df_posfreq_vs_Bcell[f"{variant_type}({Bcell})"][j], mean_vaf[j]) for j in range(reflen)]
    # log_p_value = log_continious_poisson(VariantCount_ij, mean_j * DP_ij): not VariantFrequency but VariantCount
    # prob_df[f"log_P_{variant_type}({Bcell})"] = [log_continuous_poisson(df_posfreq_vs_Bcell[f"{variant_type}({Bcell})"][j] * DP_df[f"DP({Bcell})"][j], mean_vaf[j] * DP_df[f"DP({Bcell})"][j]) for j in range(reflen)]
 
    for variant_type in ["SNV", "INS", "DEL"]:
        if variant_type == "SNV":
            mean_vaf = mean_vaf_SNV
            # prob_df[f"log_P_SNV({Bcell})"] = [log_continuous_poisson(samples_vs_SNV_count[f"{Bcell}"][j], mean_vaf[j] * samples_vs_DP[f"{Bcell}"][j]) for j in range(reflen)] 
            # above is wrong, this code is about probabily of observing VariantCount_ij, not p-value of VariantFrequency_ij
            prob_df[f"log_P_SNV({Bcell})"] = [log_continuous_cumulative_poisson(samples_vs_SNV_count[f"{Bcell}"][j], mean_vaf[j] * samples_vs_DP[f"{Bcell}"][j]) for j in range(reflen)] 
        elif variant_type == "INS":
            mean_vaf = mean_vaf_INS
            prob_df[f"log_P_INS({Bcell})"] = [log_continuous_cumulative_poisson(samples_vs_INS_count[f"{Bcell}"][j], mean_vaf[j] * samples_vs_DP[f"{Bcell}"][j]) for j in range(reflen)]
        elif variant_type == "DEL":
            mean_vaf = mean_vaf_DEL
            prob_df[f"log_P_DEL({Bcell})"] = [log_continuous_cumulative_poisson(samples_vs_DEL_count[f"{Bcell}"][j], mean_vaf[j] * samples_vs_DP[f"{Bcell}"][j]) for j in range(reflen)]

  
# ------------------------------------------- export -------------------------------------------
print(f"writing prob_df to {output_file_prefix}internal_logpvalue.csv")

# remove the rows with all NA except POS
# prob_df = prob_df.dropna(how="all", subset=prob_df.columns[1:])
prob_df.to_csv(f"{output_file_prefix}internal_logpvalue.csv", sep=',', index=False, na_rep="")

# POS vs SNV log p-value
output_snv = prob_df[["POS"] + [f"log_P_SNV({Bcell})" for Bcell in Bcell_IDlist]]
# output_snv = output_snv.dropna(how="all", subset=output_snv.columns[1:])
output_snv.to_csv(f"{output_file_prefix}internal_logpvalue_SNV.csv", sep=',', index=False, na_rep="")

# POS vs INS log p-value
output_ins = prob_df[["POS"] + [f"log_P_INS({Bcell})" for Bcell in Bcell_IDlist]]
# output_ins = output_ins.dropna(how="all", subset=output_ins.columns[1:])
output_ins.to_csv(f"{output_file_prefix}internal_logpvalue_INS.csv", sep=',', index=False, na_rep="")
# POS vs DEL log p-value
output_del = prob_df[["POS"] + [f"log_P_DEL({Bcell})" for Bcell in Bcell_IDlist]]
# output_del = output_del.dropna(how="all", subset=output_del.columns[1:])
output_del.to_csv(f"{output_file_prefix}internal_logpvalue_DEL.csv", sep=',', index=False, na_rep="")

