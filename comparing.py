# p-value の 1st quartile を PSCA と HPRC で比較する

import matplotlib.pyplot as plt
import pandas as pd

from math import log

# pvalue_df.to_csv(f"{PSCA_dir}internal_logpvalue.csv", sep=',', index=False, na_rep="")

# # POS vs SNV log p-value
# PSCA_snv = pvalue_df[["POS"] + [f"log_P_SNV({HGNA})" for HGNA in HGNA_IDlist]]
# # PSCA_snv = PSCA_snv.dropna(how="all", subset=PSCA_snv.columns[1:])
# PSCA_snv.to_csv(f"{PSCA_dir}internal_logpvalue_SNV.csv", sep=',', index=False, na_rep="")

# # POS vs INS log p-value
# PSCA_ins = pvalue_df[["POS"] + [f"log_P_INS({HGNA})" for HGNA in HGNA_IDlist]]
# # PSCA_ins = PSCA_ins.dropna(how="all", subset=PSCA_ins.columns[1:])
# PSCA_ins.to_csv(f"{PSCA_dir}internal_logpvalue_INS.csv", sep=',', index=False, na_rep="")
# # POS vs DEL log p-value
# PSCA_del = pvalue_df[["POS"] + [f"log_P_DEL({HGNA})" for HGNA in HGNA_IDlist]]
# # PSCA_del = PSCA_del.dropna(how="all", subset=PSCA_del.columns[1:])
# PSCA_del.to_csv(f"{PSCA_dir}internal_logpvalue_DEL.csv", sep=',', index=False, na_rep="")


min_rlen = 9000
reflen = 13332
alpha=0.05

PSCA_dir = f"/large/otgk/rDNA/rDNAsimpleVariantcalling/PSCA_variants/consensus/adhoc/"
HPRC_dir = f"/large/otgk/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/"
graph_dir = f"/large/otgk/rDNA/rDNAsimpleVariantcalling/graph/"

# load file
PSCA_snv = pd.read_csv(f"{PSCA_dir}internal_logpvalue_SNV.csv", sep=',', na_values=[""], header=0)
PSCA_ins = pd.read_csv(f"{PSCA_dir}internal_logpvalue_INS.csv", sep=',', na_values=[""], header=0)
PSCA_del = pd.read_csv(f"{PSCA_dir}internal_logpvalue_DEL.csv", sep=',', na_values=[""], header=0)

# drop the row with all NA
PSCA_snv = PSCA_snv.dropna(how="all", subset=PSCA_snv.columns[1:])
PSCA_ins = PSCA_ins.dropna(how="all", subset=PSCA_ins.columns[1:])
PSCA_del = PSCA_del.dropna(how="all", subset=PSCA_del.columns[1:])
# print(PSCA_snv.head())


# load the HGNA ID list
HPRC_snv = pd.read_csv(f"{HPRC_dir}internal_logpvalue_SNV.csv", sep=',', na_values=[""], header=0)
HPRC_ins = pd.read_csv(f"{HPRC_dir}internal_logpvalue_INS.csv", sep=',', na_values=[""], header=0)
HPRC_del = pd.read_csv(f"{HPRC_dir}internal_logpvalue_DEL.csv", sep=',', na_values=[""], header=0)

# drop the row with all NA
HPRC_snv = HPRC_snv.dropna(how="all", subset=HPRC_snv.columns[1:])
HPRC_ins = HPRC_ins.dropna(how="all", subset=HPRC_ins.columns[1:])
HPRC_del = HPRC_del.dropna(how="all", subset=HPRC_del.columns[1:])


print(PSCA_snv.head())
print(HPRC_snv.head())


# load the HGNA ID list
HGNA_IDlist_file = "/large/otgk/rDNA/HGNA_IDlist.txt"
HGNA_IDlist = []
with open(HGNA_IDlist_file, "r") as f:
    for line in f:
        HGNA_IDlist.append(line.strip())
PSCA_IDlist_file = "/large/otgk/rDNA/Bcell_IDlist.txt"
PSCA_IDlist = []
with open(PSCA_IDlist_file, "r") as f:
    for line in f:
        PSCA_IDlist.append(line.strip())

# ------------------------------------------- plot -------------------------------------------
# plot the pvalue
# for each position j, plot the min_{HGNAs i}(log_value_ij) vs j

# use subplot and plot 
    # 1. snv: pos vs min_{HGNAs i}(log_value_ij)
    # 2. ins: pos vs min_{HGNAs i}(log_value_ij)
    # 3. del: pos vs min_{HGNAs i}(log_value_ij)
    # 4. 3 plots in the same figure (different colors)
    # * all figure have the y = log(0.05 / len(HGNA_list) / len(variant))
    # in total, there are 4 figures
        
# PSCA, 

fig = plt.figure(figsize=(22, 12))
fig.suptitle(f"1st quartile(log_p_value) vs position:\n num of HGNAs = {len(HGNA_IDlist)}, num of PSCA = {len(PSCA_IDlist)}\n alpha = {alpha} ", fontsize=20)

# PSCA
ax4 = fig.add_subplot(2, 1, 1)
ax4.set_title("SNV, INS, DEL")
ax4.set_xlabel("position")
ax4.set_ylabel("1st quartile(log_p_value)")
ax4.scatter(PSCA_snv["POS"], [sorted([PSCA_snv[f"log_P_SNV({PSCA})"][j] for PSCA in PSCA_IDlist])[len(PSCA_IDlist) // 4] for j in PSCA_snv["POS"]], label="BcellSNV", alpha = 0.7)
ax4.scatter(PSCA_ins["POS"], [sorted([PSCA_ins[f"log_P_INS({PSCA})"][j] for PSCA in PSCA_IDlist])[len(PSCA_IDlist) // 4] for j in PSCA_ins["POS"]], label="BcellINS", alpha = 0.7)
ax4.scatter(PSCA_del["POS"], [sorted([PSCA_del[f"log_P_DEL({PSCA})"][j] for PSCA in PSCA_IDlist])[len(PSCA_IDlist) // 4] for j in PSCA_del["POS"]], label="BcellDEL", alpha = 0.7)
ax4.axhline(y=log(alpha) - log(len(PSCA_IDlist)) - log((len(PSCA_snv["POS"]) + len(PSCA_ins["POS"]) + len(PSCA_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(all_variant))")
ax4.legend(loc = "upper left")


# HPRC (HGNA)
ax4 = fig.add_subplot(2, 1, 2)
ax4.set_title("SNV, INS, DEL")
ax4.set_xlabel("position")
ax4.set_ylabel("1st quartile(log_p_value)")
ax4.scatter(HPRC_snv["POS"], [sorted([HPRC_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in HPRC_snv["POS"]], label="hprcSNV", alpha = 0.7)
ax4.scatter(HPRC_ins["POS"], [sorted([HPRC_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in HPRC_ins["POS"]], label="hprcINS", alpha = 0.7)
ax4.scatter(HPRC_del["POS"], [sorted([HPRC_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in HPRC_del["POS"]], label="hprcDEL", alpha = 0.7)
ax4.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log((len(HPRC_snv["POS"]) + len(PSCA_ins["POS"]) + len(PSCA_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(all_variant))")
ax4.legend(loc = "upper left")




# save
fig.savefig(f"{graph_dir}internal_1q.png")
print(f"{graph_dir}internal_1stquartilelogpvalue.png is saved")




# 重ね合わせる

plt.figure(figsize=(22, 12))
plt.title(f"1st quartile(log_p_value) vs position:\n num of HGNAs = {len(HGNA_IDlist)}\n alpha = {alpha} ", fontsize=20)

# PSCA
plt.scatter(PSCA_snv["POS"], [sorted([PSCA_snv[f"log_P_SNV({PSCA})"][j] for PSCA in PSCA_IDlist])[len(PSCA_IDlist) // 4] for j in PSCA_snv["POS"]], label="BcellSNV", alpha = 0.7, color="blue")
plt.scatter(PSCA_ins["POS"], [sorted([PSCA_ins[f"log_P_INS({PSCA})"][j] for PSCA in PSCA_IDlist])[len(PSCA_IDlist) // 4] for j in PSCA_ins["POS"]], label="BcellINS", alpha = 0.7, color="orange")
plt.scatter(PSCA_del["POS"], [sorted([PSCA_del[f"log_P_DEL({PSCA})"][j] for PSCA in PSCA_IDlist])[len(PSCA_IDlist) // 4] for j in PSCA_del["POS"]], label="BcellDEL", alpha = 0.7, color="green")

# HPRC
plt.scatter(HPRC_snv["POS"], [sorted([HPRC_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in HPRC_snv["POS"]], label="hprcSNV", alpha = 0.7, color="skyblue")
plt.scatter(HPRC_ins["POS"], [sorted([HPRC_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in HPRC_ins["POS"]], label="hprcINS", alpha = 0.7, color="coral")
plt.scatter(HPRC_del["POS"], [sorted([HPRC_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in HPRC_del["POS"]], label="hprcDEL", alpha = 0.7, color="limegreen")

# log(0.05 / len(HGNA_list) / len(variant))

plt.axhline(y=log(alpha) - log(len(PSCA_IDlist)) - log((len(PSCA_snv["POS"]) + len(PSCA_ins["POS"]) + len(PSCA_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(PSCA_list) / len(all_variant))")
plt.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log((len(HPRC_snv["POS"]) + len(HPRC_ins["POS"]) + len(HPRC_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(all_variant))")

plt.legend(loc = "upper left")

# save
plt.savefig(f"{graph_dir}internal_1q_overlayed.png")