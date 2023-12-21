import matplotlib.pyplot as plt
import pandas as pd

from math import log

# pvalue_df.to_csv(f"{output_file_prefix}internal_logpvalue.csv", sep=',', index=False, na_rep="")

# # POS vs SNV log p-value
# output_snv = pvalue_df[["POS"] + [f"log_P_SNV({HGNA})" for HGNA in HGNA_IDlist]]
# # output_snv = output_snv.dropna(how="all", subset=output_snv.columns[1:])
# output_snv.to_csv(f"{output_file_prefix}internal_logpvalue_SNV.csv", sep=',', index=False, na_rep="")

# # POS vs INS log p-value
# output_ins = pvalue_df[["POS"] + [f"log_P_INS({HGNA})" for HGNA in HGNA_IDlist]]
# # output_ins = output_ins.dropna(how="all", subset=output_ins.columns[1:])
# output_ins.to_csv(f"{output_file_prefix}internal_logpvalue_INS.csv", sep=',', index=False, na_rep="")
# # POS vs DEL log p-value
# output_del = pvalue_df[["POS"] + [f"log_P_DEL({HGNA})" for HGNA in HGNA_IDlist]]
# # output_del = output_del.dropna(how="all", subset=output_del.columns[1:])
# output_del.to_csv(f"{output_file_prefix}internal_logpvalue_DEL.csv", sep=',', index=False, na_rep="")


min_rlen = 9000
reflen = 13332
output_file_prefix = f"/large/otgk/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/"
graph_dir = f"/large/otgk/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/graph/"

# load file
output_snv = pd.read_csv(f"{output_file_prefix}internal_logpvalue_SNV.csv", sep=',', na_values=[""])
output_ins = pd.read_csv(f"{output_file_prefix}internal_logpvalue_INS.csv", sep=',', na_values=[""])
output_del = pd.read_csv(f"{output_file_prefix}internal_logpvalue_DEL.csv", sep=',', na_values=[""])

# drop the row with all NA
output_snv = output_snv.dropna(how="all", subset=output_snv.columns[1:])
output_ins = output_ins.dropna(how="all", subset=output_ins.columns[1:])
output_del = output_del.dropna(how="all", subset=output_del.columns[1:])
# print(output_snv.head())

print(output_snv.head())


# load the HGNA ID list
HGNA_IDlist_file = "/large/otgk/rDNA/HGNA_IDlist.txt"
HGNA_IDlist = []
with open(HGNA_IDlist_file, "r") as f:
    for line in f:
        HGNA_IDlist.append(line.strip())


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
alpha = 0.05
fig = plt.figure(figsize=(22, 12))
fig.suptitle(f"min(log_p_value) vs position:\n num of HGNAs = {len(HGNA_IDlist)}\n alpha = {alpha} ", fontsize=20)

log_alpha_limit = \
 log(alpha) - log(len(HGNA_IDlist)) - log((len(output_snv["POS"]) + len(output_ins["POS"]) + len(output_del["POS"])))

ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("SNV")
ax1.set_xlabel("position")
ax1.set_ylabel("min(log_p_value)")
# ax1.set_ylim(-10, 0)
ax1.scatter(output_snv["POS"], [min([output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist]) for j in output_snv["POS"]])
ax1.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_snv["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(snv_variant))")
ax1.legend(loc = "upper left")

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("INS")
ax2.set_xlabel("position")
ax2.set_ylabel("min(log_p_value)")
# ax2.set_ylim(-10, 0)
ax2.scatter(output_ins["POS"], [min([output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist]) for j in output_ins["POS"]])
ax2.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_ins["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(ins_variant))")
ax2.legend(loc = "upper left")

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("DEL")
ax3.set_xlabel("position")
ax3.set_ylabel("min(log_p_value)")
# ax3.set_ylim(-10, 0)
ax3.scatter(output_del["POS"], [min([output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist]) for j in output_del["POS"]])
ax3.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_del["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(del_variant))")
ax3.legend(loc = "upper left")


ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("SNV, INS, DEL", fontsize=20)
ax4.set_xlabel("position")
ax4.set_ylabel("min(log_p_value)")
# ax4.set_ylim(-10, 0)
ax4.scatter(output_snv["POS"], [min([output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist]) for j in output_snv["POS"]], label="SNV", alpha = 0.7)
ax4.scatter(output_ins["POS"], [min([output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist]) for j in output_ins["POS"]], label="INS", alpha = 0.7)
ax4.scatter(output_del["POS"], [min([output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist]) for j in output_del["POS"]], label="DEL", alpha = 0.7)
ax4.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log((len(output_snv["POS"]) + len(output_ins["POS"]) + len(output_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(all_variant))")
ax4.legend(loc = "upper left")


# save
fig.savefig(f"{graph_dir}internal_minlogpvalue.png")
print(f"{graph_dir}internal_minlogpvalue.png is saved")


# the same plot but the y should be median(log_p_value) instead of min(log_p_value)

fig = plt.figure(figsize=(22, 12))
fig.suptitle(f"median(log_p_value) vs position:\n num of HGNAs = {len(HGNA_IDlist)}\n alpha = {alpha} ")

ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("SNV")
ax1.set_xlabel("position")
ax1.set_ylabel("median(log_p_value)")
ax1.scatter(output_snv["POS"], [sorted([output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 2] for j in output_snv["POS"]])
ax1.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_snv["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(snv_variant))")
ax1.legend(loc = "upper left")

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("INS")
ax2.set_xlabel("position")
ax2.set_ylabel("median(log_p_value)")
ax2.scatter(output_ins["POS"], [sorted([output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 2] for j in output_ins["POS"]])
ax2.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_ins["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(ins_variant))")
ax2.legend(loc = "upper left")

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("DEL")
ax3.set_xlabel("position")
ax3.set_ylabel("median(log_p_value)")
ax3.scatter(output_del["POS"], [sorted([output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 2] for j in output_del["POS"]])
ax3.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_del["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(del_variant))")
ax3.legend(loc = "upper left")

ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("SNV, INS, DEL")
ax4.set_xlabel("position")
ax4.set_ylabel("median(log_p_value)")
ax4.scatter(output_snv["POS"], [sorted([output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 2] for j in output_snv["POS"]], label="SNV", alpha = 0.7)
ax4.scatter(output_ins["POS"], [sorted([output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 2] for j in output_ins["POS"]], label="INS", alpha = 0.7)
ax4.scatter(output_del["POS"], [sorted([output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 2] for j in output_del["POS"]], label="DEL", alpha = 0.7)
ax4.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log((len(output_snv["POS"]) + len(output_ins["POS"]) + len(output_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(all_variant))")
ax4.legend(loc = "upper left")

# save
fig.savefig(f"{graph_dir}internal_medianlogpvalue.png")
print(f"{graph_dir}internal_medianlogpvalue.png is saved")




# the same plot but the y should be 1st quartile.(log_p_value) instead of min(log_p_value)

fig = plt.figure(figsize=(22, 12))
fig.suptitle(f"1st quartile(log_p_value) vs position:\n num of HGNAs = {len(HGNA_IDlist)}\n alpha = {alpha} ", fontsize=20)

ax1 = fig.add_subplot(2, 2, 1)
ax1.set_title("SNV")
ax1.set_xlabel("position")
ax1.set_ylabel("1st quartile(log_p_value)")
ax1.scatter(output_snv["POS"], [sorted([output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_snv["POS"]])
ax1.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_snv["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(snv_variant))")
ax1.legend(loc = "upper left")

ax2 = fig.add_subplot(2, 2, 2)
ax2.set_title("INS")
ax2.set_xlabel("position")
ax2.set_ylabel("1st quartile(log_p_value)")
ax2.scatter(output_ins["POS"], [sorted([output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_ins["POS"]])
ax2.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_ins["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(ins_variant))")
ax2.legend(loc = "upper left")

ax3 = fig.add_subplot(2, 2, 3)
ax3.set_title("DEL")
ax3.set_xlabel("position")
ax3.set_ylabel("1st quartile(log_p_value)")
ax3.scatter(output_del["POS"], [sorted([output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_del["POS"]])
ax3.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log(len(output_del["POS"])), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(del_variant))")
ax3.legend(loc = "upper left")

ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title("SNV, INS, DEL")
ax4.set_xlabel("position")
ax4.set_ylabel("1st quartile(log_p_value)")
ax4.scatter(output_snv["POS"], [sorted([output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_snv["POS"]], label="SNV", alpha = 0.7)
ax4.scatter(output_ins["POS"], [sorted([output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_ins["POS"]], label="INS", alpha = 0.7)
ax4.scatter(output_del["POS"], [sorted([output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_del["POS"]], label="DEL", alpha = 0.7)
ax4.axhline(y=log(alpha) - log(len(HGNA_IDlist)) - log((len(output_snv["POS"]) + len(output_ins["POS"]) + len(output_del["POS"]))), color='r', linestyle='-', label="log(0.05 / len(HGNA_list) / len(all_variant))")
ax4.legend(loc = "upper left")

# save
fig.savefig(f"{graph_dir}internal_1stquartilelogpvalue.png")
print(f"{graph_dir}internal_1stquartilelogpvalue.png is saved")

