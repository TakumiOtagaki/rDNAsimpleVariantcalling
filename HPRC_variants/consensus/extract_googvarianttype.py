import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

min_rlen = 9000
reflen = 13332

HGNA_IDlist_file = "/large/otgk/rDNA/HG_IDlist.txt"
HGNA_IDlist = []
with open(HGNA_IDlist_file, "r") as f:
    for line in f:
        HGNA_IDlist.append(line.strip())

# output_file_prefix = f"/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/"
output_file_prefix = f"/large/otgk/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/"

# load file
output_snv = pd.read_csv(f"{output_file_prefix}internal_logpvalue_SNV.csv", sep=',', na_values=[""])
output_ins = pd.read_csv(f"{output_file_prefix}internal_logpvalue_INS.csv", sep=',', na_values=[""])
output_del = pd.read_csv(f"{output_file_prefix}internal_logpvalue_DEL.csv", sep=',', na_values=[""])

# drop the row with all NA
output_snv = output_snv.dropna(how="all", subset=output_snv.columns[1:])
output_ins = output_ins.dropna(how="all", subset=output_ins.columns[1:])
output_del = output_del.dropna(how="all", subset=output_del.columns[1:])
# print(output_snv.head())

# ------------------- SNV -------------------
# 1st quartile(log_p_value) vs position

# first_q_snv = [sorted([-1 * output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_snv["POS"]]
# first_q_ins = [sorted([-1 * output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_ins["POS"]]
# first_q_del = [sorted([-1 * output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4] for j in output_del["POS"]]

first_q_dict = dict()
for j in output_snv["POS"]:
    first_q_dict[(j, "SNV")] = sorted([-1 * output_snv[f"log_P_SNV({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4]
for j in output_ins["POS"]:
    first_q_dict[(j, "INS")] = sorted([-1 * output_ins[f"log_P_INS({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4]
for j in output_del["POS"]:
    first_q_dict[(j, "DEL")] = sorted([-1 * output_del[f"log_P_DEL({HGNA})"][j] for HGNA in HGNA_IDlist])[len(HGNA_IDlist) // 4]

# snv, ins, del を合わせてから、そのうち -log の値が最も大きいようなものから順番に全体の 10 % (13332 * 3 * 0.1 = 4000) を取り出す。
# 出力するのは、SNV or INS or DEL と、position j と -log の値

first_q_dict_sorted = sorted(first_q_dict.items(), key=lambda x: x[1], reverse=True)
first_q_dict_sorted = first_q_dict_sorted[:5000]


# ファイルに書き出す
with open(f"{output_file_prefix}first_q_dict_sorted.txt", "w") as f:
    for x in first_q_dict_sorted:
        f.write(f"{x[0][0]},{x[0][1]},{x[1]}\n")
