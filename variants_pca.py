# variant calling を行った後の、各サイト（j）、各バリアント種類（SNV, Ins, Del） の全ての組み合わせについて、
# 各サンプル i はそのサイト j & バリアント種類 k のバリアントをサポートするようなリードの頻度を 0 ~ 1 で正規化したものを計算してある。
# それらの値を利用して、ベクトルを作成し、Japanese と Non-Japanese の比較を行う。

import pandas as pd
import numpy as np
import os

# workdir = $workdir : 環境変数
workdir = os.environ['workdir']
print(workdir)

# path of data
path_dir = { ("Japanese", "SNV") : workdir + "/rDNA/rDNAsimpleVariantcalling/PSCA_variants/consensus/adhoc/POSvsBcellSNVfreq.csv",
            ("Japanese", "INS") : workdir + "/rDNA/rDNAsimpleVariantcalling/PSCA_variants/consensus/adhoc/POSvsBcellINSfreq.csv",
            ("Japanese", "DEL") : workdir + "/rDNA/rDNAsimpleVariantcalling/PSCA_variants/consensus/adhoc/POSvsBcellDELfreq.csv",
            ("Non-Japanese", "SNV") : workdir + "/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/POSvsHGNASNVfreq.csv",
            ("Non-Japanese", "INS") : workdir + "/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/POSvsHGNAINSfreq.csv",
            ("Non-Japanese", "DEL") : workdir + "/rDNA/rDNAsimpleVariantcalling/HPRC_variants/consensus/POSvsHGNADELfreq.csv" 
            }

# read data
df = {}

for key in path_dir.keys():
    # pos, sample1, sample2, ...,
    # 1, 0.2, 0.8, ...,
    # 2, 0.3, 0.7, ...,
    # ...

    df[key] = pd.read_csv(path_dir[key], header=0)
    df[key].set_index("POS", inplace=True)

N = len(df[("Japanese", "SNV")]) # 13332
Japanese_ids = df[("Japanese", "SNV")].columns
NonJapanese_ids = df[("Non-Japanese", "SNV")].columns

print(N, Japanese_ids, NonJapanese_ids)
# make vector: 各サンプル について、各サイト（j）、各バリアント種類（SNV, Ins, Del） の全ての組み合わせについて、頻度のベクトルを作成。
    # 長さは N(=13332) * 3 となる
# Japanese
Japanese_vecs = np.zeros((len(Japanese_ids), N*3))
for i in range(Japanese_ids):
    Japanese_vecs[i] = np.concatenate([df[("Japanese", "SNV")].iloc[:,i], df[("Japanese", "INS")].iloc[:,i], df[("Japanese", "DEL")].iloc[:,i]])

# Non-Japanese
NonJapanese_vecs = np.zeros((len(NonJapanese_ids), N*3))
for i in range(NonJapanese_ids):
    NonJapanese_vecs[i] = np.concatenate([df[("Non-Japanese", "SNV")].iloc[:,i], df[("Non-Japanese", "INS")].iloc[:,i], df[("Non-Japanese", "DEL")].iloc[:,i]])

# PCA: Japanese and Nonjapanese
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
pca.fit(Japanese_vecs + NonJapanese_vecs)
# plot
import matplotlib.pyplot as plt
plt.scatter(Japanese_pca[:,0], Japanese_pca[:,1], label="Japanese")
plt.scatter(NonJapanese_pca[:,0], NonJapanese_pca[:,1], label="Non-Japanese")


# Non-Japanese


