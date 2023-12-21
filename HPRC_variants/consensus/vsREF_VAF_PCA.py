import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

# Read the merged VAF file
merged_VAF_file = '/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/merged_VAF.csv'
df_merged = pd.read_csv(merged_VAF_file, sep=',', header=0)
# df_merged:
# POS_Variant_type,B003_VAF,B005_VAF,B007_VAF,B009_VAF,B012_VAF,B022_VAF,B025_VAF,B026_VAF,B027_VAF, ...

# each Bcell VAF, do PCA anaysis and plot the first two PCs

# ------------------------------------------------- setting -------------------------------------------------
# set the marker size
marker_size = 3

# set the font size
plt.rcParams.update({'font.size': 12})

# calculate PCA for SNV+INS+DEL SNV, INS, DEL, SNV+INS, SNV+DEL, INS+DEL
fig = plt.figure(figsize=(15, 20))




# ------------------------------------------------- PCA for SNV+INS+DEL -------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)

pca.fit(df_merged[df_merged.columns[1:]])
X = pca.transform(df_merged[df_merged.columns[1:]])
# print(X)
# print(X.shape)
# print(pca.explained_variance_ratio_)
# print(pca.singular_values_)
# print(pca.components_)


# --------------- plot ---------------
plt.figure(figsize=(8, 8))


plt.scatter(X[:, 0], X[:, 1])
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA plot of Bcell VAF for SNV+INS+DEL')
plt.savefig('/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/graph/PCA_plot_all.png')

# 寄与率
plt.figure(figsize=(8, 8))
cumulative_explained_variance_ratio = np.cumsum([0] + list(pca.explained_variance_ratio_))
plt.plot(cumulative_explained_variance_ratio)
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance')
plt.savefig('/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/graph/PCA_plot_all_explained_variance_ratio.png')


# --------------------------------------------------- PCA for SNV ---------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)
# use df_merged["POS_Variant_type"].str.contains('SNV') to select SNV
df_SNV = df_merged[df_merged["POS_Variant_type"].str.contains('SNV')]
pca.fit(df_SNV[df_SNV.columns[1:]])
X = pca.transform(df_SNV[df_SNV.columns[1:]])
# print(X)
# print(X.shape)

# --------------- plot ---------------
ax1 = fig.add_subplot(3, 2, 1)
# plt.figure(figsize=(8, 8))
# plt.scatter(X[:, 0], X[:, 1], s=marker_size)
# plt.xlabel('PC1')
# plt.ylabel('PC2')
# plt.title('PCA plot of Bcell VAF for SNV')
# plt.savefig('/nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/graph/PCA_plot_SNV.png')
ax1.scatter(X[:, 0], X[:, 1], s=marker_size)
ax1.set_xlabel('PC1')
ax1.set_ylabel('PC2')
ax1.set_title('PCA plot of Bcell VAF for SNV')




# --------------------------------------------------- PCA for INS ---------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)
# use df_merged["POS_Variant_type"].str.contains('INS') to select INS
df_INS = df_merged[df_merged["POS_Variant_type"].str.contains('INS')]
pca.fit(df_INS[df_INS.columns[1:]])
X = pca.transform(df_INS[df_INS.columns[1:]])
# print(X)
# print(X.shape)

# --------------- plot ---------------
ax2 = fig.add_subplot(3, 2, 2)
# plt.figure(figsize=(8, 8))
# plt.scatter(X[:, 0], X[:, 1], s=marker_size)
# plt.xlabel('PC1')
# plt.ylabel('PC2')
# plt.title('PCA plot of Bcell VAF for INS')
# plt.savefig('/nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/graph/PCA_plot_INS.png')
ax2.scatter(X[:, 0], X[:, 1], s=marker_size)
ax2.set_xlabel('PC1')
ax2.set_ylabel('PC2')
ax2.set_title('PCA plot of Bcell VAF for INS')



# --------------------------------------------------- PCA for DEL ---------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)
# use df_merged["POS_Variant_type"].str.contains('DEL') to select DEL
df_DEL = df_merged[df_merged["POS_Variant_type"].str.contains('DEL')]
pca.fit(df_DEL[df_DEL.columns[1:]])
X = pca.transform(df_DEL[df_DEL.columns[1:]])
# print(X)
# print(X.shape)

# --------------- plot ---------------
# plt.figure(figsize=(8, 8))
# plt.scatter(X[:, 0], X[:, 1])
# plt.xlabel('PC1')
# plt.ylabel('PC2')
# plt.title('PCA plot of Bcell VAF for DEL')
# plt.savefig('/nfs/data05/otgk/rDNA/PSCA_variants/consensus/adhoc/graph/PCA_plot_DEL.png')

ax3 = fig.add_subplot(3, 2, 3)
ax3.scatter(X[:, 0], X[:, 1], s=marker_size)
ax3.set_xlabel('PC1')
ax3.set_ylabel('PC2')
ax3.set_title('PCA plot of Bcell VAF for DEL')


# --------------------------------------------------- PCA for SNV+INS ---------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)
# use df_merged["POS_Variant_type"].str.contains('SNV') to select SNV
df_SNV_INS = df_merged[df_merged["POS_Variant_type"].str.contains('SNV|INS')]
pca.fit(df_SNV_INS[df_SNV_INS.columns[1:]])
X = pca.transform(df_SNV_INS[df_SNV_INS.columns[1:]])
# print(X)
# print(X.shape)

# --------------- plot ---------------
ax4 = fig.add_subplot(3, 2, 4)
# plt.figure(figsize=(8, 8))
# plt.scatter(X[:, 0], X[:, 1], s=marker_size)
# plt.xlabel('PC1')
# plt.ylabel('PC2')
# plt.title('PCA plot of Bcell VAF for SNV+INS')

ax4.scatter(X[:, 0], X[:, 1], s=marker_size)
ax4.set_xlabel('PC1')
ax4.set_ylabel('PC2')
ax4.set_title('PCA plot of Bcell VAF for SNV+INS')


# --------------------------------------------------- PCA for SNV+DEL ---------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)
# use df_merged["POS_Variant_type"].str.contains('SNV') to select SNV
df_SNV_DEL = df_merged[df_merged["POS_Variant_type"].str.contains('SNV|DEL')]
pca.fit(df_SNV_DEL[df_SNV_DEL.columns[1:]])
X = pca.transform(df_SNV_DEL[df_SNV_DEL.columns[1:]])
# print(X)
# print(X.shape)


# --------------- plot ---------------
ax5 = fig.add_subplot(3, 2, 5)
ax5.scatter(X[:, 0], X[:, 1], s=marker_size)
ax5.set_xlabel('PC1')
ax5.set_ylabel('PC2')
ax5.set_title('PCA plot of Bcell VAF for SNV+DEL')


# --------------------------------------------------- PCA for INS+DEL ---------------------------------------------------
# --------------- PCA ---------------
pca = PCA(n_components=2)
# use df_merged["POS_Variant_type"].str.contains('SNV') to select SNV
df_INS_DEL = df_merged[df_merged["POS_Variant_type"].str.contains('INS|DEL')]
pca.fit(df_INS_DEL[df_INS_DEL.columns[1:]])
X = pca.transform(df_INS_DEL[df_INS_DEL.columns[1:]])
# print(X)
# print(X.shape)


# --------------- plot ---------------
ax6 = fig.add_subplot(3, 2, 6)
ax6.scatter(X[:, 0], X[:, 1], s=marker_size)
ax6.set_xlabel('PC1')
ax6.set_ylabel('PC2')
ax6.set_title('PCA plot of Bcell VAF for INS+DEL')





# ------------------------------------------------- save and show -------------------------------------------------
fig.tight_layout()
fig.savefig('/nfs/data05/otgk/rDNA/variantcall/PSCA_variants/consensus/adhoc/graph/PCA_plot_comb.png')




