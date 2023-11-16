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

# To run this script:
# conda activate rDNAVariantCall
# python3 merge_vaf.py

