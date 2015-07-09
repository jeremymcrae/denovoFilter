# estimated error rate at 0.0012 from DNM calls in parents in DDG2P genes, set
# slightly higher to be conservative
ERROR_RATE = 0.002

# threshold for removing sites with high strand bias, or parental alt frequency
P_CUTOFF = 1e-3

devtools::use_data(ERROR_RATE, P_CUTOFF, internal=TRUE, overwrite=TRUE)
