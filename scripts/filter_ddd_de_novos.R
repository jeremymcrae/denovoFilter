# filter DNG callset initially on basis of

library(denovoFilter)
library(reshape)
library(ggplot2)

# DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
# DE_NOVOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/denovo_gear_trios_extracted_passed_variants_19.02.15.tsv"
DE_NOVOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/denovo_gear_trios_extracted_passed_variants_11.05.15.tsv"
TRIOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"
VALIDATIONS_PATH = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv"

open_validations <- function(path) {
    validations = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    validations = validations[, c("person_stable_id", "chr", "pos", "ref",
        "alt", "validation_result")]
    names(validations) = c("person_stable_id", "chrom", "pos", "ref", "alt", "result")
    
    # convert the result codes to human-interpretable strings
    codes = data.frame(result=c("DNM", "FP", "INH", "P/U", "R", "U"),
        status=c("de_novo", "false_positive", "inherited", "partially_uncertain", "repeat",
        "uncertain"))
    validations = merge(validations, codes, by="result")
    validations$result = NULL
    validations$status = as.character(validations$status)
    
    return(validations)
}

main <- function() {
    validations = open_validations(VALIDATIONS_PATH)
    
    de_novos = read.table(DE_NOVOS_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    de_novos = preliminary_filtering(de_novos)
    de_novos = fix_missing_gene_symbols(de_novos)
    de_novos = extract_alt_and_ref_counts(de_novos)
    de_novos = count_gene_recurrence(de_novos)
    
    de_novos = test_sites(de_novos)
    de_novos = test_genes(de_novos)
    
    pass_status = get_filter_status(de_novos)
    passed = de_novos[pass_status & de_novos$coding, ]
    passed = get_independent_de_novos(passed, TRIOS_PATH)
    
    # check_validations(de_novos, validations)
    # check_validations(passed, validations)
    
    de_novos = subset_de_novos(de_novos)
    write.table(passed, file="~/de_novos.ddd_4k.ddd_only.txt",
        quote=FALSE, row.names=FALSE, sep="\t")
}

main()
