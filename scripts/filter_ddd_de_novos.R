# filter DNG callset initially on basis of

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
DE_NOVOS_PATH = file.path(DATAFREEZE_DIR, "denovo_gear_trios_extracted_passed_variants_11.05.15.tsv")
TRIOS_PATH = file.path(DATAFREEZE_DIR, "family_relationships.txt")

main <- function() {
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
    
    de_novos = subset_de_novos(de_novos)
    write.table(passed, file="~/de_novos.ddd_4k.ddd_only.txt",
        quote=FALSE, row.names=FALSE, sep="\t")
}

main()
