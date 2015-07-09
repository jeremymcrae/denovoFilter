
#' adds gene symbols to variants lacking them, using a mupit function
#'
#' @param de_novos dataframe of de novo variants
#' @export
#'
#' @return data frame of de novos, but with additional annotations for many
#'         variants previously lacking a HGNC symbol.
fix_missing_gene_symbols <- function(de_novos) {
    
    # get the variants with no gene annotation, ensure chrom, start and stop
    # positions columns exist
    missing_genes = de_novos[de_novos$symbol == "", ]
    missing_genes$start_pos = as.character(missing_genes$pos)
    missing_genes$end_pos = as.character(as.numeric(missing_genes$start_pos) +
        (nchar(missing_genes$ref) - 1))
    
    # find the HGNC symbols (if any) for the variants
    hgnc_symbols = apply(missing_genes, 1, mupit::get_gene_id_for_variant, verbose=TRUE)
    
    # add the HGNC symbols to the rows that need it
    de_novos$symbol[de_novos$symbol == ""] = hgnc_symbols
    
    # 360 out of 17000 de novos still lack HGNC symbols. Their consequences are:
    #
    #   consequence                count
    #   ========================   =====
    #   downstream_gene_variant      17
    #   intergenic_variant          259
    #   regulatory_region_variant    63
    #   upstream_gene_variant        27
    #
    # In spot checks, these are sufficiently distant from genes that we can't
    # add them to the analysis of their nearest gene. We shall analyse these
    # per site by giving them mock gene symbols.
    missing_genes = de_novos[de_novos$symbol == "", ]
    de_novos$symbol[de_novos$symbol == ""] = paste("fake_symbol.",
        missing_genes$chrom, "_", missing_genes$pos, sep="")
    
    return(de_novos)
}
