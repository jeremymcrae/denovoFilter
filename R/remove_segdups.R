
#' excludes de novo calls within segdup regions
#'
#' @param de_novos dataframe of candidate de novo calls
#' @param segdup_path path to file listing segmantal duplication regions
#'        throughout the genome.
#' @export
#'
#' @return dataframe of candidate de novo calls with sites in segmental
#'         duplications removed.
exclude_segdups <- function(de_novos, segdup_path) {
    
    # load all the segdup regions
    segdups = read.table(segdup_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    segdups$chrom = gsub("chr", "", segdups$chrom)
    
    in_segdups = apply(de_novos, 1, function(x) check_in_segdup(x, segdups))
    
    return(de_novos[!in_segdups, ])
}

#' check if a de novo call is within the segdup regions
#'
#' @param variant dataframe/list for single variant
#' @param segdups dataframe of segmental duplication regions in the genome
#' @export
#'
#' @return true/false for whwther the site is in  a segmental duplication.
check_in_segdup <- function(variant, segdups) {
    
    chrom = variant[["chrom"]]
    pos = as.numeric(variant[["pos"]])
    
    regions = segdups[segdups$chrom == chrom, ]
    in_segdup = any((pos >= regions$chromStart) & (pos <= regions$chromEnd))
    
    return(in_segdup)
}
