# script to obtain independent de novo events (ie exclude de novos that recur 
# within families).

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
TRIOS_PATH = file.path(DATAFREEZE_DIR, "family_relationships.txt")

#' opens a file defining the DDD individuals, and their families
open_ddd_families <- function(path) {
    families = read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    families$affected = NULL
    families$path_to_vcf = NULL
    families$sex = NULL
    
    return(families)
}

#' remove de novos that are shared between multiple probands of a family
get_independent_de_novos <- function(de_novos) {
    families = open_ddd_families(TRIOS_PATH)
    
    # merge the family IDs with the de novo events, sort the de novos first, so
    # that we can correctly match to the duplicate status
    if ("person_stable_id" %in% names(de_novos)) {
        de_novos = de_novos[order(de_novos$person_stable_id), ]
        dups = merge(de_novos, families, by.x="person_stable_id", by.y="individual_id", all.x=TRUE)
    } else {
        de_novos = de_novos[order(de_novos$person_id), ]
        dups = merge(de_novos, families, by.x="person_id", by.y="individual_id", all.x=TRUE)
    }
    
    # get the family ID, and de novo coordinates, which are sufficient to 
    # identify duplicates within families
    if ("pos" %in% names(dups)) {
        dups = dups[, c("family_id", "chrom", "pos", "ref", "alt")]
    } else {
        dups = dups[, c("family_id", "chrom", "start_pos", "ref_allele", "alt_allele")]
    }
    
    # restrict ourselves to the non-duplicates (this retains the first de novo 
    # for each family)
    de_novos = de_novos[!duplicated(dups), ]
    
    return(de_novos)
}

