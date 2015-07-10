
#' check filtering performance against validated de novos
#'
#' @param de_novos dataframe of candidate de novo sites
#' @param validations dataframe of de novo sites with validation status
#' @export
#'
#' @return something
check_validations <- function(de_novos, validations, plot_pp_dnm=FALSE) {
    
    # annotate the de novos with their validation status
    de_novos = merge(de_novos, validations, by=c("person_stable_id", "chrom", "pos", "ref", "alt"), all.x=TRUE)
    de_novos$status[is.na(de_novos$status)] = "new_unknown"
    
    # determine the number of variants in each validation category
    frequencies = data.frame(table(de_novos$status))
    names(frequencies) = c("status", "frequency")
    
    if (plot_pp_dnm) {
        by_status = split(de_novos, de_novos$status)
        for (status in by_status) {
            plot(density(status$pp_dnm), main=unique(status$status))
        }
        dev.off()
    }
    
    return(frequencies)
}
