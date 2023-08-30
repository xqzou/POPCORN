#' Apply hardfilters of QUAL for DEL and DUP separately to annotated unique CNV segment list
#' @param cnv Annotated CNV data.frame Required columns: CHROM, ID, POS, ALT, SVLEN, END, QUAL_median, BSR_median
#' @param cutoff.qual.dup cutoff of DUP QUAL_median cross population
#' @param cutoff.qual.del cutoff of DEL QUAL_median cross population
#' @export
hardfilter_QUAL_type <- function(cnv,cutoff.qual.del, cutoff.qual.dup){

  print(paste0("cutoff.svlen: ", cutoff.svlen))
  print(paste0("cutoff.qual.del: ", cutoff.qual.del))
  print(paste0("cutoff.qual.dup: ", cutoff.qual.dup))


  cnvs_filter_del <- cnv %>%
    filter(ALT=="DEL" & QUAL_median>=cutoff.qual.del)

  cnvs_filter_dup <- cnv %>%
    filter(ALT=="DUP" & QUAL_median>=cutoff.qual.dup)



  return(rbind(cnvs_filter_del, cnvs_filter_dup))
}
