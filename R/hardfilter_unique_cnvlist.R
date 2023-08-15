#' Apply hardfilters to annotated CNV list
#' @param cnv Annotated CNV data.frame Required columns: SampleID, CHROM, ID, POS, ALT, SVLEN, END, QUAL_median, BSR_median
#' @param cutoff.svlen cutoff of CNV size
#' @param cutoff.qual cutoff of CNV QUAL_median cross population
#' @param cutoff.bsr cutoff of CNV BSR_median cross population
#' @param cutoff.tlmdis cutoff of distance to telomere
#' @param cutoff.ctmdis cutoff of distance to centromere
#' @param cutoff.olgap cutoff of overlaping of CNV with gap region
#' @param cutoff.olsg cutoff of overlaping of CNV with segmental duplication
#' @param cutoff.olsrep cutoff of overlaping of CNV with simple repeats
#' @param cutoff.olmhc cutoff of overlaping of CNV with MHC
#' @export
hardfilter_unique_cnvlist <- function(cnv,cutoff.svlen, cutoff.qual, cutoff.bsr, cutoff.tlmdis, cutoff.ctmdis, cutoff.olgap, cutoff.olsg, cutoff.olsrep, cutoff.olmhc){

  print(paste0("cutoff.svlen: ", cutoff.svlen))
  print(paste0("cutoff.qual: ", cutoff.qual))
  print(paste0("cutoff.bsr: ", cutoff.bsr))
  print(paste0("cutoff.tlmdis: ", cutoff.tlmdis))
  print(paste0("cutoff.ctmdis: ", cutoff.ctmdis))
  print(paste0("cutoff.olgap: ", cutoff.olgap))
  print(paste0("cutoff.olsg: ", cutoff.olsg))
  print(paste0("cutoff.olsrep: ", cutoff.olsrep))
  print(paste0("cutoff.olmhc: ", cutoff.olmhc))

  cnvs_filter <- cnv %>%
    filter(abs(SVLEN)>=cutoff.svlen) %>%
    filter(QUAL_median>=cutoff.qual) %>%
    filter(BSR_median>=cutoff.bsr) %>%
    filter(tlm_dis_nearest>=cutoff.tlmdis) %>%
    filter(ctm_dis>=cutoff.ctmdis) %>%
    filter(CNV_prop_Gap<cutoff.olgap) %>%
    filter(CNV_prop_sg<cutoff.olsg) %>%
    filter(CNV_prop_srep<cutoff.olsrep)%>%
    filter(CNV_prop_MHC<cutoff.olmhc)

  return(cnvs_filter)
}
