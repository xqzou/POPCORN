#' Generate a list of unique CNVs across all given samples
#' Re-calibrate the CNV QUAL, BSR 
#' @param sample_cnvs Sample-level CNV calls. Required columns: SampleID, CHROM, ID, POS, ALT, SVLEN, END, QUAL, BSR, CN
#' @export
gen_unique_cnvlist <- function(sample_cnvs){
  sample_cnvs$BSR <- 1000*sample_cnvs[,"BC"]/abs(sample_cnvs$SVLEN)
  cnv_uniqlist <- sample_cnvs %>%
    group_by(CHROM, ID, POS, ALT, SVLEN, END) %>%
    summarise(QUAL_median=median(QUAL),
              QUAL_max=max(QUAL),
              BSR_median=median(BSR),
              CN_median=median(CN),
              CN_mean=mean(CN),
              CN_min=min(CN),
              CN_max=max(CN),
              N=n()) %>%
    mutate('CHROM'=paste0("chr",CHROM))
  
  return(cnv_uniqlist)
}