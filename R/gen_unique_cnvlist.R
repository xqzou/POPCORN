#' Generate a list of unique CNVs across all given samples
#' Re-calibrate the CNV QUAL, BSR
#' @param sample_cnvs A data.frame of Sample-level CNV calls. Required columns: SampleID, CHROM, ID, POS, ALT, SVLEN, END, QUAL, BC, CN
#' @export
gen_unique_cnvlist <- function(sample_cnvs){
  sample_cnvs <- sample_cnvs %>%
    mutate(BSR=1000*BC/abs(SVLEN))
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
    mutate('CHROM'=paste0("chr",CHROM)) %>%
    ungroup()

  return(cnv_uniqlist)
}
