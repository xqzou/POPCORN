#' Annotate overlapping of a list of unique autosomal CNVs with the following regions
#' 1) telomere / centromere
#' 2) Gap (N-masked)
#' 3) Segmental duplications
#' 4) Simple Repeats
#' 5) MHC
#' @param cnv uniq CNV calls. Required columns: SampleID, CHROM, ID, POS, ALT, SVLEN, END, QUAL, BSR, CN
#' @export
anno_cnv <- function(cnv){

  cnv$p_tlm_dis <- 0
  cnv$q_tlm_dis <- 0
  cnv$tlm_dis_nearest <- 0
  cnv$ctm_dis_p <- 0
  cnv$ctm_dis_q <- 0
  cnv$ctm_dis <- 0
  chr_seq <- cnv %>% select(CHROM) %>% unique()

  chr_seq_all <- c(paste0("chr",seq(1,22)),"chrX","chrY")
  if(!all(chr_seq %in% chr_seq_all)){
    stop("The CHROM column should contain entries from the following set:\n",
        "(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY)\n",
         "Exiting...",call.=FALSE)
  }

  ##############################################################
  # Annotate relative distance to telomere and centromere
  # Normalize each arm to 1 according to the arm length
  ##############################################################
  # loop by chr
  for (j in 1:length(chr_seq)) {

    chr_j <- chr_seq[j]
    #   cnv_chr <- cnv[cnv$CHROM==chr_j,]
    #  cnv_chr$center <- round(cnv[cnv$CHROM==chr_seq[j],]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2)

    p_tlm_start <- telomere[telomere$chrom==chr_j & telomere$ix==1, "chromStart"]
    p_tlm_end <- telomere[telomere$chrom==chr_j & telomere$ix==1, "chromEnd"]
    q_tlm_start <- telomere[telomere$chrom==chr_j & telomere$ix!=1, "chromStart"]
    q_tlm_end <- telomere[telomere$chrom==chr_j & telomere$ix!=1, "chromEnd"]

    ctm_center <- round(centromere[centromere$chrom==chr_j, "start"] + (centromere[centromere$chrom==chr_j, "end"] - centromere[centromere$chrom==chr_j, "start"])/2)
    #chr_len <- telomere[telomere$chrom==chr_j & telomere$ix!=1, "chromEnd"] - telomere[telomere$chrom==chr_j & telomere$ix==1, "chromStart"]
    p_len <- ctm_center - telomere[telomere$chrom==chr_j & telomere$ix==1, "chromStart"]
    q_len <- telomere[telomere$chrom==chr_j & telomere$ix!=1, "chromEnd"] - ctm_center

    #cnv[cnv$CHROM==chr_j,]$p_tlm_dis <- 100*(round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2) - p_tlm_end)/chr_len
    #cnv[cnv$CHROM==chr_j,]$q_tlm_dis <- 100*(q_tlm_start - round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2))/chr_len
    #cnv[cnv$CHROM==chr_j,]$ctm_dis <- 100*(ctm_center - round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2))/chr_len
    #cnv[cnv$CHROM==chr_j,]$tlm_dis_nearest <- ifelse(cnv[cnv$CHROM==chr_j,]$p_tlm_dis<cnv[cnv$CHROM==chr_j,]$q_tlm_dis, cnv[cnv$CHROM==chr_j,]$p_tlm_dis, cnv[cnv$CHROM==chr_j,]$q_tlm_dis)

    cnv[cnv$CHROM==chr_j,]$p_tlm_dis <- (round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2) - p_tlm_start)/p_len
    cnv[cnv$CHROM==chr_j,]$q_tlm_dis <- (q_tlm_end - round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2))/q_len
    cnv[cnv$CHROM==chr_j,]$ctm_dis_p <- (ctm_center - round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2))/p_len
    cnv[cnv$CHROM==chr_j,]$ctm_dis_q <- (round(cnv[cnv$CHROM==chr_j,]$POS + abs(cnv[cnv$CHROM==chr_j,]$SVLEN)/2)-ctm_center)/q_len

    cnv[cnv$CHROM==chr_j,]$tlm_dis_nearest <- ifelse(cnv[cnv$CHROM==chr_j,]$p_tlm_dis<cnv[cnv$CHROM==chr_j,]$q_tlm_dis, cnv[cnv$CHROM==chr_j,]$p_tlm_dis, cnv[cnv$CHROM==chr_j,]$q_tlm_dis)
    cnv[cnv$CHROM==chr_j,]$ctm_dis <- ifelse(cnv[cnv$CHROM==chr_j,]$ctm_dis_p > cnv[cnv$CHROM==chr_j,]$ctm_dis_q, cnv[cnv$CHROM==chr_j,]$ctm_dis_p, cnv[cnv$CHROM==chr_j,]$ctm_dis_q)

  }

  ##############################################################
  # Annotate overlaps with gap regions
  ##############################################################
  cnv_gap <- computeOl_2any(cnv, gap) %>%
    dplyr::rename(ID=test, Gap=base, CNV_prop_Gap=test.prop, Gap_prop=base.prop, overlap_width_Gap=overlap.width)

  cnv_gap <- cnv_gap %>%
    group_by(ID) %>%
    summarise(CNV_prop_Gap=sum(CNV_prop_Gap),
              overlap_width_Gap=sum(overlap_width_Gap))

  cnv <- cnv %>%
    left_join(cnv_gap,by="ID") %>%
    replace_na(list(CNV_prop_Gap=0, overlap_width_Gap=0))

  ##############################################################
  # Annotate MHC CNVs
  # MHC (HLA): chr6:28,510,120-33,480,577
  ##############################################################
  mhc <- data.frame("CHROM" = c("chr6"), "POS"=c(28510120), "END"=c(33480577), "ID"=c("chr6:28510120-33480577"))
  cnv_mhc <- computeOl_2any(cnv, mhc) %>%
    dplyr::rename(ID=test, MHC=base, CNV_prop_MHC=test.prop, MHC_prop=base.prop, overlap_width_MHC=overlap.width)

  cnv_mhc <- cnv_mhc %>%
    group_by(ID) %>%
    summarise(CNV_prop_MHC=sum(CNV_prop_MHC),
              overlap_width_MHC=sum(overlap_width_MHC))

  cnv <- cnv %>%
    left_join(cnv_mhc,by="ID") %>%
    replace_na(list(CNV_prop_MHC=0, overlap_width_MHC=0))

  # cnv <- cnv %>%
  #   mutate(MHC=case_when(
  #     CHROM=="chr6" & (END > 28510120 & END < 33480577) ~ "mhc",
  #     CHROM=="chr6" & (POS > 28510120 & POS < 33480577) ~ "mhc",
  #     TRUE ~ "non_mhc"))

  ##############################################################
  # Annotate overlaps with segmental duplications regions
  ##############################################################

  cnv_sg <- computeOl_2any(cnv, sg) %>%
    dplyr::rename(ID=test, sg=base, CNV_prop_sg=test.prop, sg_prop=base.prop, overlap_width_sg=overlap.width)

  cnv_sg <- cnv_sg %>%
    group_by(ID) %>%
    summarise(CNV_prop_sg=max(CNV_prop_sg),
              overlap_width_sg=max(overlap_width_sg))

  cnv <- cnv %>%
    left_join(cnv_sg,by="ID") %>%
    replace_na(list(CNV_prop_sg=0, overlap_width_sg=0))

  ##############################################################
  # Annotate overlaps with simple repeats
  ##############################################################
  cnv_srep <- computeOl_2any(cnv, srep) %>%
    dplyr::rename(ID=test, srep=base, CNV_prop_srep=test.prop, srep_prop=base.prop, overlap_width_srep=overlap.width)

  cnv_srep <- cnv_srep %>%
    group_by(ID) %>%
    summarise(CNV_prop_srep=max(CNV_prop_srep),
              overlap_width_srep=max(overlap_width_srep))

  cnv <- cnv %>%
    left_join(cnv_srep,by="ID") %>%
    replace_na(list(CNV_prop_srep=0, overlap_width_srep=0))

  return(cnv)

}
