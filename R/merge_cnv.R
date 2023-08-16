# Merging CNVs includes two steps: cluster CNVs and reassign CNVs to samples

#' Cluster same filtered CNVs from different Samples
#' @param cnvs_filter CNV list (usually filtered). Required columns:CHROM, ID, POS, ALT, END, SVLEN, QUAL_median, BSR_median, N
#' @param sample_cnv sample-level CNV list
#' @param cutoff.olcluster cutoff of overlaping of two CNVs
#' @param outputname output file name of merged CNV list
#' @export
merge_cnv <- function(cnvs_filter, sample_cnv, cutoff.olcluster, outputname){

  ####################################
  # Cluster CNVs
  ####################################

  cnv_overlap <- computeOl_2df(cnvs_filter, cnvs_filter) %>%
    filter(!(test.prop==1&base.prop==1)) %>%
    rowwise() %>%
    mutate(least.prop=min(test.prop, base.prop)) %>%
    mutate(highest.prop=max(test.prop, base.prop))

  # 2) baseID is the CNV with the smaller starting position
  cnv_overlap <- cnv_overlap %>%
    rowwise() %>%
    mutate(testID=max(test,base)) %>%
    mutate(baseID=min(test,base)) %>%
    select(-one_of("test","base", "test.prop","base.prop")) %>%
    unique()

  # 3) choose the pairs with at least overlap proportion >= cutoff.olcluster
  cnv_overlap <-cnv_overlap %>%
    filter(least.prop>=cutoff.olcluster)

  if(dim(cnv_overlap)[1]==0){

    cnv_final <- cnvs_filter %>%
      select(CHROM, ID, POS, ALT, END, SVLEN, QUAL_median, BSR_median, N) %>%
      rename(QUAL=QUAL_median, BSR=BSR_median) %>%
      mutate(groupid=0)

    #  write.table(cnv_final, paste0(outputname,"_indexcnvs_final.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

    # 8.3 ) add index_cnv column to cnv_uniqlist (the original unique CNV list)
    cnv_uniqlist_index <- cnvs_filter %>%
      select(CHROM, ID, POS, ALT, END, SVLEN) %>%
      mutate(index_CNV=ID)
    #  write.table(cnv_uniqlist_index, paste0(outputname,"_uniq_hardfiltcnvslist_index.txt"), sep = "\t", col.names = T, row.names = F, quote = F)


  }else{


    # 4) order overlap pairs by baseID
    cnv_overlap <- cnv_overlap[order(cnv_overlap$baseID),]

    # 5) group CNVs if they have >X overlap with at least one other CNV in the cluster
    cnv_overlap$groupid=0
    cnv_overlap$pairid=rownames(cnv_overlap)
    cnv_overlap[1,"groupid"]=1
    current_groupid=1
    unassigned_pair <- cnv_overlap[cnv_overlap$groupid==0,]
    current_list <- c(as.character(cnv_overlap[1,"testID"]), as.character(cnv_overlap[1,"baseID"]))
    assigned_list <- data.frame("assigned_cnv"=c(as.character(cnv_overlap[1,"testID"]), as.character(cnv_overlap[1,"baseID"])), "groupid"=1)

    while(dim(unassigned_pair)[1]>0){

      new_pairs <- unassigned_pair[((unassigned_pair$testID %in% current_list) | (unassigned_pair$baseID %in% current_list)),]
      if(dim(new_pairs)[1]>0){
        current_list <- unique(c(current_list, as.character(new_pairs$testID), as.character(new_pairs$baseID)))
        cnv_overlap[cnv_overlap$pairid %in% new_pairs$pairid, "groupid"] <- current_groupid
      }else{
        current_groupid <- current_groupid+1
        cnv_overlap[cnv_overlap$pairid==unassigned_pair[1,"pairid"],"groupid"] <- current_groupid
        current_list <- c(as.character(unassigned_pair[1,"testID"]), as.character(unassigned_pair[1,"baseID"]))
      }
      unassigned_pair <- cnv_overlap[cnv_overlap$groupid==0,]
      #  print(dim(unassigned_pair)[1])

    }

    # Grouping result of overlaping CNVs pairs
    #  write.table(cnv_overlap, paste0(outputname,"_olcnv_groups",".txt"), sep = "\t", col.names = T, row.names = F, quote = F)

    # 6) select index CNV for each group
    # Remember CNV's POS and END do not represent the true breakpoints, so it does not matter if we do not choose the real breakpoints for the index CNVs
    # Here the index CNV is rather a ID for the group of CNVs.
    #cnv_overlap <- read.table("cnv_overlap.txt", sep = "\t", header = T, as.is = T)
    # list of redundant CNVs
    cnvlist <- cnv_overlap[,c("testID","baseID","groupid")] %>%
      pivot_longer(1:2,names_to="type",values_to="ID") %>%
      select(ID, groupid) %>%
      unique()

    # add POS and END information for these CNVs.
    # for each group of CNVs, choose the minumum of POS and maximum of END to represent the index CNV.
    groupid_indexcnv <- cnvlist %>%
      mutate(TYPE=sub("\\-.*","",ID)) %>%
      left_join(cnvs_filter[,c("CHROM","ID","POS", "END", "QUAL_median","BSR_median", "N")], by="ID") %>%
      group_by(groupid,TYPE,CHROM) %>%
      summarise(POS=min(POS),
                END=max(END),
                QUAL=median(QUAL_median),
                BSR=median(BSR_median),
                N=n()) %>%
      mutate(index_CNV=paste0(TYPE, ":",CHROM,"-",POS, "-",END))


    # 7) Using index CNV for CNVs overlapping with other CNVs
    #
    olcnvlist_index <- cnvlist %>%
      left_join(groupid_indexcnv, by="groupid")%>%
      mutate(ALT=case_when(
        TYPE=="G" ~ "DUP",
        TYPE=="L" ~ "DEL"))
    #  write.table(olcnvlist_index, paste0(outputname,"_olcnvlist_index.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

    # 8) Using groupid_indexcnv to replace the CNVs in the group
    # 8.1) rename group CNVs and generate ALT column
    groupid_indexcnv_rename <- groupid_indexcnv %>%
      dplyr::rename(ID=index_CNV) %>%
      mutate(ALT=case_when(
        TYPE=="G" ~ "DUP",
        TYPE=="L" ~ "DEL"))


    # 8.2) Remove CNVs that are in olcnvlist_index and add group CNVs (multiple overlapping CNVs will be replaced by one index CNV)
    cnv_final <- cnvs_filter %>%
      ungroup() %>%
      filter(!ID %in% olcnvlist_index$ID) %>%
      select(CHROM, ID, POS, ALT, END, QUAL_median, BSR_median, N) %>%
      dplyr::rename(QUAL=QUAL_median, BSR=BSR_median) %>%
      mutate(groupid=0) %>%
      add_row(groupid_indexcnv_rename[,c("CHROM","ID","POS","ALT","END","QUAL","BSR","N","groupid")]) %>%
      mutate(SVLEN=END-POS)

    print(paste0("Number of CNVs after merging:", dim(cnv_final)[1]))

    write.table(cnv_final, paste0(outputname,"_merged_uniq_cnvlist.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

    # 8.3 ) add index_cnv column to cnv_uniqlist (the original unique CNV list)
    cnv_uniqlist_index <- cnvs_filter %>%
      ungroup() %>%
      filter(!ID %in% olcnvlist_index$ID) %>%
      select(CHROM, ID, POS, ALT, END) %>%
      mutate(index_CNV=ID) %>%
      add_row(olcnvlist_index[,c("CHROM", "ID", "POS", "ALT", "END", "index_CNV")])
    write.table(cnv_uniqlist_index, paste0(outputname,"_indexCNVmap.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

  }

  ####################################
  # reassign CNVs to samples
  ####################################

  # Sample-level CNV list

  # keep the cnvs that are in hardfiltered_CNV_index
  cnv_fl <- cnv_uniqlist_index[,c("ID","index_CNV")] %>%
    inner_join(sample_cnv,by="ID")

  # choosing cnvs from the finalCNV calls for unique representation of index_CNVs across different samples
  sample_cnv_final <- cnv_final %>%
    left_join(cnv_fl[,c("index_CNV","CN","SampleID")], by=c("ID"="index_CNV"))

  write.table(sample_cnv_final, paste0(outputname,"_merged_sample_cnvlist.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

}
