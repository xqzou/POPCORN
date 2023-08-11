#' Find overlaping regions between two any data.frames of 2 lists
#' @param df.test test data.frame. Required columns:"CHROM","POS","ID","END"
#' @param df.base base data.frame. Required columns:"CHROM","POS","ID","END"
#' @export
computeOl_2any <- function(df.test,df.base){
  in.dat.test <- df.test
  test.gr <- with(in.dat.test,GRanges(seqnames=Rle(CHROM),ranges=IRanges(start=as.numeric(POS),end=as.numeric(END)),uid=ID))

  in.dat.base <- df.base
  base.gr <- with(in.dat.base,GRanges(seqnames=Rle(CHROM),ranges=IRanges(start=as.numeric(POS),end=as.numeric(END)),uid=ID))

  # message(f.test)
  hits <- suppressWarnings(findOverlaps(test.gr, base.gr))
  overlaps <- suppressWarnings(pintersect(test.gr[queryHits(hits)], base.gr[subjectHits(hits)]))
  percentOverlap.test <- width(overlaps) / width(test.gr[queryHits(hits)])
  percentOverlap.base <- width(overlaps) / width(base.gr[subjectHits(hits)])
  ## one.way results
  ol.dt <- data.frame(test=test.gr[queryHits(hits),]$uid,
                      base=base.gr[subjectHits(hits),]$uid,
                      test.prop=percentOverlap.test,
                      base.prop=percentOverlap.base,
                      overlap.width=width(overlaps))

  ol.dt <- ol.dt %>%
    filter(overlap.width>1) %>%
    unique()

  return(ol.dt)
}


