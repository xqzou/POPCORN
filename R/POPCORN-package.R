#' @keywords internal
"_PACKAGE"

## usethis namespace: start
utils::globalVariables(c("telomere", "centromere","gap","sg","srep"))
#' @importFrom dplyr mutate filter group_by ungroup summarise select left_join n
#' @importFrom tidyr replace_na
#' @importFrom stats median
#' @importFrom GenomicRanges findOverlaps pintersect width
#' @importFrom utils write.table
## usethis namespace: end
NULL
