#' @keywords internal
"_PACKAGE"

## usethis namespace: start
utils::globalVariables(c("telomere", "centromere","gap","sg","srep"))
#' @importFrom dplyr mutate filter group_by ungroup summarise select left_join inner_join n one_of rowwise rename case_when add_row
#' @importFrom tidyr replace_na pivot_longer
#' @importFrom stats median
#' @importFrom GenomicRanges GRanges findOverlaps pintersect width
#' @importFrom utils write.table
## usethis namespace: end
NULL
