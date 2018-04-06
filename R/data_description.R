#' Wolf network data set
#'
#' A data set measuring dominance and its behavioral measures in a captive wolf pack.
#' 
#' 
#' @name wolf
#' @docType data
#'
#' @format A data set with three variables
#' \describe{
#'	\item{wolf}{16 x 16 numeric matrix of dominance measures}
#' 	\item{wolf_age_diff}{16 x 16 numeric matrix of difference in ages (column less row)}
#' 	\item{wolf_same_sex}{16 x 16 numeric matrix of indicators of same sex}
#'}
#'
#' @details  This is data on a captive family of wolves in Arnheim, Germany. The 16 wolves studied here were housed in a large wooded enclosure and observed in 1978. This matrix displays deference acts. The number in a cell represents the number of occasions on which the row wolf was seen to exhibit a "low posture" display directed toward the column wolf. The behavior could involve approach or retreat, but the fact that it was performed in "low posture" suggests that it was deferent. Data obtained March 20, 2018 from \url{http://moreno.ss.uci.edu/data.html#wolf}. 
#'
#' @source \url{http://moreno.ss.uci.edu/data.html#wolf}
#'
#' @references Jan A. R. A. M. van Hooff and Joep A. B. Wensing, "Dominance and its behavioral measures in a captive wolf pack," Chapter 11 in Harry Frank, ed., Man and Wolf. Dordrecht: Junk, 1987, pp. 219-252. 
#'
#' @keywords datasets
#'
#' @examples
#' data("wolf")
#' 
NULL
