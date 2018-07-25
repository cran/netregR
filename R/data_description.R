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




#' Social interaction data set
#'
#' A synthetic data set of standardized, directed interactions between 25 students in a seventh grade class. 
#' 
#' 
#' @name interactions
#' @docType data
#'
#' @format A data set with four variables. Includes the true parameters.
#' \describe{
#'	\item{interactions}{25 x 25 x 5 numeric array of directed relations}
#' 	\item{xbinary}{25 x 25 x 5 numeric array of  binary indictors}
#' 	\item{xabs}{25 x 25 x 5 numeric array of standardized absolute difference in indicated interest in each subject area}
#' 	\item{betatrue}{Numeric vector of length 7 that contains true coefficients. The first two (value 1) pertain to \code{shared_project} and \code{grade_difference_abs}. The last three are separate intercepts for each observation of the network.}
#' 	\item{Omegatrue}{3000 x 3000 numeric matrix: the true covariance matrix of the errors.}
#' 	\item{phitrue}{2x6 numeric matrix: true parameters of covariance matrix.}
#' 	}
#'
#' @details We generated a symthetic data set form a true linear model with jointly exchangeable errors. The interactions (the outcomes) between 25 students represent normalized, directed relations between them in 5 different contexts (texts). The observation could be, for instance, the standardized number of characters texted from one student to another over a month pertaining to five subjects: school, friends, family, significant others, and popular culture. The first covariate, \code{xbinary}, indicates whether both students indicated in a survey that they were interested in each topic. The second covariate, \code{xabs}, measures the absolute, standardized difference in number of characters in total texts of each student of each subject area. 
#'
#'
#' @keywords datasets
#'
#' @examples
#' data("interactions")
#' 
NULL

