#' The Combinatory Method (s) for the construction of rectangular PBIB designs
#'
#' The application of the Combinatory Method (s), with \eqn{s} chosen in \eqn{[2, l-1]},
#' on rectangular association scheme to obtain the configuration and the
#' parameters of the \code{PBIB} design associated.
#' @usage CombS(n, l, s)
#' @param n Number of lines of the association schemes array.
#' @param l Number of columns of the association schemes array.
#' @param s Number of the token treatments from the same row of the association scheme.
#' @return A LIST :
#'  \itemize{
#'   \item \code{PBIB } The configuration of the PBIB.
#'   \item \code{Type } The type of the design
#'   \item \code{V } Number of treatments.
#'   \item \code{B } Number of blocs.
#'   \item \code{R } Repetition of each treatment.
#'   \item \code{K } Size of blocs.
#'   \item \code{lamda } Vector of m-lambda.
#'   \item \code{Resolvable } Is the design Resolvable ?
#'   }
#' @details \itemize{
#'          \item For \eqn{2 < s < l}, we obtain a rectangular PBIB design.
#'          \item For \eqn{s = l}, we obtain a singular group divisible designs.
#'  }
#'
#' @author Mohamed Laib, Imane Rezgui, Zebida Gheribi-Aoulmi and Herve Monod
#' @references
#' Imane Rezgui, Z. Gheribi-Aoulmi (2014). New construction method
#' of rectangular partially balanced incomplete block designs and
#' singular group divisible designs,
#' \href{http://www.thescipub.com/jmss.toc}{Journal of Mathematics and Statistics, 10, 45- 48.}
#'
#' M.N. Vartak 1955. On an application of Kronecker product of Matrices to Statistical designs. Ann. Math. Stat.,26(420-438).
#'
#'
#' @seealso \code{\link{UType}}
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' n<-3
#' l<-3
#' s<-2
#' CombS(l,n,s)
#' }
#' @export
CombS <-function (n, l, s){
    V <- n * l
    if (s > l || s <= 1) {
        stop("Choose : 1 < s <= l")
    }

reso<-(n*l)%%(2*s)
bbo<-ifelse(reso==0,"Yes","No")
        lamda <- NULL

        BIB <- Op(n, l, s)

        T <- BIB[1, 1]
        R <- length(which(T == BIB))
        if (l ==s) {
            lamda[[1]] <- (n - 1) * choose(l - 2, s - 2)
            lamda[[2]] <- choose(l - 1, s - 1)
            lamda[[3]] <- choose(l - 2, s - 2)
            return(list(PBIB = BIB, Type = "Singular group divisible design",
                V = V, B = dim(BIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
        }
        else {
            lamda[[1]] <- R
            lamda[[2]] <- 1

            return(list(PBIB = BIB, Type = "Rectangular PBIB design",
                V = V, B = dim(BIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
        }

}

Op <- function(n, l, s) {
  V <- n * l
  mt <- matrix(1:V, nrow = n, byrow = TRUE)
  n <- dim(mt)[1]
  l <- dim(mt)[2]
  V <- n * l
  a <- combn(1:l, s)
  b <- combn(1:n, 2)
  v <- dim(b)[2]
  vv <- dim(a)[2]
  MAT <- NULL
  A <- 1
  y <- 1
  while (A <= vv) {
    for (k in 1:v) {
      zz <- a[, A]
      ss <- b[, k]
      MAT[[y]] <- as.vector(t(mt[ss, zz]))
      y <- y + 1
    }
    A <- A + 1
  }
  return(Reduce("rbind", MAT))
}
