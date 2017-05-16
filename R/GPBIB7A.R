#'  Generalized rectangular right angular (7) design with \eqn{\lambda_{i}}{lambda(i)}
#'  equal to \eqn{\lambda_{i+4}}{lamda(i+4)} \eqn{(i=1,...,4)}
#'
#' gives the configuration and the parametres of the design obtained by
#' the first construction method of \code{GPBIB_7} (see 3.3.1 of the paper
#' rezgui et al (2015))
#' @usage GPBIB7A(n, l, s, w)
#' @param n Number of lines of the association schemes array.
#' @param l Number of columns of the association schemes array.
#' @param s Number of the token treatments from the same row of the association scheme.
#' @param w Number of the association scheme arrays.
#' @return A LIST :
#'  \itemize{
#'   \item \code{PBIB } The configuration of the PBIB.
#'   \item \code{Type } The type of the design
#'   \item \code{V } Number of treatments.
#'   \item \code{B } Number of blocs.
#'   \item \code{R } Repetition of each treatment.
#'   \item \code{K } Size of blocs.
#'   \item \code{lambda } Vector of m-lambda.
#'   \item \code{Resolvable } Is the design Resolvable ?
#'   }
#'
#' @author Mohamed Laib, Imane Rezgui, Zebida Gheribi-Aoulmi and Herve Monod
#' @references
#' Imane Rezgui, Z. Gheribi-Aoulmi and H. Monod (2015). U-type Designs via
#' New Generalized Partially Balanced Incomplete Block Designs with m = 4, 5
#' and 7 Associated Classes,
#' \href{http://dx.doi.org/10.4236/am.2015.62024}{Applied mathematics, 6, 242-264.}
#'
#' Imane Rezgui, Z.Gheribi-Aoulmi and H. Monod, New association schemes
#' with 4, 5 and 7 associated classes and their associated partially
#' balanced incomplete block designs; Advances and Applications in Discrete
#' Mathematics Vol.12 Issue 2 197-206.
#'
#' @seealso \code{\link{GPBIB7B}} and \code{\link{UType}}
#' @note For \eqn{w=2}, the \code{GPBIB_7} is a rectangular right angular (7) (PBIB_7).
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' n<-3
#' l<-3
#' s<-3
#' w<-3
#' GPBIB7A(n, l, s, w)
#' }
#' @export

GPBIB7A <-function(n,l,s,w){
  if (s<2 & l<2 & n<2) stop("n,l,s should be greater than 1")
  V<-n*l
  reso<-(n*l)%%(2*s)
  bbo<-ifelse(reso==0,"Yes","No")

  A<-NULL;mat<-NULL;lamda<-NULL
  for (i in 1:w){
    A[[i]]<-matrix(1:V, ncol=l, byrow=TRUE)
    z<-(i-1)*V
    A[[i]]<-A[[i]]+z
  }

  Bp<-NULL
  for (j in 1:w) {
    X<-Opn(A[[j]],s)
    Bp<-cbind(Bp,X)
  }
  BIB<-Bp
  T <- BIB[1, 1]
  R <- length(which(T == BIB))
  lamda[1] <- (n - 1) * choose(l - 2, s - 2)
  lamda[2] <- choose(l - 1, s - 1)
  lamda[3] <- choose(l - 2, s - 2)
  lamda[4] <- R
  lamda[5] <- lamda[1]
  lamda[6] <- lamda[2]
  lamda[7] <- lamda[3]
  return(list(PBIB = BIB, Type = "Generalized rectangular right angular (7) (GPBIB_7) designs", V = w * V, B = dim(BIB)[1], R = R, K = w*2 * s, lamda = lamda, Resolvable=bbo))
}
