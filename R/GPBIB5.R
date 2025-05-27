#' Generalized rectangular right angular (5) design.
#'
#' gives the configuration and the parametres of the design obtained
#' by the construction method of GPBIB_5 (see 3.2 of the paper rezgui
#' et al (2015)).
#' @usage GPBIB5(n, l, s, w)
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
#'   \item \code{lamda } Vector of m-lambda.
#'   \item \code{Resolvable } Is the design Resolvable ?
#'   }
#'
#' @author Mohamed Laib, Imane Rezgui, Zebida Gheribi-Aoulmi and Herve Monod
#' @references
#' Imane Rezgui, Z. Gheribi-Aoulmi and H. Monod (2015). U-type Designs via
#' New Generalized Partially Balanced Incomplete Block Designs with m = 4, 5
#' and 7 Associated Classes, \doi{10.4236/am.2015.62024}, Applied mathematics, 6, 242-264.
#'
#' Imane Rezgui, Z.Gheribi-Aoulmi and H. Monod, New association schemes
#' with 4, 5 and 7 associated classes and their associated partially
#' balanced incomplete block designs; Advances and Applications in Discrete
#' Mathematics Vol.12 Issue 2 197-206.
#'
#' @seealso \code{\link{UType}}
#' @note For \eqn{w=2}, the \code{GPBIB_5} is a rectangular right angular (5) (PBIB_5).
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' n<-3
#' l<-3
#' s<-3
#' w<-3
#' GPBIB5(n, l, s, w)
#' }
#' @export

GPBIB5 <-function(n,l,s,w){
if (s<3 & l<2 & n<2) stop("n,l should be greater than 1 and s greater than 2")
V<-n*l
reso<-(n*l)%%(2*s)
bbo<-ifelse(reso==0,"Yes","No")

A<-NULL;mat<-NULL;lamda<-NULL
for (i in 1:w){
A[[i]]<-matrix(1:V, ncol=l, byrow=TRUE)
z<-(i-1)*V
A[[i]]<-A[[i]]+z}

Bp<-NULL
for (j in 1:w) {
AA<-A[[j]]
AB<-A[-j]
M<-length(AB)


for (m in 1:M){
AS<-AB[[m]]

for (k in 1:l) {
co<-AA[,k]



mt<-cbind(co,AS)
X<-Opn(mt,s)
y<-dim(X)[1]
for (x in y:1){
if (any(X[x,1]==co)==FALSE){
X<-X[-x,]}}
Bp<-rbind(Bp,X)}}}
PBIB<-Bp
T <- PBIB[1, 1]
R <- length(which(T == PBIB))
lamda[1] <- l * (n - 1) * choose(l - 2, s - 3)*(w-1)
lamda[2] <- choose(l, s - 1) + (l * choose(l - 1, s - 2))*(w-1)
lamda[3] <- l * choose(l - 2, s - 3)*(w-1)
lamda[4] <- 2 * (n - 1) * choose(l - 1, s - 2)
lamda[5] <- 2 * choose(l - 1, s - 2)
return(list(PBIB = PBIB, Type = "Generalized rectangular right angular (5) (GPBIB_5) design", V = w * V, B = dim(PBIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
}
