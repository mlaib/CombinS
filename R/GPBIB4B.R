#' Generalized rectangular right angular (4) design with \eqn{\lambda_4}{lambda4} not equal to 0
#'
#' Gives the configuration and the parametres of the design obtained by
#' the seconde construction method of GPBIB_4 (see 3.1.2 of the paper
#' rezgui et al (2015)).
#' @usage GPBIB4B(n, l, s, w)
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
#' and 7 Associated Classes,
#' \href{http://dx.doi.org/10.4236/am.2015.62024}{Applied mathematics, 6, 242-264.}
#'
#' Imane Rezgui, Z.Gheribi-Aoulmi and H. Monod, New association schemes
#' with 4, 5 and 7 associated classes and their associated partially
#' balanced incomplete block designs; Advances and Applications in Discrete
#' Mathematics Vol.12 Issue 2 197-206.
#'
#' @seealso \code{\link{GPBIB4A}} and \code{\link{UType}}
#' @note For \eqn{w=2}, the \code{GPBIB_4} is a rectangular right angular (4) (PBIB_4)
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' n<-3
#' l<-3
#' s<-3
#' w<-3
#' GPBIB4B(n, l, s, w)
#' }
#' @export

GPBIB4B <-function(n,l,s,w){
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
vec<-NULL

for (m in 1:M){
AS<-AB[[m]]

for (k in 1:l) {
co<-AA[,k]
vec[[1]]<-AA[,k]

for (p in 2:n) {
vec[[p]]<-c(co[-1],co[1])
co<-c(co[-1],co[1])}
N<-length(vec)

for (p in 1:N) {

mt<-cbind(vec[[p]],AS)
X<-Opn(mt,s)
y<-dim(X)[1]
nn<-vec[[p]]
for (x in y:1){
if (any(X[x,1]==nn)==FALSE){
X<-X[-x,]}}
Bp<-rbind(Bp,X)}}}}
PBIB<-Bp
T <- PBIB[1, 1]
R <- length(which(T == PBIB))
lamda[1] <- l * n * (n - 1) * choose(l - 2, s - 3)*(w-1)
lamda[2] <- n * (choose(l, s - 1) + (l * choose(l - 1,s - 2)))*(w-1)
lamda[3] <- n * l * choose(l - 2, s - 3)*(w-1)
lamda[4] <- 4 * (n - 1) * choose(l - 1, s - 2)
return(list(PBIB = PBIB, Type = "Generalized rectangular right angular (4) (GPBIB_4) design (lamda not equal to 0)", V = w * V, B = dim(PBIB)[1], R = R, K = 2 * s, lamda = lamda, Resolvable=bbo))
}
