#' Generalized rectangular right angular (4) design with \eqn{\lambda_4}{lambda4} = 0
#'
#' Gives the configuration and the parametres of the design obtained by
#' the first construction method of \code{GPBIB_4} (see 3.1.1 of the paper
#' rezgui et al (2015)).
#' @usage GPBIB4A(n, l, s, w)
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
#' @details \itemize{
#'          \item For \eqn{s = l}, the previous method gives configuration of nested group divisible designs.
#'  }
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
#' @seealso \code{\link{GPBIB4B}} and \code{\link{UType}}
#' @note For \eqn{w=2}, the \code{GPBIB_4} is a rectangular right angular (4) (PBIB_4)
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' n<-3
#' l<-3
#' s<-3
#' w<-3
#' GPBIB4A(n, l, s, w)
#' }
#' @export
GPBIB4A <-function(n,l,s,w){
if (s<2 & l<2 & n<2) stop("n,l,s should be great than 1")
V<-n*l
reso<-(n*l)%%(2*s)
bbo<-ifelse(reso==0,"Yes","No")

A<-NULL;mat<-NULL;lamda<-NULL
for (i in 1:w){
A[[i]]<-matrix(1:V, ncol=l, byrow=TRUE)
z<-(i-1)*V
A[[i]]<-A[[i]]+z}


for (i in 1:w){
mw<-A[[i]]
mat[[i]]<-Opn(mw,s)}

BIB<-Reduce("rbind",mat)
T<-BIB[1,1];R<-length(which(T==BIB))
lamda[1]<-(n-1)*choose(l-2,s-2);lamda[2]<-choose(l-1,s-1);lamda[3]<-choose(l-2,s-2);lamda[4]<-0
return(list(PBIB=BIB,Type="Generalized rectangular right angular (4) (GPBIB_4) design",V=w*V,B=dim(BIB)[1],R=R,K=2*s,lamda=lamda, Resolvable=bbo))}


Opn<-function(mt,s){
  n<-dim(mt)[1];l<-dim(mt)[2];V<-n*l
  a<-combn(1:l, s);b<-combn(1:n, 2)
  v<-dim(b)[2];vv<-dim(a)[2]
  MAT<-NULL;A<-1;y<-1
  while(A<=vv) {
    for (k in 1:v){
      s<-a[,A];ss<-b[,k]
      MAT[[y]]<-as.vector(t(mt[ss,s]))
      y<-y+1}
    A<-A+1}
  return(Reduce("rbind",MAT))}


