#' U-type design via some PBIB designs
#'
#' Applies the Fang algorithm on our constructed designs to obtain
#' the configuration and the parameters of the associated U-type design.
#' @usage UType(lst)
#' @param lst The output of one of our package functions.
#' @return A LIST :
#'  \itemize{
#'   \item \code{v} Number of runs.
#'   \item \code{r } Number of factors.
#'   \item \code{UtypeDesign } The configuration of the U-type design..
#'   }
#'
#' @author Mohamed Laib, Imane Rezgui, Zebida Gheribi-Aoulmi and Herve Monod
#' @references
#'K.T. Fang, R.Li and A.Sudjanto (2006). Design ans Modeling for Computer
#'Experiments. Taylor & Francis Group, LLC London.
#'
#' @importFrom stats na.omit
#' @examples
#' \dontrun{
#' M<-GPBIB4A(4,4,2,2)
#' UType(M)
#' }
#' @export

UType <-function (lst) {
bob<-lst$Resolvable

if (bob=="Yes") {
mat<-lst$PBIB
    w <- mat
    v <- sort(unique(as.vector(w)))
    W <- w
    malist <- NULL
    TRI <- function(a, b) all(a %in% b)
    k <- 1
    while (TRUE) {
        vv <- v
        malist[[k]] <- c(0, 0)
        it <- 0
        while (TRUE) {
            it <- it + 1
            bool <- apply(W, 1, TRI, vv)
            if (!any(bool)) {
                break
            }
            else {
                ind <- which(bool)[1]
                u <- W[ind, ]
                vv <- vv[!(vv %in% u)]
                W <- W[-ind, , drop = FALSE]
                malist[[k]] <- rbind(malist[[k]], u)
            }
        }
        malist[[k]] <- malist[[k]][-1, ]
        if (!all(as.vector(W) %in% v) | length(W) == 0)
            break
        k <- k + 1
    }
    x <- Reduce("rbind", malist)
    a <- max(x)
    b <- length(malist)
    c <- dim(x)[1]
    lev <- c/b
    ex <- length(v)
    UD <- matrix(nrow = a, ncol = b)
    v <- sort(unique(as.vector(x)))
    for (i in 1:b) {
        q <- malist[[i]]
        e <- c()
        for (j in v) {
            e[j] <- which(q == j)
            if (e[j] > lev) {
                e[j] <- (e[j]%%lev)
            }
            if (e[j] == 0) {
                e[j] <- lev
            }
        }
        UD[, i] <- e
    }
    ud <- na.omit(UD)
    N <- dim(ud)[1]
    F <- dim(ud)[2]
    return(list(v = N, r = F, UtypeDesign = ud))
}
else {"The PBIB is not resolvable"}}
