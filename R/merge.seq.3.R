#' An internal function.
#'
#' This is a internal function. It shouldn't be called by the end user.
#'
#'
#' @param a vector of retention time.
#' @param mz vector of m/z ratio.
#' @param inte vector of signal strength.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @references Bioinformatics. 25(15):1930-36.  BMC Bioinformatics. 11:559.
#' @keywords models
merge.seq.3 <-
function(a, mz, inte)             ### the input need to be pre-ordered by a
{
    l <- length(a)
    breaks <- c(0, which(a[1:(l - 1)] != a[2:l]), l)
    new.int <- new.mz <- rep(0, length(breaks)-1)

    for (i in 1:(length(breaks) - 1)) {
        this.int<-inte[(breaks[i] + 1):breaks[i + 1]]
        this.mz<-mz[(breaks[i] + 1):breaks[i + 1]]
        new.int[i] <- sum(this.int)
        new.mz[i] <- median(this.mz[which(this.int==max(this.int))])
    }
    new.a <- unique(a)
    return(cbind(new.mz, new.a, new.int))
}
