#' An internal function that is not supposed to be directly accessed by the
#' user. Find m/z tolerance level.
#' 
#' The function finds the tolerance level in m/z from a given vector of
#' observed m/z values.
#' 
#' The method assumes a mixture model: an unknown distribution of m/z
#' variations in the same peak, and an exponential distribution of between-peak
#' diffs. The parameter of the exponential distribution is estimated by the
#' upper 75% of the sorted data, and the cutoff is selected where the density
#' of the empirical distribution is >1.5 times the density of the exponential
#' distribution.
#' 
#' @param a The vector of observed m/z values.
#' @param uppermost Consider only m/z diffs smaller than this value.
#' @param aver.bin.size The average bin size to determine the number of equally
#' spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density
#' estimation. It overrides aver.bin.size when too few observations are
#' present.
#' @param max.bins the maximum number of bins to use in the kernel density
#' estimation. It overrides aver.bin.size when too many observations are
#' present.
#' @return The tolerance level is returned.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @keywords models
#' @examples
#' 
#' 
#' data(prof)
#' find.tol(prof[[1]][,1])
#' 
#' 
find.tol <-
function(a, uppermost=1e-4, aver.bin.size=4000, min.bins=50, max.bins=200)
{
    a<-a[order(a)]
    l<-length(a)
    da<-(a[2:l]-a[1:(l-1)])/((a[2:l]+a[1:(l-1)])/2)
    da<-da[da < uppermost]
    n<-min(max.bins, max(round(length(da)/aver.bin.size), min.bins))
    des<-density(da,kernel="gaussian",n=n, bw=uppermost/n*2,from=0)
    y<-des$y[des$x>0]
    x<-des$x[des$x>0]
    
    to.use<-da[da>max(da)/4]-max(da)/4
    library(MASS)
    this.rate<-fitdistr(to.use, "exponential")$estimate
    exp.y<-dexp(x,rate=this.rate)
    exp.y<-exp.y*sum(y[x>max(da)/4])/sum(exp.y[x>max(da)/4])
    
    plot(x,y, xlab="Delta",ylab="Density",main="find m/z tolerance",cex=.25)
    lines(x,exp.y,col="red")
    
    yy<-cumsum(y)
    y2<-cumsum(exp.y)
    yy<-cumsum(y>1.5*exp.y)
    yi<-1:length(yy)
    sel<-min(which(yy<yi))-1
    
    abline(v=x[sel],col="blue")
    return(x[sel])
    
}
