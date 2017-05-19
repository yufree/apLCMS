#' Plot the data in the m/z and retention time plane.
#' 
#' This is a diagnostic function. It takes the original CDF file, as well as
#' the detected feature table, and plots the data in the m/z - retention time
#' plane, using a user-defined range. The entire data is too big to plot, thus
#' the main purpose is to focus on small subregions of the data and check the
#' peak detection results.
#' 
#' 
#' @param rawname The CDF file name.
#' @param f The output object of prof.to.feature().
#' @param mzlim The m/z range to plot.
#' @param timelim The retention time range to plot.
#' @param lwd Line width parameter, to be passed on to the function line().
#' @return There is no return value.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @references Bioinformatics. 25(15):1930-36.  BMC Bioinformatics. 11:559.
#' @keywords models
#' @examples
#' 
plot.cdf.2d <-
function(rawname, f, mzlim,timelim, lwd=1)
### rawname is the cdf file name
### f is the output object of prof.to.features()
### lwd is line width
{
    ########### read the raw table
    this.raw<-load.lcms(rawname)
    masses<-this.raw$masses
    intensi<-this.raw$intensi
    labels<-this.raw$labels
    times<-this.raw$times
    rm(this.raw)
    
    times<-times[order(times)]
    base.curve<-unique(times)
    base.curve<-base.curve[order(base.curve)]
    base.curve<-cbind(base.curve, base.curve*0)
    
    curr.order<-order(masses)
    intensi<-intensi[curr.order]
    labels<-labels[curr.order]
    masses<-masses[curr.order]
    
    
    #####################
    sel<-which(labels>=timelim[1] & labels<=timelim[2] & masses>=mzlim[1] & masses <=mzlim[2])
    ex.nc<-cbind(labels[sel], masses[sel])
    plot(ex.nc[,1:2],cex=.1)
    
    f<-f[f[,1]>=mzlim[1] & f[,1]<=mzlim[2],]
    for(i in 1:nrow(f))
    {
        this.mz<-f[i,1]
        this.time<-c(f[i,2]-f[i,3], f[i,2]+f[i,3])
        lines(this.time,rep(this.mz,2),col="red",lwd=lwd)
    }
}
