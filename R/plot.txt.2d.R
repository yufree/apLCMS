#' Plot the data in the m/z and retention time plane.
#' 
#' This is a diagnostic function. It takes the original text file, as well as
#' the detected feature table, and plots the data in the m/z - retention time
#' plane, using a user-defined range. The entire data is too big to plot, thus
#' the main purpose is to focus on small subregions of the data and check the
#' peak detection results.
#' 
#' The columns in the text file need to be separated by tab. The first column
#' needs to be the retention time, the second column the m/z values, and the
#' third column the intensity values. The first row needs to be the column
#' labels, rather than values.
#' 
#' @param rawname The text file name.
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
plot.txt.2d <-
function(rawname, f, mzlim,timelim, lwd=1)
### rawname is the text file name
### f is the output object of prof.to.features()
### lwd is line width
{
    ########### read the raw table
    ex.nc<-read.table(rawname, header=T, sep="\t")
    ex.nc<-as.matrix(ex.nc)
    ex.nc<-ex.nc[ex.nc[,3] != 0,]
    
    #####################
    
    ex.nc<-ex.nc[ex.nc[,1]>=timelim[1] & ex.nc[,1]<=timelim[2] & ex.nc[,2]>=mzlim[1] & ex.nc[,2] <=mzlim[2],]
    plot(ex.nc[,1:2],cex=.1)
    
    f<-f[f[,1]>=mzlim[1] & f[,1]<=mzlim[2],]
    for(i in 1:nrow(f))
    {
        this.mz<-f[i,1]
        this.time<-c(f[i,2]-f[i,3], f[i,2]+f[i,3])
        lines(this.time,rep(this.mz,2),col="red", lwd=lwd)
    }
}
