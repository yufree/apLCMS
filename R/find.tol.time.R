#' An internal function that is not supposed to be directly accessed by the
#' user. Find elution time tolerance level.
#' 
#' This function finds the time tolerance level. Also, it returns the grouping
#' information given the time tolerance.
#' 
#' The peaks are first ordered by m/z, and split into groups by the m/z
#' tolerance. Then within every peak group, the pairwise elution time
#' difference is calculated. All the pairwise elution time differences within
#' groups are merged into a single vector. A mixture model (unknown
#' distribution for distance between peaks from the same feature, and a
#' triangle-shaped distribution for distance between peaks from different
#' features) is fit to find the elution time tolerance level. The elution times
#' within each peak group are then ordered. If a gap between consecutive
#' retention times is larger than the elution time tolerance level, the group
#' is further split at the gap. Grouping information is returned, as well as
#' the elution time tolerance level.
#' 
#' @param mz mz value of all peaks in all profiles in the study.
#' @param chr retention time of all peaks in all profiles in the study.
#' @param lab label of all peaks in all profiles in the study.
#' @param num.exp The number of spectra in this analysis.
#' @param mz.tol m/z tolerance level for the grouping of signals into peaks.
#' This value is expressed as the percentage of the m/z value. This value,
#' multiplied by the m/z value, becomes the cutoff level.
#' @param chr.tol the elution time tolerance. If NA, the function finds the
#' tolerance level first. If a numerical value is given, the function directly
#' goes to the second step - grouping peaks based on the tolerance.
#' @param aver.bin.size The average bin size to determine the number of equally
#' spaced points in the kernel density estimation.
#' @param min.bins the minimum number of bins to use in the kernel density
#' estimation. It overrides aver.bin.size when too few observations are
#' present.
#' @param max.bins the maximum number of bins to use in the kernel density
#' estimation. It overrides aver.bin.size when too many observations are
#' present.
#' @param max.mz.diff As the m/z tolerance in alignment is expressed in
#' relative terms (ppm), it may not be suitable when the m/z range is wide.
#' This parameter limits the tolerance in absolute terms. It mostly influences
#' feature matching in higher m/z range.
#' @param max.num.segments the maximum number of segments.
#' @return A list object is returned: \item{chr.tol}{ The elution time
#' tolerance level.} \item{comp2 }{A matrix with six columns. Every row
#' corrsponds to a peak in one of the spectrum. The columns are: m/z, elution
#' time, spread, signal strength, spectrum label, and peak group label. The
#' rows are ordered by the median m/z of each peak group, and with each peak
#' group the rows are ordered by the elution time.}
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @keywords models
#' @examples
#' 
#' 
#' 
find.tol.time <-
function(mz, chr, lab, num.exp, mz.tol=2e-5, chr.tol=NA, aver.bin.size=200, min.bins=50, max.bins=100, max.mz.diff=0.01, max.num.segments=10000)
{
    o<-order(mz)
    mz<-mz[o]
    chr<-chr[o]
    lab<-lab[o]
    rm(o)
    
    l<-length(mz)
    
    breaks<-c(0, which((mz[2:l]-mz[1:(l-1)]) > min(max.mz.diff, mz.tol*((mz[2:l]+mz[1:(l-1)])/2))), l)
    
    for(i in 2:length(breaks))
    {
        this.o<-order(chr[(breaks[i-1]+1):breaks[i]])
        this.o<-this.o+breaks[i-1]
        mz[(breaks[i-1]+1):breaks[i]]<-mz[this.o]
        chr[(breaks[i-1]+1):breaks[i]]<-chr[this.o]
        lab[(breaks[i-1]+1):breaks[i]]<-lab[this.o]
    }
    
    breaks<-breaks[c(-1,-length(breaks))]
    if(is.na(chr.tol))
    {
        da<-0
        if(length(breaks)>max.num.segments)
        {
            s<-floor(seq(2, length(breaks), length.out=max.num.segments))
        }else{
            s<-2:length(breaks)
        }
        
        for(i in s)
        {
            this.sel<-(breaks[i-1]+1):breaks[i]
            
            if(length(this.sel) <= 3*num.exp)
            {
                this.d<-as.vector(dist(chr[this.sel]))
                if(length(this.d)>100) this.d<-sample(this.d,100)
                da<-c(da,this.d)
            }
        }
        
        da<-da[!is.na(da)]
        uppermost<-max(da)
        n=min(max.bins,max(min.bins, round(length(da)/aver.bin.size)))
        des<-density(da,kernel="gaussian",n=n, bw=uppermost/n*2,from=0)
        y<-des$y[des$x>0]
        x<-des$x[des$x>0]
        
        this.l<-lm(y[x>uppermost/4]~x[x>uppermost/4])
        
        exp.y<-this.l$coef[1]+this.l$coef[2]*x
        
        plot(x,y,main="find retention time tolerance", xlab="Delta", ylab="Density",cex=0.25)
        lines(x,exp.y,col="red")
        y2<-y[1:(length(y)-1)]
        y3<-y[2:(length(y))]
        y2[which(y2<y3)]<-y3[which(y2<y3)]
        y[1:(length(y)-1)]<-y2
        
        yy<-cumsum(y>1.5*exp.y)
        yi<-1:length(yy)
        sel<-min(which(yy<yi))-1
        
        abline(v=x[sel],col="blue")
        chr.tol<-x[sel]
    }
    
    da<-chr[2:l]-chr[1:(l-1)]
    breaks.2<-which(da>chr.tol)
    all.breaks<-c(0, unique(c(breaks, breaks.2)), l)
    all.breaks<-all.breaks[order(all.breaks)]
    
    grps<-1:length(mz)
    for(i in 2:length(all.breaks))
    {
        grps[(all.breaks[i-1]+1):all.breaks[i]]<-i
    }
    
    to.return<-list(mz=mz, chr=chr, lab=lab, grps=grps, chr.tol=chr.tol, mz.tol=mz.tol)
    return(to.return)
}
