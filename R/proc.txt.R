#' Filter noise and detect peaks from LC/MS data in text format
#'
#' This function applies the run filter to remove noise. Data points are
#' grouped into EICs in this step.
#'
#' The columns in the text file need to be separated by tab. The first column
#' needs to be the retention time, the second column the m/z values, and the
#' third column the intensity values. The first row needs to be the column
#' labels, rather than values.  The m/z are grouped based on the tolerance
#' level using multi-stage smoothing and peak finding. Non-parametric density
#' estimation is used in both m/z dimension and elution time dimension to
#' fine-tune the signal grouping. A run filter is applied, which requires a
#' "true peak" to have a minimum length in the retention time dimension
#' (parameter: min.run), as well as being detected at or higher than a
#' proportion of the time points within the time period (parameter: min.pres).
#'
#' @param filename The text file name. If the file is not in the working
#' directory, the path needs to be given.
#' @param min.pres Run filter parameter. The minimum proportion of presence in
#' the time period for a series of signals grouped by m/z to be considered a
#' peak.
#' @param min.run Run filter parameter. The minimum length of elution time for
#' a series of signals grouped by m/z to be considered a peak.
#' @param tol m/z tolerance level for the grouping of data points. This value
#' is expressed as the fraction of the m/z value. This value, multiplied by the
#' m/z value, becomes the cutoff level. The recommended value is the machine's
#' nominal accuracy level. Divide the ppm value by 1e6. For FTMS, 1e-5 is
#' recommended.
#' @param find.tol.maxd maximum distance between datapoints that are allowed in
#' the procedure to find tolerance.
#' @param baseline.correct After grouping the observations, the highest
#' intensity in each group is found. If the highest is lower than this value,
#' the entire group will be deleted. The default value is NA, in which case the
#' program uses the 75th percentile of the height of the noise groups.
#' @param baseline.correct.noise.percentile The perenctile of signal strength
#' of those EIC that don't pass the run filter, to be used as the baseline
#' threshold of signal strength.
#' @return A matrix with four columns: m/z value, retention time, intensity,
#' and group number.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @keywords models
#' @export
proc.txt <-
function(filename, min.pres=0.5, min.run=12,tol=NA, find.tol.maxd=1e-4, baseline.correct.noise.percentile=0.25, baseline.correct=0)
{
    par(mfrow=c(2,2))
    this.name<-paste(strsplit(tolower(filename),".txt")[[1]],"_",tol,"_",find.tol.maxd,".rawprof",sep="")
    all.files<-dir()
    is.there<-all.files[which(all.files == this.name)]

    if(length(is.there) > 0)
    {
        load(this.name)
        plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",main="tolerance level loaded",axes=FALSE)
        text(x=0,y=0,tol,cex=1.2)
    }else{
        ex.nc<-read.table(filename, header=TRUE, sep="\t")
        ex.nc<-as.matrix(ex.nc)
        ex.nc<-ex.nc[ex.nc[,3] != 0,]
        ex.nc<-ex.nc[order(ex.nc[,1], ex.nc[,2]),]

        this<-new("list")
        this$masses<-ex.nc[,2]
        this$intensi<-ex.nc[,3]
        this$labels<-ex.nc[,1]
        this$times<-unique(ex.nc[,1])

        raw.prof<-adaptive.bin(this, min.run=min.run, min.pres=min.pres, tol=tol, baseline.correct=baseline.correct)
        save(raw.prof, file=this.name)
    }

    newprof<-cbind(raw.prof$masses, raw.prof$labels, raw.prof$intensi, raw.prof$grps)
    h.1<-log10(raw.prof$height.rec[raw.prof$height.rec[,2]<= max(2, raw.prof$min.count.run*min.pres/2),3])
    h.2<-log10(raw.prof$height.rec[raw.prof$height.rec[,2]>= raw.prof$min.count.run*min.pres,3])
    if(length(h.1)>50)
    {
        plot(density(h.1),xlab="maximum height of group (log scale)", xlim=range(c(h.1, h.2)),main="Black - noise groups \n Blue - selected groups")
    }else{
        plot(NA,NA,xlab="maximum height of group (log scale)", xlim=range(c(h.1, h.2)), ylim=c(0,1), main="Black - noise groups \n Blue - selected groups")
        if(length(h.1)>0)	abline(v=h.1)
    }

    if(is.na(baseline.correct))
    {
        baseline.correct<-10^quantile(h.1, baseline.correct.noise.percentile)
        message(c("maximal height cut is automatically set at the", baseline.correct.noise.percentile, "percentile of noise group heights: ", baseline.correct))
    }else{
        message(c("maximal height cut is provided by user: ", baseline.correct))
    }
    abline(v=log10(baseline.correct), col="red")

    if(is.na(baseline.correct)) baseline.correct<-0
    run.sel<-raw.prof$height.rec[which(raw.prof$height.rec[,2] >= raw.prof$min.count.run*min.pres & raw.prof$height.rec[,3]>baseline.correct),1]
    newprof<-newprof[newprof[,4] %in% run.sel,]
    new.prof<-cont.index(newprof,  min.pres=min.pres, min.run=min.run)
    lines(density(log10(new.prof$height.rec)),col="blue")
    time.range.rec<-new.prof$time.range.rec
    mz.pres.rec<-new.prof$mz.pres.rec
    hist(time.range.rec, xlab="Range of retention time in the same group", ylab="Density",freq=FALSE,nclass=100,main="Group retention time range distribution")
    hist(mz.pres.rec, xlab="% signal present in the same group", ylab="Density" ,freq=FALSE, nclass=20, main="Group % present signal distribution")

    return(new.prof$new.rec)
}
