#' Peak detection using the machine learning approach.
#' 
#' The procedure uses information of known metabolites, and constructs
#' prediction models to differentiate EICs.
#' 
#' The subroutine takes CDF, mxml etc LC/MS profile.  First the profile is
#' sliced into EICs using adaptive binning. Then data features are extracted
#' from each EIC. The EICs are classified into two groups: those that have m/z
#' values that match to known m/z values, and those that don’t. Classification
#' models are built to separate the two classes, and each EIC is given a score
#' by the classification model. Those with better scores are selected to enter
#' the feature quantification step.
#' 
#' @param filename The cdf file name. If the file is not in the working
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
#' @param baseline.correct After grouping the observations, the highest
#' intensity in each group is found. If the highest is lower than this value,
#' the entire group will be deleted. The default value is NA, in which case the
#' program uses the 75th percentile of the height of the noise groups.
#' @param ridge.smoother.window The size of the smoother window used by the
#' kernel smoother to remove long ridge noise from each EIC.
#' @param smoother.window The smoother windows to use in data feature
#' generation.
#' @param known.mz The m/z values of the known metabolites.
#' @param match.tol.ppm The ppm tolerance to match identified features to known
#' metabolites/features.
#' @param do.plot Whether to produce diagnostic plots.
#' @param pos.confidence The confidence level for the features matched to the
#' known feature list.
#' @param neg.confidence The confidence level for the features not matching to
#' the known feature list.
#' @param max.ftrs.to.use The maximum number of data features to use in a
#' predictive model.
#' @param do.grp.reduce Whether to reduce data features that are similar. It is
#' based on data feature predictability.
#' @param remove.bottom.ftrs The number of worst performing data features to
#' remove before model building.
#' @param max.fpr The proportion of unmatched features to be selected in the
#' feature detection step.
#' @param min.tpr The proportion of matched features to be selected in the
#' feature detection step.
#' @param intensity.weighted Whether to weight the local density by signal
#' intensities.
#' @return A matrix with four columns: m/z value, retention time, intensity,
#' and group number.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @keywords models
#' @examples
#' 
learn.cdf <- function(filename, tol=2e-5, min.run=4, min.pres=0.3, baseline.correct=0, ridge.smoother.window=50, smoother.window=c(1, 5, 10), known.mz, match.tol.ppm=5, do.plot=FALSE, pos.confidence=0.99, neg.confidence=0.99, max.ftrs.to.use=10, do.grp.reduce=TRUE, remove.bottom.ftrs=0, max.fpr=seq(0,0.6, by=0.1), min.tpr=seq(0.8, 1,by=0.1), intensity.weighted=FALSE)
{
    if(do.plot) par(mfrow=c(2,2))

    this.name<-paste(strsplit(tolower(filename),"\\.")[[1]][1], "_", tol,"_",ridge.smoother.window, "_",baseline.correct, ".rawlearn",sep="")
    all.files<-dir()
    is.there<-all.files[which(all.files == this.name)]

    if(length(is.there) > 0)
    {
        load(this.name)
        if(do.plot)
        {
            plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",main="tolerance level loaded",axes=FALSE)
            text(x=0,y=0,tol,cex=1.2)
        }
    }else{
        this<-load.lcms(filename)
        raw.prof<-adaptive.bin.2(this, tol=tol, ridge.smoother.window=ridge.smoother.window, baseline.correct=baseline.correct, weighted=intensity.weighted)
        save(raw.prof, file=this.name)
    }

    if(is.na(baseline.correct)) baseline.correct<-0
    run.sel<-raw.prof$height.rec[which(raw.prof$height.rec[,2] >= min.run*min.pres & raw.prof$height.rec[,3]>baseline.correct),1]
    newprof<-cbind(raw.prof$masses, raw.prof$labels, raw.prof$intensi, raw.prof$grps)
    newprof<-newprof[newprof[,4] %in% run.sel,]
    raw.prof$height.rec<-raw.prof$height.rec[raw.prof$height.rec[,1] %in% newprof[,4],]

    new.prof<-cont.index(newprof,  min.pres=min.pres, min.run=min.run)
    raw.prof$height.rec<-raw.prof$height.rec[raw.prof$height.rec[,1] %in%  new.prof$new.rec[,4],]
    raw.prof$masses<-new.prof$new.rec[,1]
    raw.prof$labels<-new.prof$new.rec[,2]
    raw.prof$intensi<-new.prof$new.rec[,3]
    raw.prof$grps<-new.prof$new.rec[,4]
    raw.prof$height.rec[,3]<-new.prof$height.rec

    eic.rec.0<-eic.disect(raw.prof, smoother.window=smoother.window)
    prof<-eic.rec.0$prof
    eic.rec.0<-eic.rec.0$eic.ftrs

    a<-eic.pred(eic.rec.0, known.mz,to.use=max.ftrs.to.use, match.tol.ppm=match.tol.ppm, do.grp.reduce=do.grp.reduce, remove.bottom=remove.bottom.ftrs, max.fpr=max.fpr[1], min.tpr=min.tpr[1], do.plot=do.plot)


    fpr.tpr.combo<-expand.grid(max.fpr, min.tpr)
    prof.sel<-new("list")

    for(i in 1:nrow(fpr.tpr.combo))
    {
        chosen<-a$matched*0
        chosen[a$matched==1 & a$tpr<=fpr.tpr.combo[i,2]]<-1
        chosen[a$matched==0 & a$fpr<=fpr.tpr.combo[i,1]]<-1

        selected<-which(chosen == 1)
        prof.sel[[i]]<-prof[which(prof[,4] %in% raw.prof$height.rec[selected,1]),]
        names(prof.sel)[[i]]<-paste("fpr", fpr.tpr.combo[i,1], "tpr", fpr.tpr.combo[i,2],sep="_")
    }
    prof.sel$fpr.tpr.combo<-fpr.tpr.combo
    prof.sel$model.detail<-a

    return(prof.sel)
}
