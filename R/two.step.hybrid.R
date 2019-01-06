#' Two step hybrid feature detection.
#'
#' A two-stage hybrid feature detection and alignment procedure, for data
#' generated in multiple batches.
#'
#' The function first conducts hybrid feature detection and alignment in each
#' batch separately. Then a between-batch RT correction and feature alignment
#' is conducted. Weak signal recovery is conducted at the single feature table
#' level.
#'
#' @param folder The folder where all CDF files to be processed are located.
#' For example “C:/CDF/this_experiment”
#' @param info A table with two columns. The first column is the file names,
#' and the second column is the batch label of each file.
#' @param min.within.batch.prop.detect A feature needs to be present in at
#' least this proportion of the files, for it to be initially detected as a
#' feature for a batch. This parameter replaces the "min.exp" parameter in
#' semi.sup().
#' @param min.within.batch.prop.report A feature needs to be present in at
#' least this proportion of the files, in a proportion of batches controlled by
#' "min.batch.prop", to be included in the final feature table. This parameter
#' replaces the "min.exp" parameter in semi.sup().
#' @param min.batch.prop A feature needs to be present in at least this
#' proportion of the batches, for it to be considered in the entire data.
#' @param batch.align.mz.tol The m/z tolerance in ppm for between-batch
#' alignment.
#' @param batch.align.chr.tol The RT tolerance for between-batch alignment.
#' @param file.pattern The pattern in the names of the files to be processed.
#' The default is ".cdf". Other formats supported by mzR package can also be
#' used, e.g. "mzML" etc.
#' @param known.table A data frame containing the known metabolite ions and
#' previously found features. It contains 18 columns: "chemical_formula": the
#' chemical formula if knonw; "HMDB_ID": HMDB ID if known; "KEGG_compound_ID":
#' KEGG compound ID if known; "neutral.mass": the neutral mass if known:
#' "ion.type": the ion form, such as H+, Na+, ..., if known; "m.z": m/z value,
#' either theoretical for known metabolites, or mean observed value for unknown
#' but previously found features; "Number_profiles_processed": the total number
#' of LC/MS profiles that were used to build this database; "Percent_found": in
#' what percentage was this feature found historically amount all data
#' processed in building this database; "mz_min": the minimum m/z value
#' observed for this feature; "mz_max": the maximum m/z value observed for this
#' feature; "RT_mean": the mean retention time observed for this feature;
#' "RT_sd": the standard deviation of retention time observed for this feature;
#' "RT_min": the minimum retention time observed for this feature; "RT_max":
#' the maximum retention time observed for this feature; "int_mean.log.": the
#' mean log intensity observed for this feature; "int_sd.log.": the standard
#' deviation of log intensity observed for this feature; "int_min.log.": the
#' minimum log intensity observed for this feature; "int_max.log.": the maximum
#' log intensity observed for this feature;
#' @param n.nodes The number of CPU cores to be used through doSNOW.
#' @param min.pres This is a parameter of thr run filter, to be passed to the
#' function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param min.run This is a parameter of thr run filter, to be passed to the
#' function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param mz.tol The user can provide the m/z tolerance level for peak
#' identification. This value is expressed as the percentage of the m/z value.
#' This value, multiplied by the m/z value, becomes the cutoff level. Please
#' see the help for proc.cdf() for details.
#' @param baseline.correct.noise.percentile The perenctile of signal strength
#' of those EIC that don't pass the run filter, to be used as the baseline
#' threshold of signal strength. This parameter is passed to proc.cdf()
#' @param shape.model The mathematical model for the shape of a peak. There are
#' two choices - “bi-Gaussian” and “Gaussian”. When the peaks are asymmetric,
#' the bi-Gaussian is better. The default is “bi-Gaussian”.
#' @param baseline.correct This is a parameter in peak detection. After
#' grouping the observations, the highest observation in each group is found.
#' If the highest is lower than this value, the entire group will be deleted.
#' The default value is NA, which allows the program to search for the cutoff
#' level. Please see the help for proc.cdf() for details.
#' @param peak.estim.method the bi-Gaussian peak parameter estimation method,
#' to be passed to subroutine prof.to.features. Two possible values: moment and
#' EM.
#' @param min.bw The minimum bandwidth in the smoother in prof.to.features().
#' Please see the help file for prof.to.features() for details.
#' @param max.bw The maximum bandwidth in the smoother in prof.to.features().
#' Please see the help file for prof.to.features() for details.
#' @param sd.cut A parameter for the prof.to.features() function. A vector of
#' two. Features with standard deviation outside the range defined by the two
#' numbers are eliminated.
#' @param sigma.ratio.lim A parameter for the prof.to.features() function. A
#' vector of two. It enforces the belief of the range of the ratio between the
#' left-standard deviation and the righ-standard deviation of the bi-Gaussian
#' fuction used to fit the data.
#' @param component.eliminate In fitting mixture of bi-Gaussian (or Gaussian)
#' model of an EIC, when a component accounts for a proportion of intensities
#' less than this value, the component will be ignored.
#' @param moment.power The power parameter for data transformation when fitting
#' the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param align.chr.tol The user can provide the elution time tolerance level
#' to override the program’s selection. This value is in the same unit as the
#' elution time, normaly seconds. Please see the help for match.time() for
#' details.
#' @param align.mz.tol The user can provide the m/z tolerance level for peak
#' alignment to override the program’s selection.  This value is expressed as
#' the percentage of the m/z value. This value, multiplied by the m/z value,
#' becomes the cutoff level.Please see the help for feature.align() for
#' details.
#' @param max.align.mz.diff As the m/z tolerance in alignment is expressed in
#' relative terms (ppm), it may not be suitable when the m/z range is wide.
#' This parameter limits the tolerance in absolute terms. It mostly influences
#' feature matching in higher m/z range.
#' @param pre.process Logical. If true, the program will not perform time
#' correction and alignment. It will only generate peak tables for each spectra
#' and save the files. It allows manually dividing the task to multiple
#' machines.
#' @param recover.mz.range A parameter of the recover.weaker() function. The
#' m/z around the feature m/z to search for observations. The default value is
#' NA, in which case 1.5 times the m/z tolerance in the aligned object will be
#' used.
#' @param recover.chr.range A parameter of the recover.weaker() function. The
#' retention time around the feature retention time to search for observations.
#' The default value is NA, in which case 0.5 times the retention time
#' tolerance in the aligned object will be used.
#' @param use.observed.range A parameter of the recover.weaker() function. If
#' the value is TRUE, the actual range of the observed locations of the feature
#' in all the spectra will be used.
#' @param match.tol.ppm The ppm tolerance to match identified features to known
#' metabolites/features.
#' @param new.feature.min.count The number of profiles a new feature must be
#' present for it to be added to the database.
#' @param recover.min.count The minimum time point count for a series of point
#' in the EIC for it to be considered a true feature.
#' @return A list is returned.  \item{batchwise.results}{A list. Each item in
#' the list is the product of semi.sup() from a single batch.}
#' \item{final.ftrs}{Feature table. This is the end product of the function.}
#' @author Tianwei Yu < tianwei.yu@@emory.edu>
#' @seealso semi.sup, cdf.to.ftrs, proc.cdf, prof.to.feature, adjust.time,
#' feature.align, recover.weaker
#' @keywords models
#' @export
two.step.hybrid <-function(folder, info, min.within.batch.prop.detect=0.1, min.within.batch.prop.report=0.5, min.batch.prop=0.5, batch.align.mz.tol=1e-5, batch.align.chr.tol=50, file.pattern=".cdf", known.table=NA, n.nodes=4, min.pres=0.5, min.run=12, mz.tol=1e-5, baseline.correct.noise.percentile=0.05, shape.model="bi-Gaussian",  baseline.correct=0, peak.estim.method="moment", min.bw=NA, max.bw=NA, sd.cut=c(0.1,100), sigma.ratio.lim=c(0.05, 20), component.eliminate=0.01, moment.power=2, align.mz.tol=NA, align.chr.tol=NA, max.align.mz.diff=0.01, pre.process=FALSE, recover.mz.range=NA, recover.chr.range=NA, use.observed.range=TRUE, match.tol.ppm=NA, new.feature.min.count=2, recover.min.count=3)
{
    setwd(folder)
    info<-as.matrix(as.data.frame(info))
    batches<-unique(info[,2])

    batchwise<-new("list")
    message("total number of batches: ", length(batches))

    for(batch.i in 1:length(batches))
    {
        message("working on batch number ", batch.i)
        batch<-batches[batch.i]
        files<-dir(path = ".", pattern=file.pattern, ignore.case = TRUE)
        this.subs<-which(files %in% info[info[,2]==batch,1])
		message("total number of files in this batch ", length(this.subs))

        b<-semi.sup(folder,n.nodes=n.nodes, subs=this.subs, file.pattern=file.pattern,known.table=known.table, sd.cut=sd.cut,sigma.ratio.lim=sigma.ratio.lim, component.eliminate=component.eliminate, moment.power=moment.power, min.pres=min.pres, min.run=min.run, min.exp=ceiling(min.within.batch.prop.detect*length(this.subs)), mz.tol=mz.tol, baseline.correct.noise.percentile=baseline.correct.noise.percentile, align.mz.tol=align.mz.tol, align.chr.tol=align.chr.tol, max.align.mz.diff=max.align.mz.diff, recover.mz.range=recover.mz.range, recover.chr.range=recover.chr.range,use.observed.range=use.observed.range, shape.model=shape.model,new.feature.min.count=new.feature.min.count, recover.min.count=recover.min.count)

        b$final.ftrs<-b$final.ftrs[order(b$final.ftrs[,1], b$final.ftrs[,2]),]
        b$final.times<-b$final.times[order(b$final.times[,1], b$final.times[,2]),]

        batchwise[[batch.i]]<-b
    }

    fake.features<-new("list")
    for(batch.i in 1:length(batches))
    {
        this.fake<-batchwise[[batch.i]]$final.ftrs
        this.fake[,3:4]<-NA
        this.fake[,5]<-apply(this.fake[,-1:-4],1,median)
        this.fake<-this.fake[,1:5]
        fake.features[[batch.i]]<-this.fake
    }

    cl <- makeCluster(n.nodes,type='SOCK')
    registerDoParallel(cl)
    #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
    clusterEvalQ(cl, library(apLCMS))

    fake2<-adjust.time(fake.features,mz.tol=batch.align.mz.tol, chr.tol=batch.align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)

    message("Alignment across batches")
    fake3<-feature.align(fake2, min.exp=ceiling(min.batch.prop*length(batches)), mz.tol=batch.align.mz.tol,chr.tol=batch.align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)

    stopCluster(cl)

    message("Recovery across batches")

    cl <- makeCluster(n.nodes,type='SOCK')
    registerDoParallel(cl)
    #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
    clusterEvalQ(cl, library(apLCMS))

    aligned<-foreach(batch.i=1:length(batches), .combine=cbind) %dopar%
    {
        this.fake<-batchwise[[batch.i]]$final.ftrs
        this.features<-batchwise[[batch.i]]$features2
        this.medians<-apply(this.fake[,-1:-4],1,median)

        orig.time<-this.fake[,2]
        adjusted.time<-fake2[[batch.i]][,2]

        this.aligned<-matrix(0, nrow=nrow(fake3$aligned.ftrs), ncol=ncol(this.fake)-4)

        # adjusting the time (already within batch adjusted)

        for(j in 1:length(this.features))
        {
            for(i in 1:nrow(this.features[[j]]))
            {
                diff.time<- abs(orig.time - this.features[[j]][i,2])
                sel<-which(diff.time == min(diff.time))[1]
                this.features[[j]][i,2] <- adjusted.time[sel]
            }
        }

        for(i in 1:nrow(this.aligned))
        {
            if(fake3$aligned[i, batch.i+4] != 0)
            {
                sel<-which(fake3$aligned[i,3]<=this.fake[,1] & fake3$aligned[i,4]>=this.fake[,1] & abs(this.medians - fake3$aligned[i, batch.i+4]) < 1)
                if(length(sel) == 0) sel<-which(fake3$aligned[i,3]<=this.fake[,1] & fake3$aligned[i,4]>=this.fake[,1])
                if(length(sel) == 0)
                {
                    message("batch", batch.i, " row ", i, " match issue.")
                }else{
                    if(length(sel) == 1)
                    {
                        this.aligned[i, ]<-this.fake[sel,-1:-4]
                    }else{
                        this.aligned[i, ]<-apply(this.fake[sel,-1:-4],2,sum)
                    }
                }
            }else{
                ### go into individual feature tables to find a match
                recaptured<-rep(0, ncol(this.aligned))
                for(j in 1:length(this.features))
                {
                    diff.mz<-abs(this.features[[j]][,1]-fake3$aligned[i,1])
                    diff.time<-abs(this.features[[j]][,2]-fake3$aligned[i,2])
                    sel<-which(diff.mz < fake3$aligned[i,1]*batch.align.mz.tol & diff.time <= batch.align.chr.tol)
                    if(length(sel)>1) sel<-sel[which(diff.time[sel] == min(diff.time[sel]))[1]]

                    if(length(sel) ==1) recaptured[j]<-this.features[[j]][sel,5]
                }
                this.aligned[i, ]<-recaptured
            }
        }

        colnames(this.aligned)<-colnames(this.fake)[-1:-4]

        this.aligned
    }

    aligned<-cbind(fake3$aligned.ftrs[,1:4], aligned)

	batch.presence.mat<-matrix(0, nrow=nrow(aligned), ncol=length(batches))
	for(batch.i in 1:length(batches))
	{
		batch<-batches[batch.i]
        this.mat<-aligned[,which(colnames(aligned) %in% info[info[,2]==batch,1])]
		this.mat<- 1*(this.mat != 0)
		this.presence<-apply(this.mat,1,sum)/ncol(this.mat)
		batch.presence.mat[,batch.i]<-1*(this.presence >= min.within.batch.prop.report)
	}
	batch.presence<-apply(batch.presence.mat,1,sum)/ncol(batch.presence.mat)
	final.aligned<-aligned[which(batch.presence >= min.batch.prop),]

    stopCluster(cl)
    to.return<-new("list")
    to.return$batchwise.results<-batchwise
	to.return$all.detected.ftrs<-aligned
    to.return$final.ftrs<-final.aligned
    return(to.return)
}
