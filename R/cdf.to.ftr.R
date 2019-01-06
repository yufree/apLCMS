#' Convert a number of cdf files in the same directory to a feature table
#'
#' This is a wrapper function, which calls four other functions to convert a
#' number of cdf files to a feature table. All cdf files to be processed must
#' be in a single folder.
#'
#' The wrapper function calls five other functions to perform the feature table
#' generation. Every spectrum (cdf file) first goes through proc.cdf() and
#' prof.to.feature() to generate a spectrum-level peak table. The eluction time
#' correction is done by match.time(). Then the peaks are aligned across
#' spectra by feature.align(). For features deteced in a portion of the
#' spectra, weaker signals in other spectra are recovered by recover.weaker().
#' From version 4, the parameter mz.tol can no longer be NA. This is to allow
#' the program better process data other than FTLCMS. It is recommended that
#' the user use the machine's claimed accuracy. For FTMS, 1e-5 is recommended.
#'
#' @param folder The folder where all CDF files to be processed are located.
#' For example ÒC:/CDF/this_experimentÓ
#' @param file.pattern The pattern in the names of the files to be processed.
#' The default is ".cdf". Other formats supported by mzR package can also be
#' used, e.g. "mzML" etc.
#' @param n.nodes The number of CPU cores to be used through doSNOW.
#' @param min.exp If a feature is to be included in the final feature table, it
#' must be present in at least this number of spectra.
#' @param min.pres This is a parameter of thr run filter, to be passed to the
#' function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param min.run This is a parameter of thr run filter, to be passed to the
#' function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param subs If not all the CDF files in the folder are to be processed, the
#' user can define a subset using this parameter. For example, subs=15:30, or
#' subs=c(2,4,6,8)
#' @param mz.tol The user can provide the m/z tolerance level for peak
#' identification. This value is expressed as the percentage of the m/z value.
#' This value, multiplied by the m/z value, becomes the cutoff level. Please
#' see the help for proc.cdf() for details.
#' @param baseline.correct.noise.percentile The perenctile of signal strength
#' of those EIC that don't pass the run filter, to be used as the baseline
#' threshold of signal strength. This parameter is passed to proc.cdf()
#' @param shape.model The mathematical model for the shape of a peak. There are
#' two choices - bi-Gaussian and Gaussian. When the peaks are asymmetric, the
#' bi-Gaussian is better. The default is bi-Gaussian.
#' @param BIC.factor the factor that is multiplied on the number of parameters
#' to modify the BIC criterion. If larger than 1, models with more peaks are
#' penalized more.
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
#' to override the programÕs selection. This value is in the same unit as the
#' elution time, normaly seconds. Please see the help for match.time() for
#' details.
#' @param align.mz.tol The user can provide the m/z tolerance level for peak
#' alignment to override the programÕs selection.  This value is expressed as
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
#' @param recover.min.count The minimum time point count for a series of point
#' in the EIC for it to be considered a true feature.
#' @param intensity.weighted Whether to weight the local density by signal
#' intensities in the initial peak detection.
#' @return A list is returned.  \item{features}{ A list object, each component
#' of which being the peak table from a single spectrum.} \item{features2}{A
#' list object, each component of which being the peak table from a single
#' spectrum, after elution time correction.} \item{aligned.ftrs}{ Feature table
#' BEFORE weak signal recovery.} \item{final.ftrs}{Feature table after weak
#' signal recovery. This is the end product of the function.} \item{pk.times}{
#' Table of feature elution time BEFORE weak signal recovery.}
#' \item{final.times}{Table of feature elution time after weak signal
#' recovery.} \item{mz.tol}{The input mz.tol value by the user.}
#' \item{align.mz.tol}{The m/z tolerance level in the alignment across spectra,
#' either input from the user or automatically selected when the user input is
#' NA.} \item{align.chr.tol}{The retention time tolerance level in the
#' alignment across spectra, either input from the user or automatically
#' selected when the user input is NA.}
#' @author Tianwei Yu <tyu8@@sph.emory.edu>
#' @seealso proc.cdf, prof.to.feature, adjust.time, feature.align,
#' recover.weaker
#' @keywords models
#' @export
#'
cdf.to.ftr <- function(folder, file.pattern=".cdf", n.nodes=4, min.exp=2, min.pres=0.5, min.run=12, mz.tol=1e-5, baseline.correct.noise.percentile=0.05, shape.model="bi-Gaussian",  BIC.factor=2, baseline.correct=0, peak.estim.method="moment", min.bw=NA, max.bw=NA, sd.cut=c(0.01,500), sigma.ratio.lim=c(0.01, 100), component.eliminate=0.01, moment.power=1, subs=NULL, align.mz.tol=NA, align.chr.tol=NA, max.align.mz.diff=0.01, pre.process=FALSE, recover.mz.range=NA, recover.chr.range=NA, use.observed.range=TRUE,recover.min.count=3, intensity.weighted=FALSE)
{
    library(mzR)
    library(doParallel)
    setwd(folder)

    files<-dir(pattern=file.pattern, ignore.case = TRUE)
    files<-files[order(files)]
    if(!is.null(subs))
    {
        if(!is.na(subs[1])) files<-files[subs]
    }

    ###############################################################################################

    dir.create("error_files")
    message("***************************** prifiles --> feature lists *****************************")
    suf.prof<-paste(min.pres,min.run,mz.tol,baseline.correct,sep="_")
    suf<-paste(suf.prof, shape.model, sd.cut[1], sd.cut[2],component.eliminate, moment.power, sep="_")
    if(shape.model=="bi-Gaussian") suf<-paste(suf, sigma.ratio.lim[1], sigma.ratio.lim[2],sep="_")

    to.do<-paste(matrix(unlist(strsplit(tolower(files),"\\.")),nrow=2)[1,],suf, min.bw, max.bw,".feature",sep="_")
    to.do<-which(!(to.do %in% dir()))
    message(c("number of files to process: ", length(to.do)))

    if(length(to.do)>0)
    {
        grps<-round(seq(0, length(to.do), length=n.nodes+1))
        grps<-unique(grps)

        cl <- makeCluster(n.nodes,type='SOCK')
        registerDoParallel(cl)
        #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
        clusterEvalQ(cl, library(apLCMS))


        features<-foreach(i=2:length(grps)) %dopar%
        {
            this.subset<-to.do[(grps[i-1]+1):grps[i]]
            for(j in this.subset)
            {
                this.name<-paste(strsplit(tolower(files[j]),"\\.")[[1]][1],suf, min.bw, max.bw,".feature",sep="_")

                this.feature<-NA
                that.name<-paste(strsplit(tolower(files[j]),"\\.")[[1]][1],suf.prof,".profile",sep="_")

                processable<-"goodgood"
                processable<-try(this.prof<-proc.cdf(files[j], min.pres=min.pres, min.run=min.run, tol=mz.tol, baseline.correct=baseline.correct, baseline.correct.noise.percentile=baseline.correct.noise.percentile, do.plot=FALSE, intensity.weighted=intensity.weighted))
                if(substr(processable,1,5)=="Error")
                {
                    file.copy(from=files[j], to="error_files")
                    file.remove(files[j])
                }else{
                    save(this.prof,file=that.name)
                }

                if(substr(processable,1,5)!="Error")
                {
                    processable.2<-"goodgood"
                    processable.2<-try(this.feature<-prof.to.features(this.prof, min.bw=min.bw, max.bw=max.bw, sd.cut=sd.cut, shape.model=shape.model, estim.method=peak.estim.method, do.plot=FALSE, component.eliminate=component.eliminate, power=moment.power, BIC.factor=BIC.factor))

                    if(substr(processable.2,1,5)=="Error")
                    {
                        file.copy(from=files[j], to="error_files")
                        file.remove(files[j])
                        this.feature<-NA
                    }else{
                        save(this.feature, file=this.name)
                    }
                }
            }
            1
        }
        stopCluster(cl)

    }

    all.files<-dir()
    sel<-which(files %in% all.files)
    files<-files[sel]

    features<-new("list")
    for(i in 1:length(files))
    {
        this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf, min.bw, max.bw,".feature",sep="_")
        cat(this.name, " ")
        load(this.name)
        features[[i]]<-this.feature
    }

    gc()

    if(!pre.process)
    {
        ###############################################################################################
        message("****************************** time correction ***************************************")
        suf<-paste(suf,align.mz.tol,align.chr.tol,subs[1],subs[length(subs)],sep="_")
        this.name<-paste("time_correct_done_",suf,".bin",sep="")

        all.files<-dir()
        is.done<-all.files[which(all.files == this.name)]

        if(length(is.done)==0)
        {
            cl <- makeCluster(n.nodes,type='SOCK')
            registerDoParallel(cl)
            #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
            clusterEvalQ(cl, library(apLCMS))

            message(c("***** correcting time, CPU time (seconds) ",as.vector(system.time(f2<-adjust.time(features,mz.tol=align.mz.tol, chr.tol=align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)))[1]))
            save(f2,file=this.name)
            stopCluster(cl)

        }else{
            load(this.name)
        }
        gc()

        ###############################################################################################
        message("****************************  aligning features **************************************")
        suf<-paste(suf,min.exp,sep="_")
        this.name<-paste("aligned_done_",suf,".bin",sep="")
        all.files<-dir()
        is.done<-all.files[which(all.files == this.name)]
        if(length(is.done)==0)
        {
            cl <- makeCluster(n.nodes,type='SOCK')
            registerDoParallel(cl)
            #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
            clusterEvalQ(cl, library(apLCMS))

            message(c("***** aligning features, CPU time (seconds): ", as.vector(system.time(aligned<-feature.align(f2, min.exp=min.exp,mz.tol=align.mz.tol,chr.tol=align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)))[1]))
            save(aligned,file=this.name)
            stopCluster(cl)

        }else{
            load(this.name)
        }
        gc()

        ###############################################################################################
        message("**************************** recovering weaker signals *******************************")
        suf<-paste(suf,recover.mz.range, recover.chr.range, use.observed.range,sep="_")

        worklist<-paste(matrix(unlist(strsplit(tolower(files),"\\.")),nrow=2)[1,],suf,".recover",sep="_")
        to.do<-which(!(worklist %in% dir()))
        grps<-round(seq(0, length(to.do), length=n.nodes+1))
        grps<-unique(grps)

        message(c("number of files to process: ", length(to.do)))

        if(length(to.do)>0)
        {
            cl <- makeCluster(n.nodes,type='SOCK')
            registerDoParallel(cl)
            #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
            clusterEvalQ(cl, library(apLCMS))


            features.recov<-foreach(i=2:length(grps)) %dopar%
            {
                this.subset<-to.do[(grps[i-1]+1):grps[i]]
                for(j in this.subset)
                {
                    this.name<-paste(strsplit(tolower(files[j]),"\\.")[[1]][1],suf,".recover",sep="_")
                    this.recovered<-recover.weaker(filename=files[j], loc=j, aligned.ftrs=aligned$aligned.ftrs, pk.times=aligned$pk.times, align.mz.tol=aligned$mz.tol, align.chr.tol=aligned$chr.tol, this.f1=features[[j]], this.f2=f2[[j]], mz.range=recover.mz.range, chr.range=recover.chr.range, use.observed.range=use.observed.range, orig.tol=mz.tol, min.bw=min.bw, max.bw=max.bw, bandwidth=.5, recover.min.count=recover.min.count)
                    save(this.recovered, file=this.name)
                }
            }
            stopCluster(cl)
            gc()
        }

        new.aligned<-aligned
        for(i in 1:length(files))
        {
            this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf,".recover",sep="_")
            load(this.name)
            new.aligned$aligned.ftrs[,i+4]<-this.recovered$this.ftrs
            new.aligned$pk.times[,i+4]<-this.recovered$this.times
            new.aligned$features[[i]]<-this.recovered$this.f1
            new.aligned$f2[[i]]<-this.recovered$this.f2
            gc()
        }

        #################################################################################################
        rec<-new("list")
        colnames(aligned$aligned.ftrs)<-colnames(aligned$pk.times)<-colnames(new.aligned$aligned.ftrs)<-colnames(new.aligned$pk.times)<-c("mz","time","mz.min","mz.max",files)
        rec$features<-new.aligned$features
        rec$features2<-new.aligned$f2
        rec$aligned.ftrs<-aligned$aligned.ftrs
        rec$pk.times<-aligned$pk.times
        rec$final.ftrs<-new.aligned$aligned.ftrs
        rec$final.times<-new.aligned$pk.times
        rec$align.mz.tol<-new.aligned$mz.tol
        rec$align.chr.tol<-new.aligned$chr.tol
        rec$mz.tol<-mz.tol

        return(rec)
    }
}
