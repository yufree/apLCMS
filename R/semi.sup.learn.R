#' Semi-supervised feature detection using machine learning approach.
#'
#' The semi-supervised procedure utilizes a database of known metabolites and
#' previously detected features to identify features in a new dataset. It is
#' recommended ONLY for experienced users. The user may need to construct the
#' known feature database that strictly follows the format described below.
#'
#' The function first conducts a machine-learning feature detection in the new
#' dataset. And the conducts the regular feature alignment, retention time
#' adjustment and weak signal recovery.
#'
#' @param folder The folder where all CDF files to be processed are located.
#' For example “C:/CDF/this_experiment”
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
#' @param min.exp If a feature is to be included in the final feature table, it
#' must be present in at least this number of spectra.
#' @param min.pres This is a parameter of thr run filter, to be passed to the
#' function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param min.run This is a parameter of thr run filter, to be passed to the
#' function proc.cdf(). Please see the help for proc.cdf() for details.
#' @param mz.tol The user can provide the m/z tolerance level for peak
#' identification. This value is expressed as the percentage of the m/z value.
#' This value, multiplied by the m/z value, becomes the cutoff level. Please
#' see the help for proc.cdf() for details.
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
#' @param subs If not all the CDF files in the folder are to be processed, the
#' user can define a subset using this parameter. For example, subs=15:30, or
#' subs=c(2,4,6,8)
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
#' @param use.learn whether to use machine learning approach. The default is
#' TRUE.
#' @param ridge.smoother.window The size of the smoother window used by the
#' kernel smoother to remove long ridge noise from each EIC.
#' @param smoother.window The smoother windows to use in data feature
#' generation.
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
#' intensities in the initial peak detection stage.
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
#' selected when the user input is NA.} \item{updated.known.table}{ The known
#' table updated using the newly processed data. It should be used for future
#' datasets generated using the same machine and LC column. }
#' \item{ftrs.known.table.pairing}{ The paring information between the feature
#' table of the current dataset and the known feature tabel. }
#' @author Tianwei Yu < tianwei.yu@@emory.edu>
#' @seealso cdf.to.ftrs, semi.sup, learn.cdf, prof.to.feature, adjust.time,
#' feature.align, recover.weaker
#' @keywords models
#' @export
semi.sup.learn <-
function(folder, file.pattern=".cdf", known.table=NA, n.nodes=4, min.exp=2, min.pres=0.3, min.run=4, mz.tol=1e-5, shape.model="bi-Gaussian",  baseline.correct=0, peak.estim.method="moment", min.bw=NA, max.bw=NA, sd.cut=c(0.01,500), component.eliminate=0.01, moment.power=1, sigma.ratio.lim=c(0.01, 100), subs=NULL, align.mz.tol=NA, align.chr.tol=NA, max.align.mz.diff=0.01, pre.process=FALSE, recover.mz.range=NA, recover.chr.range=NA, use.observed.range=TRUE, match.tol.ppm=5, new.feature.min.count=2, recover.min.count=3, use.learn=TRUE, ridge.smoother.window=50, smoother.window=c(1, 5, 10),pos.confidence=0.99, neg.confidence=0.99, max.ftrs.to.use=10, do.grp.reduce=TRUE, remove.bottom.ftrs=0, max.fpr=0.5, min.tpr=0.9, intensity.weighted=FALSE)
{


    library(mzR)
    library(doparallel)
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
    if(use.learn) suf.prof<-paste(min.run, min.pres, ridge.smoother.window, mz.tol, baseline.correct, max.ftrs.to.use, do.grp.reduce, remove.bottom.ftrs, max.fpr, min.tpr, sep="_")
    suf<-paste(suf.prof, shape.model, min.run, min.pres, sd.cut[1], sd.cut[2],component.eliminate, moment.power, sep="_")
    if(shape.model=="bi-Gaussian") suf<-paste(suf, sigma.ratio.lim[1], sigma.ratio.lim[2],sep="_")

    cl <- makeCluster(n.nodes,type='SOCK')
    registerDoParallel(cl)
    #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
    clusterEvalQ(cl, library(apLCMS))

    features<-foreach(i=1:length(files)) %dopar%
    {
        this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf, min.bw, max.bw,".feature",sep="_")
        all.files<-dir()
        do.exist<-all.files[which(all.files == this.name)]

        if(length(do.exist)==0)
        {
            this.feature<-NA
            that.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf.prof,".profile",sep="_")
            all.files<-dir()
            that.exist<-all.files[which(all.files == that.name)]

            processable<-"goodgood"
            if(length(that.exist)==0)
            {
                if(use.learn)
                {
                    message(c(i, "  ", that.name, "  ", as.vector(system.time(processable<-try(this.prof<-learn.cdf(files[i], known.mz=known.table[,6], tol=mz.tol, min.run=min.run, min.pres=min.pres, match.tol.ppm=match.tol.ppm, baseline.correct=baseline.correct, ridge.smoother.window=ridge.smoother.window, smoother.window=smoother.window,pos.confidence=pos.confidence, neg.confidence=neg.confidence, max.ftrs.to.use=max.ftrs.to.use, do.grp.reduce=do.grp.reduce, remove.bottom.ftrs=remove.bottom.ftrs, max.fpr=max.fpr, min.tpr=min.tpr, intensity.weighted=intensity.weighted))[1]))[1]))
                }else{
                    message(c(i, "  ", that.name, "  ", as.vector(system.time(processable<-try(this.prof<-proc.cdf(files[i], min.pres=min.pres, min.run=min.run, tol=mz.tol, baseline.correct=baseline.correct, do.plot=FALSE, intensity.weighted=intensity.weighted))[1]))[1]))
                }
                if(substr(processable[1],1,5)=="Error")
                {
                    file.copy(from=files[i], to="error_files")
                    file.remove(files[i])
                }else{
                    save(this.prof,file=that.name)
                }
            }else{
                message(paste(i, "  loading: ", that.name))
                load(that.name)
            }

            if(substr(processable[1],1,5)!="Error")
            {
                processable.2<-"goodgood"
                message(c(i, "  ", this.name,"  ",as.vector(system.time(processable.2<-try(this.feature<-prof.to.features(this.prof[[1]], min.bw=min.bw, max.bw=max.bw, sd.cut=sd.cut, shape.model=shape.model, estim.method=peak.estim.method, do.plot=FALSE, component.eliminate=component.eliminate, power=moment.power))))[1]))

                if(substr(processable.2,1,5)=="Error")
                {
                    file.copy(from=files[i], to="error_files")
                    file.remove(files[i])
                    this.feature<-NA
                }else{
                    save(this.feature, file=this.name)
                }
            }
        }else{
            message(paste(i, "  loading: ", this.name))
            load(this.name)
        }
        this.feature
    }
    stopCluster(cl)

    all.files<-dir()
    sel<-which(files %in% all.files)
    files<-files[sel]
    features<-features[sel]


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

        if(is.na(align.mz.tol)) align.mz.tol=2 * match.tol.ppm *1e-6
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

        if(is.na(align.mz.tol)) align.mz.tol=2 * match.tol.ppm *1e-6
        message(c("***** aligning features, CPU time (seconds): ", as.vector(system.time(aligned<-feature.align(f2, min.exp=min.exp,mz.tol=align.mz.tol,chr.tol=align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)))[1]))
        save(aligned,file=this.name)
        stopCluster(cl)
    }else{
        load(this.name)
    }
    gc()

    ###############################################################################################
    message("merging to known peak table")
    if(is.na(match.tol.ppm)) match.tol.ppm<-aligned$mz.tol*1e6

    mass.d2<-mass.match(aligned$aligned.ftrs[,1], known.table[,6],match.tol.ppm)
    mass.matched.pos<-which(mass.d2>0)

    known.assigned<-rep(0, nrow(known.table))
    new.assigned<-rep(0, nrow(aligned$aligned.ftrs))
    new.known.pairing<-c(0,0)

    for(i in mass.matched.pos)
    {
        if(new.assigned[i] == 0)
        {
            #find all potentially related known/newly found peaks
            old.sel.new<-i
            this.mz.thres<-aligned$aligned.ftrs[i,1]*match.tol.ppm/1e6
            sel.known<-which(abs(known.table[,6]-aligned$aligned.ftrs[i,1]) < this.mz.thres)
            sel.new<-NULL
            for(m in 1:length(sel.known)) sel.new<-c(sel.new, which(abs(aligned$aligned.ftrs[,1]-known.table[sel.known[m], 6]) < this.mz.thres))
            sel.known<-unique(sel.known)
            sel.new<-unique(sel.new)

            while(length(sel.new)>length(old.sel.new))
            {
                old.sel.new<-sel.new
                sel.known<-NULL
                for(m in 1:length(sel.new)) sel.known<-c(sel.known, which(abs(known.table[,6]-aligned$aligned.ftrs[sel.new[m],1]) < this.mz.thres))
                sel.new<-NULL
                for(m in 1:length(sel.known)) sel.new<-c(sel.new, which(abs(aligned$aligned.ftrs[,1]-known.table[sel.known[m], 6]) < this.mz.thres))
                sel.known<-unique(sel.known)
                sel.new<-unique(sel.new)
            }

            time.matched<-mass.matched<-matrix(0, ncol=length(sel.new), nrow=length(sel.known))

            for(k in 1:length(sel.known))
            {
                time.matched[k,]<-abs(aligned$aligned.ftrs[sel.new,2]-known.table[sel.known[k],11])
                mass.matched[k,]<-abs(known.table[sel.known[k],6]-aligned$aligned.ftrs[sel.new,1])
            }
            mass.matched<-1*(mass.matched/median(known.table[sel.known,6]) <= match.tol.ppm*1e-6)
            time.matched[mass.matched == 0] <- 1e10


            time.matched[is.na(time.matched)]<-aligned$chr.tol/2
            both.matched<-find.match(time.matched, unacceptable=aligned$chr.tol/2)

            for(m in 1:length(sel.new))
            {
                k<-which(both.matched[,m]==1)
                if(length(k)==1)
                {
                    if(known.assigned[sel.known[k]]==0)
                    {
                        new.assigned[sel.new[m]]<-1
                        known.assigned[sel.known[k]]<-1
                        new.known.pairing<-rbind(new.known.pairing, c(sel.new[m], sel.known[k]))
                    }
                }
            }
        }
    }
    colnames(new.known.pairing)<-c("new","known")
    new.known.pairing<-new.known.pairing[-1,]

    to.add.ftrs<-matrix(0, ncol=ncol(aligned$aligned.ftrs), nrow=nrow(known.table)-nrow(new.known.pairing))
    to.add.times<-matrix(NA, ncol=ncol(aligned$aligned.ftrs), nrow=nrow(known.table)-nrow(new.known.pairing))
    sel<-1:nrow(known.table)
    sel<-sel[-(new.known.pairing[,2])]

    to.add.ftrs[,1]<-to.add.times[,1]<-known.table[sel, 6]
    to.add.ftrs[,2]<-to.add.times[,2]<-known.table[sel, 11]
    to.add.ftrs[,3]<-to.add.times[,3]<-known.table[sel, 9]
    to.add.ftrs[,4]<-to.add.times[,4]<-known.table[sel, 10]

    aligned.ftrs<-rbind(aligned$aligned.ftrs, to.add.ftrs)
    pk.times<-rbind(aligned$pk.times, to.add.times)
    new.known.pairing<-rbind(new.known.pairing, cbind(1:nrow(to.add.ftrs)+nrow(aligned$aligned.ftrs), sel))
    gc()


    ###############################################################################################
    message("**************************** recovering weaker signals *******************************")
    suf<-paste(suf,recover.mz.range, recover.chr.range, use.observed.range,match.tol.ppm,new.feature.min.count,recover.min.count,sep="_")

    cl <- makeCluster(n.nodes,type='SOCK')
    registerDoParallel(cl)
    #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
    clusterEvalQ(cl, library(apLCMS))

    foreach(i=1:length(files)) %dopar%
    {
        this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf,"semi_sup.recover",sep="_")
        all.files<-dir()
        do.exist<-all.files[which(all.files == this.name)]

        if(length(do.exist)==0)
        {
            this.recovered<-recover.weaker(filename=files[i], loc=i, aligned.ftrs=aligned.ftrs, pk.times=pk.times, align.mz.tol=aligned$mz.tol, align.chr.tol=aligned$chr.tol, this.f1=features[[i]], this.f2=f2[[i]], mz.range=recover.mz.range, chr.range=recover.chr.range, use.observed.range=use.observed.range, orig.tol=mz.tol, min.bw=min.bw, max.bw=max.bw, bandwidth=.5, recover.min.count=recover.min.count)
            save(this.recovered, file=this.name)
        }
    }
    stopCluster(cl)

    new.aligned<-aligned
    new.aligned$orig.known.table<-known.table

    for(i in 1:length(files))
    {
        this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf,"semi_sup.recover",sep="_")
        load(this.name)
        aligned.ftrs[,i+4]<-this.recovered$this.ftrs
        pk.times[,i+4]<-this.recovered$this.times
        new.aligned$features[[i]]<-this.recovered$this.f1
        new.aligned$f2[[i]]<-this.recovered$this.f2
    }

    new.pk.times<-apply(pk.times[,-1:-4], 1, mean,na.rm=T)
    pk.times[is.na(pk.times[,2]),2]<-new.pk.times[is.na(pk.times[,2])]
    aligned.ftrs[,2]<-pk.times[,2]
    new.aligned$aligned.ftrs<-aligned.ftrs
    new.aligned$pk.times<-pk.times


    ################## updating known.table ############
    ### notice aligned.ftrs contains all features from the known table and new data

    known.2<-known.table
    num.exp.found<-apply(aligned.ftrs[,5:ncol(aligned.ftrs)]!=0, 1, sum)
    known.num.experiments<-unique(known.2[,7])
    new.num.experiments<-ncol(aligned.ftrs)-4+known.num.experiments

    for(i in 1:nrow(new.known.pairing))
    {
        if(num.exp.found[new.known.pairing[i,1]] >= 1)
        {
            known.2[new.known.pairing[i,2],]<-peak.characterize(existing.row=known.2[new.known.pairing[i,2],],ftrs.row=aligned.ftrs[new.known.pairing[i,1],], chr.row=pk.times[new.known.pairing[i,1],])
        }else{
            known.2[new.known.pairing[i,2],c(7,8)]<-c(new.num.experiments, known.num.experiments*known.2[new.known.pairing[i,2],8]/new.num.experiments)
        }
    }

    newly.found.ftrs<-which(!(1:nrow(aligned.ftrs) %in% new.known.pairing[,1]))
    for(i in newly.found.ftrs)
    {
        if(num.exp.found[i] >= new.feature.min.count)
        {
            this.row<-peak.characterize(existing.row=NA,ftrs.row=aligned.ftrs[i,], chr.row=pk.times[i,])
            #			this.row[c(7,8)]<-c(new.num.experiments, this.row[7]*this.row[8]/new.num.experiments)
            known.2<-rbind(known.2, this.row)
            new.known.pairing<-rbind(new.known.pairing, c(i,nrow(known.2)))
        }
    }

    #################################################################################################
    sel<-which(num.exp.found >= min.exp)
    new.aligned$aligned.ftrs<-new.aligned$aligned.ftrs[sel,]
    new.aligned$pk.times<-new.aligned$pk.times[sel,]
    missed.loc<-(1:nrow(aligned.ftrs))[!(1:nrow(aligned.ftrs) %in% new.known.pairing[,1])]
    missed.loc<-cbind(missed.loc, rep(NA, length(missed.loc)))
    new.known.pairing<-rbind(new.known.pairing, missed.loc)
    new.known.pairing<-new.known.pairing[order(new.known.pairing[,1]),]
    new.known.pairing<-new.known.pairing[sel,]
    new.known.pairing[,1]<-1:nrow(new.known.pairing)

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
    rec$updated.known.table<-known.2
    rec$ftrs.known.table.pairing<-new.known.pairing

    return(rec)

}
