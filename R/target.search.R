#' Targeted search of metabolites with given m/z and (optional) retention time
#'
#' The function conducts targeted search only. The search is based on m/z and
#' (optionally) retention time. If there are sufficient number of peaks (>=100)
#' in each profile, the function will conduct retention time correction and
#' peak alignment, in order to reduce potential redundancies.
#'
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
#' @param min.bw The minimum bandwidth in the smoother in prof.to.features().
#' Please see the help file for prof.to.features() for details.
#' @param max.bw The maximum bandwidth in the smoother in prof.to.features().
#' Please see the help file for prof.to.features() for details.
#' @param subs If not all the CDF files in the folder are to be processed, the
#' user can define a subset using this parameter. For example, subs=15:30, or
#' subs=c(2,4,6,8)
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
#' @return \item{features}{ A list object, each component of which being the
#' peak table from a single spectrum.} \item{filled.ftrs}{The target features
#' are filled one by one. Notice this table may contain duplicates if some
#' target features are too close. } \item{reduced.ftrs}{If the number of target
#' features are big enough (>=100 detected in each profile), retention time
#' correction and peak alignments are conducted to generate this feature table
#' without redundancy.} \item{filled.times}{The target features are filled one
#' by one. This is the retention time table. Notice this table may contain
#' duplicates if some target features are too close. } \item{reduced.times}{If
#' the number of target features are big enough (>=100 detected in each
#' profile), retention time correction and peak alignments are conducted to
#' generate this feature table without redundancy. This is the retention time
#' table of the aligned features. }
#' @author Tianwei Yu < tianwei.yu@@emory.edu>
#' @seealso cdf.to.ftrs, proc.cdf, prof.to.feature, adjust.time, feature.align,
#' recover.weaker
#' @keywords models
#' @export

target.search <- function(folder, file.pattern=".cdf", known.table=NA, n.nodes=4, min.exp=2, min.bw=NA, max.bw=NA, subs=NULL, align.mz.tol=2e-5, align.chr.tol=150, max.align.mz.diff=0.01, recover.mz.range=NA, recover.chr.range=NA, use.observed.range=TRUE, match.tol.ppm=5, new.feature.min.count=2, recover.min.count=3)
{
    target.onefile<-function(filename, loc, aligned.ftrs, pk.times, align.mz.tol, align.chr.tol, mz.range=NA, chr.range=NA, use.observed.range=TRUE, orig.tol=1e-5,min.bw=NA,max.bw=NA,bandwidth=.5, recover.min.count=3)
    {
        duplicate.row.remove<-function(new.table)
        {
            new.table<-new.table[order(new.table[,1], new.table[,2], new.table[,5]),]
            n<-1
            m<-2
            to.remove<-rep(0, nrow(new.table))

            while(m <= nrow(new.table))
            {
                if(abs(new.table[m,1]-new.table[n,1])<1e-10 & abs(new.table[m,2]-new.table[n,2])<1e-10 & abs(new.table[m,5]-new.table[n,5])<1e-10)
                {
                    to.remove[m]<-1
                    m<-m+1
                }else{
                    n<-m
                    m<-m+1
                }
                #cat("*(", n, m, ")")
            }

            if(sum(to.remove)>0) new.table<-new.table[-which(to.remove==1),]
            new.table
        }

        library(splines)
        library(mzR)
        if(is.na(mz.range)) mz.range<-1.5*align.mz.tol
        if(is.na(chr.range)) chr.range<-align.chr.tol/2

        this.raw<-load.lcms(filename)
        masses<-this.raw$masses
        intensi<-this.raw$intensi
        labels<-this.raw$labels
        times<-this.raw$times
        rm(this.raw)

        times<-times[order(times)]
        base.curve<-unique(times)
        base.curve<-base.curve[order(base.curve)]
        base.curve<-cbind(base.curve, base.curve*0)


        masses<-c(masses, -100000)
        mass.breaks<-which(masses[1:(length(masses)-1)] > masses[2:length(masses)])
        mass.breaks<-c(0,mass.breaks)
        masses<-masses[1:(length(masses)-1)]

        curr.order<-order(masses)
        intensi<-intensi[curr.order]
        labels<-labels[curr.order]
        masses<-masses[curr.order]


        if(is.na(min.bw)) min.bw<-diff(range(times, na.rm=TRUE))/60
        if(is.na(max.bw)) max.bw<-diff(range(times, na.rm=TRUE))/15
        if(min.bw >= max.bw) min.bw<-max.bw/4

        base.curve<-times
        aver.diff<-mean(diff(base.curve))
        base.curve<-cbind(base.curve, base.curve*0)
        all.times<-base.curve[,1]
        if(all.times[1]>0) all.times<-c(0,all.times)
        all.times<-c(all.times, 2*all.times[length(all.times)]-all.times[length(all.times)-1])
        all.times<-(all.times[1:(length(all.times)-1)]+all.times[2:length(all.times)])/2
        all.times<-all.times[2:length(all.times)]-all.times[1:(length(all.times)-1)]

        this.ftrs<-aligned.ftrs[, (loc+4)]
        this.times<-pk.times[,(loc+4)]
        custom.mz.tol<-mz.range*aligned.ftrs[,1]
        observed.mz.range<-(aligned.ftrs[,4]-aligned.ftrs[,3])/2
        #	if(use.observed.range) custom.mz.tol[which(custom.mz.tol < observed.mz.range)]<-observed.mz.range[which(custom.mz.tol < observed.mz.range)]

        custom.chr.tol<-rep(chr.range, nrow(aligned.ftrs))

        if(use.observed.range)
        {
            observed.chr.range<-(apply(pk.times[,5:ncol(pk.times)],1,max)-apply(pk.times[,5:ncol(pk.times)],1,min))/2
            num.present<-apply(!is.na(pk.times[,5:ncol(pk.times)]),1,sum)
            custom.chr.tol[which(num.present>=5 & custom.chr.tol > observed.chr.range)]<-observed.chr.range[which(num.present>=5 & custom.chr.tol > observed.chr.range)]
        }

        l<-length(masses)
        curr.bw<-0.5*orig.tol*max(masses)
        all.mass.den<-density(masses, weights=intensi/sum(intensi), bw=curr.bw, n=2^min(15, floor(log2(l))-2))
        all.mass.turns<-find.turn.point(all.mass.den$y)
        all.mass.vlys<-all.mass.den$x[all.mass.turns$vlys]
        breaks<-c(0, unique(round(approx(masses,1:l,xout=all.mass.vlys,rule=2,ties='ordered')$y))[-1])
        this.mz<-rep(NA, length(this.ftrs))
        this.f1<-matrix(rep(NA,5),nrow=1)
        target.time<-aligned.ftrs[,2]

        for(i in 1:length(this.ftrs))
        {
            if(this.ftrs[i] == 0 & aligned.ftrs[i,1] < masses[breaks[length(breaks)]])
            {
                if(aligned.ftrs[i,1] <= masses[breaks[2]])
                {
                    this.found<-c(1,2)
                }else{
                    this.found<-c(which(abs(masses[breaks]-aligned.ftrs[i,1]) < custom.mz.tol[i]), min(which(masses[breaks] > aligned.ftrs[i,1])), max(which(masses[breaks] < aligned.ftrs[i,1])))+1
                    this.found<-c(min(this.found), max(this.found))
                }

                if(length(this.found)>1)
                {
                    this.sel<-(breaks[this.found[1]]+1):breaks[this.found[2]]
                    this.masses<-masses[this.sel]
                    this.labels<-labels[this.sel]
                    this.intensi<-intensi[this.sel]

                    this.bw=0.5*orig.tol*aligned.ftrs[i,1]
                    mass.den<-density(this.masses, weights=this.intensi/sum(this.intensi), bw=this.bw)
                    mass.den$y[mass.den$y < min(this.intensi)/10]<-0
                    mass.turns<-find.turn.point(mass.den$y)
                    mass.pks<-mass.den$x[mass.turns$pks]
                    mass.vlys<-c(-Inf, mass.den$x[mass.turns$vlys], Inf)
                    mass.pks<-mass.pks[which(abs(mass.pks-aligned.ftrs[i,1]) < custom.mz.tol[i]/1.5)]

                    if(length(mass.pks) > 0)
                    {
                        this.rec<-matrix(c(Inf, Inf, Inf),nrow=1)
                        for(k in 1:length(mass.pks))
                        {
                            mass.lower<-max(mass.vlys[mass.vlys < mass.pks[k]])
                            mass.upper<-min(mass.vlys[mass.vlys > mass.pks[k]])

                            that.sel<-which(this.masses > mass.lower & this.masses <= mass.upper)
                            if(length(that.sel) > recover.min.count)
                            {
                                that.labels<-this.labels[that.sel]
                                that.masses<-this.masses[that.sel]
                                that.intensi<-this.intensi[that.sel]
                                that.order<-order(that.labels)
                                that.labels<-that.labels[that.order]
                                that.masses<-that.masses[that.order]
                                that.intensi<-that.intensi[that.order]

                                that.prof<-merge.seq.3(that.labels, that.masses, that.intensi)

                                that.mass<-sum(that.prof[,1]*that.prof[,3])/sum(that.prof[,3])
                                curr.rec<-c(that.mass, NA,NA)
                                if(nrow(that.prof) < 10)
                                {

                                    if(!is.na(target.time[i]))
                                    {
                                        thee.sel<-which(abs(that.prof[,2]-target.time[i]) < custom.chr.tol[i]*2)
                                    }else{
                                        thee.sel<-1:nrow(that.prof)
                                    }
                                    if(length(thee.sel)>recover.min.count)
                                    {
                                        if(length(thee.sel)>1)
                                        {
                                            that.inte<-interpol.area(that.prof[thee.sel,2], that.prof[thee.sel,3], base.curve[,1], all.times)
                                        }else{
                                            that.inte<-that.prof[thee.sel,3]*aver.diff
                                        }
                                        curr.rec[3]<-that.inte
                                        curr.rec[2]<-median(that.prof[thee.sel,2])
                                        this.rec<-rbind(this.rec, curr.rec)
                                    }
                                }else{
                                    this<-that.prof[,2:3]
                                    this<-this[order(this[,1]),]
                                    this.span<-range(this[,1])
                                    this.curve<-base.curve[base.curve[,1]>=this.span[1] & base.curve[,1] <=this.span[2],]
                                    this.curve[this.curve[,1] %in% this[,1],2]<-this[,2]

                                    bw<-min(max(bandwidth*(max(this[,1])-min(this[,1])),min.bw), max.bw)
                                    this.smooth<-ksmooth(this.curve[,1],this.curve[,2], kernel="normal",bandwidth=bw)
                                    smooth.y<-this.smooth$y
                                    turns<-find.turn.point(smooth.y)
                                    pks<-this.smooth$x[turns$pks]
                                    vlys<-this.smooth$x[turns$vlys]
                                    vlys<-c(-Inf, vlys, Inf)

                                    pks.n<-pks
                                    for(m in 1:length(pks))
                                    {
                                        this.vlys<-c(max(vlys[which(vlys<pks[m])]), min(vlys[which(vlys>pks[m])]))
                                        pks.n[m]<-sum(this[,1]>= this.vlys[1] & this[,1] <= this.vlys[2])
                                    }

                                    if(!is.na(target.time[i]))
                                    {
                                        pks.d<-abs(pks-target.time[i])    # distance from the target peak location
                                        pks.d[pks.n==0]<-Inf
                                        pks<-pks[which(pks.d==min(pks.d))[1]]
                                    }else{
                                        pks<-pks[pks.n>recover.min.count]
                                    }

                                    all.pks<-pks
                                    all.vlys<-vlys
                                    all.this<-this

                                    if(length(all.pks)>0)
                                    {
                                        for(pks.i in 1:length(all.pks))
                                        {
                                            pks<-all.pks[pks.i]
                                            vlys<-c(max(all.vlys[which(all.vlys<pks)]), min(all.vlys[which(all.vlys>pks)]))

                                            this<-all.this[which(all.this[,1] >= vlys[1] & all.this[,1] <= vlys[2]),]
                                            if(is.null(nrow(this)))
                                            {
                                                curr.rec[3]<-this[2]*aver.diff
                                                curr.rec[2]<-this[1]
                                            }else{
                                                x<-this[,1]
                                                y<-this[,2]

                                                if(nrow(this)>=10)
                                                {
                                                    miu<-sum(x*y)/sum(y)
                                                    sigma<-sqrt(sum(y*(x-miu)^2)/sum(y))
                                                    if(sigma==0)
                                                    {
                                                        curr.rec[3]<-sum(y)*aver.diff
                                                        curr.rec[2]<-miu
                                                    }else{
                                                        fitted<-dnorm(x, mean=miu, sd=sigma)
                                                        this.sel<-y>0 & fitted/dnorm(miu,mean=miu,sd=sigma) >1e-2
                                                        sc<-exp(sum(fitted[this.sel]^2*log(y[this.sel]/fitted[this.sel])/sum(fitted[this.sel]^2)))
                                                    }
                                                }else{
                                                    sc<-interpol.area(x, y, base.curve[,1], all.times)
                                                    miu<-median(x)
                                                }
                                                curr.rec[3]<-sc
                                                curr.rec[2]<-miu
                                            }
                                            this.rec<-rbind(this.rec, curr.rec)
                                        }
                                    }
                                }
                            }
                        }

                        if(!is.na(target.time[i]))
                        {
                            this.sel<-which(abs(this.rec[,2]-target.time[i])<custom.chr.tol[i])
                        }else{
                            this.sel<-1:nrow(this.rec)
                            this.sel<-this.sel[this.sel != 1]
                        }


                        if(length(this.sel)>0)
                        {
                            if(length(this.sel)>1)
                            {
                                if(!is.na(target.time[i]))
                                {
                                    this.d<-(this.rec[,2]-target.time[i])^2/custom.chr.tol[i]^2+(this.rec[,1]-aligned.ftrs[i,1])^2/custom.mz.tol[i]^2
                                    this.sel<-which(this.d==min(this.d))[1]
                                }else{
                                    this.d<-abs(this.rec[,1]-aligned.ftrs[i,1])
                                    this.sel<-which(this.d==min(this.d))[1]
                                }
                            }
                            this.f1<-rbind(this.f1,c(this.rec[this.sel,1],this.rec[this.sel,2],NA, NA, this.rec[this.sel,3]))
                            this.ftrs[i]<-this.rec[this.sel, 3]
                            this.times[i]<-this.rec[this.sel,2]
                            this.mz[i]<-this.rec[this.sel,1]
                        }
                    }
                }
            }
        }
        this.f1<-this.f1[-1,]
        to.return<-new("list")
        to.return$this.mz<-this.mz
        to.return$this.ftrs<-this.ftrs
        to.return$this.times<-this.times
        if(!is.null(nrow(this.f1))) this.f1<-duplicate.row.remove(this.f1)
        to.return$this.f1<-this.f1

        return(to.return)
    }

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
    message("Making target peak table")

    aligned.ftrs<-matrix(0, ncol=4+length(files), nrow=nrow(known.table))
    pk.times<-matrix(NA, ncol=4+length(files), nrow=nrow(known.table))

    aligned.ftrs[,1]<-pk.times[,1]<-known.table[, 6]
    aligned.ftrs[,2]<-pk.times[,2]<-known.table[, 11]
    aligned.ftrs[,3]<-pk.times[,3]<-known.table[, 9]
    aligned.ftrs[,4]<-pk.times[,4]<-known.table[, 10]

    ###############################################################################################
    message("**************************** recovering target signals *******************************")
    suf<-paste("recover", recover.mz.range, recover.chr.range, use.observed.range,match.tol.ppm,new.feature.min.count,recover.min.count,nrow(aligned.ftrs), sep="_")


    cl <- makeCluster(n.nodes,type='SOCK')
    registerDoParallel(cl)
    #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
    clusterEvalQ(cl, library(apLCMS))

    foreach(i=1:length(files)) %dopar%
    {
        this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf,"targeted.recover",sep="_")
        all.files<-dir()
        do.exist<-all.files[which(all.files == this.name)]

        if(length(do.exist)==0)
        {
            this.recovered<-target.onefile(filename=files[i], loc=i, aligned.ftrs=aligned.ftrs, pk.times=pk.times, align.mz.tol=align.mz.tol, align.chr.tol=align.chr.tol, mz.range=recover.mz.range, chr.range=recover.chr.range, use.observed.range=use.observed.range, orig.tol=align.mz.tol, min.bw=min.bw, max.bw=max.bw, bandwidth=.5, recover.min.count=recover.min.count)
            save(this.recovered, file=this.name)
            gc()
        }
    }
    stopCluster(cl)

    ##############################################################################################
    message("loading feature tables after target search")
    features.recov<-new("list")

    for(i in 1:length(files))
    {
        this.name<-paste(strsplit(tolower(files[i]),"\\.")[[1]][1],suf,"targeted.recover",sep="_")
        load(this.name)
        aligned.ftrs[,i+4]<-this.recovered$this.ftrs
        pk.times[,i+4]<-this.recovered$this.times
        features.recov[[i]]<-this.recovered$this.f1
        gc()
    }

    ##############################################################################################

    if(min(unlist(lapply(features.recov, nrow))) >= 100)
    {
        message("****************************** time correction *************************")
        suf<-paste(suf,"round 2",sep="_")
        this.name<-paste("time_correct_done_",suf,".bin",sep="")

        all.files<-dir()
        is.done<-all.files[which(all.files == this.name)]

        if(length(is.done)==0)
        {
            cl <- makeCluster(n.nodes,type='SOCK')
            registerDoParallel(cl)
            #clusterEvalQ(cl, source("~/Desktop/Dropbox/1-work/apLCMS_code/new_proc_cdf.r"))
            clusterEvalQ(cl, library(apLCMS))

            message(c("***** correcting time, CPU time (seconds) ",as.vector(system.time(f2.recov<-adjust.time(features.recov, mz.tol=align.mz.tol, chr.tol=align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)))[1]))
            save(f2.recov,file=this.name)
            stopCluster(cl)

        }else{
            load(this.name)
        }
        gc()
        ###############################################################################################
        message("**************************** aligning features *************************")
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

            message(c("***** aligning features, CPU time (seconds): ", as.vector(system.time(aligned.recov<-feature.align(f2.recov, min.exp=min.exp,mz.tol=align.mz.tol,chr.tol=align.chr.tol, find.tol.max.d=10*mz.tol, max.align.mz.diff=max.align.mz.diff)))[1]))
            save(aligned.recov,file=this.name)
            stopCluster(cl)
        }else{
            load(this.name)
        }
        gc()
    }

    #################################################################################################

    rec<-new("list")
    colnames(aligned.ftrs)<-colnames(pk.times)<-c("mz","time","mz.min","mz.max",files)
    rec$features<-features.recov
    rec$filled.ftrs<-aligned.ftrs
    rec$filled.times<-pk.times
    if(min(unlist(lapply(features.recov, nrow))) >= 100)
    {
        rec$reduced.ftrs<-aligned.recov$aligned.ftrs
        rec$reduced.times<-aligned.recov$pk.times
    }

    return(rec)
}
