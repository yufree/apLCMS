#' Adjust retention time across spectra.
#'
#' This function adjusts the retention time in each LC/MS profile to achieve
#' better between-profile agreement.
#'
#' The function first searches for the m/z tolerance level using a mixture
#' model. After the mz.tol is obtained, the peaks are grouped based on it. The
#' function then searches for the retention time tolerance level. Because the
#' peaks are grouped using m/z, only metabolites that share m/z require this
#' parameter. A rather lenient retention time tolerance level is found using a
#' mixture model.
#'
#' The profile with the highest number of peaks is selected as the template and
#' every other spetrum is adjusted to it one at a time. At every m/z value, if
#' each of the two spetra has just one peak, and the peaks are within the
#' retention time tolerance range, the pair of retention time values are used
#' in the curve fitting. A kernel smoother is fitted using the difference in
#' retention time against the retention time in the profile to be adjusted.
#'
#' @param features A list object. Each component is a matrix which is the
#' output from proc.to.feature().
#' @param mz.tol The m/z tolerance level for peak alignment. The default is NA,
#' which allows the program to search for the tolerance level based on the
#' data. This value is expressed as the percentage of the m/z value. This
#' value, multiplied by the m/z value, becomes the cutoff level.
#' @param chr.tol The retention time tolerance level for peak alignment. The
#' default is NA, which allows the program to search for the tolerance level
#' based on the data.
#' @param colors The vector of colors to be used for the line plots of time
#' adjustments. The default is NA, in which case the program uses a set of
#' default color set.
#' @param find.tol.max.d Argument passed to find.tol(). Consider only m/z diffs
#' smaller than this value. This is only used when the mz.tol is NA.
#' @param max.align.mz.diff As the m/z tolerance is expressed in relative terms
#' (ppm), it may not be suitable when the m/z range is wide. This parameter
#' limits the tolerance in absolute terms. It mostly influences feature
#' matching in higher m/z range.
#' @return A list object with the exact same structure as the input object
#' features, i.e. one matrix per profile being processed. The only difference
#' this output object has with the input object is that the retention time
#' column in each of the matrices is changed to new adjusted values.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @seealso feature.align
#' @keywords models
#' @examples
#' data(features)
#' adjusted<-adjust.time(features, colors=c("red","blue","green","cyan"))
#' @export
adjust.time <-
function(features,mz.tol=NA, chr.tol=NA,colors=NA,find.tol.max.d=1e-4, max.align.mz.diff=0.01)  # features is a list project, each sub-object is a matrix as identified by prof.to.features
{

    num.exp<-nrow(summary(features))
    if(num.exp>1)
    {
        par(mfrow=c(2,2))
        plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",main="",axes=FALSE)
        text(x=0,y=0,"Retention time \n adjustment",cex=2)

        a<-summary(features)
        sizes<-as.numeric(a[,1])/ncol(features[[1]])
        sizes<-cumsum(sizes)
        #		sel<-max(which(sizes<=5e6))
        sel<-length(sizes)

        mz<-chr<-lab<-rep(0, sizes[sel])
        sizes<-c(0, sizes)

        for(i in 1:sel)
        {
            mz[(sizes[i]+1):sizes[i+1]]<-features[[i]][,1]
            chr[(sizes[i]+1):sizes[i+1]]<-features[[i]][,2]
            lab[(sizes[i]+1):sizes[i+1]]<-i
        }

        if(is.na(mz.tol))
        {
            mz.tol<-find.tol(mz,uppermost=find.tol.max.d)
        }else{
            plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",main="m/z tolerance level given",axes=FALSE)
            text(x=0,y=0,mz.tol,cex=1.2)
        }

        if(!is.na(chr.tol))
        {
            plot(c(-1,1),c(-1,1),type="n",xlab="",ylab="",main="retention time \n tolerance level given",axes=FALSE)
            text(x=0,y=0,chr.tol,cex=1.2)
        }

        all.ft<-find.tol.time(mz, chr, lab, num.exp=num.exp, mz.tol=mz.tol, chr.tol=chr.tol, max.mz.diff=max.align.mz.diff)
        chr.tol<-all.ft$chr.tol

        message("**** performing time correction ****")
        message(paste("m/z tolerance level: ", mz.tol))
        message(paste("time tolerance level:", chr.tol))

        for(i in 1:num.exp)
        {
            this<-features[[i]]
            sel<-which(all.ft$lab == i)
            that<-cbind(all.ft$mz[sel], all.ft$chr[sel], all.ft$grps[sel])
            this<-this[order(this[,1], this[,2]), ]
            that<-that[order(that[,1], that[,2]), ]
            this<-cbind(this, rep(i, nrow(this)), that[,3])
            features[[i]]<-this
        }

        num.ftrs<-as.vector(table(all.ft$lab))
        template<-which(num.ftrs==max(num.ftrs))[1]
        message(paste("the template is sample",template))
        if(is.na(colors[1])) colors<-c("red","blue","dark blue","orange","green","yellow","cyan","pink","violet","bisque","azure","brown","chocolate",rep("grey",num.exp))

        candi<-features[[template]][,1:2]

        features.2<-foreach(j=1:num.exp) %dopar%
        {
            this.feature<-features[[j]]
            if(j != template)
            {
                this.comb<-rbind(cbind(candi, rep(template,nrow(candi))),cbind(this.feature[,1:2],rep(j,nrow(this.feature))))
                this.comb<-this.comb[order(this.comb[,1]),]
                l<-nrow(this.comb)

                sel<-which(this.comb[2:l,1]-this.comb[1:(l-1),1]< mz.tol * this.comb[1:(l-1),1] * 2 &
                abs(this.comb[2:l,2]-this.comb[1:(l-1),2])< chr.tol &
                this.comb[2:l,3] != this.comb[1:(l-1),3]
                )
                if(length(sel)<20)
                {
                    cat("too few, aborted")
                }else{
                    all.ftr.table <- cbind(this.comb[sel,2], this.comb[sel+1,2])
                    to.flip <- which(this.comb[sel,3] == j)
                    temp<-all.ftr.table[to.flip,2]
                    all.ftr.table[to.flip,2]<-all.ftr.table[to.flip,1]
                    all.ftr.table[to.flip,1]<-temp

                    # now the first column is the template retention time. the second column is the to-be-adjusted retention time

                    cat(c("sample",j,"using",nrow(all.ftr.table),","))
                    if(j %% 3 ==0) cat("\n")

                    all.ftr.table<-all.ftr.table[order(all.ftr.table[,2]),]

                    this.dev<-all.ftr.table[,2]             ## the to be adjusted time
                    aver.time<-all.ftr.table[,1]-this.dev   ## the difference between the true time and the to-be-adjusted time

                    this.feature<-this.feature[order(this.feature[,2], this.feature[,1]),]
                    this.corrected<-this.old<-this.feature[,2]
                    to.correct<-this.old[this.old>=min(this.dev) & this.old<=max(this.dev)]

                    this.smooth<-ksmooth(this.dev, aver.time, kernel="normal", bandwidth=(max(this.dev)-min(this.dev))/5, x.points=to.correct)

                    this.corrected[this.old>=min(this.dev) & this.old<=max(this.dev)] <- this.smooth$y + to.correct
                    this.corrected[this.old < min(this.dev)]<-  this.corrected[this.old < min(this.dev)] + mean(this.smooth$y[this.smooth$x==min(this.smooth$x)])
                    this.corrected[this.old > max(this.dev)]<-  this.corrected[this.old > max(this.dev)] + mean(this.smooth$y[this.smooth$x==max(this.smooth$x)])
                    this.feature[,2]<-this.corrected
                    this.feature<-this.feature[order(this.feature[,1], this.feature[,2]),]
                }
            }
            this.feature
        }
    }else{
        message("Only one sample.  No need to correct for time.")
    }

    plot(range(features[[1]][,2]), c(-chr.tol, chr.tol),type="n", ylab="Retention time deviation", xlab="Original Retention time")
    for(i in 1:num.exp)
    {
        features[[i]]<-features[[i]][order(features[[i]][,1], features[[i]][,2]),]
        points(features[[i]][,2], features.2[[i]][,2]-features[[i]][,2],col=colors[i], cex=.2)
    }
    return(features.2)
}
