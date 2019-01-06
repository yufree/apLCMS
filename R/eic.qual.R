#' Internal function: Calculate the single predictor quality.
#'
#' For each column of an EIC data feature matrix, find its predictive power on
#' the m/z match to known metabolites.
#'
#'
#' @param eic.rec The EIC data feature matrix. Each row is an EIC. Each column
#' is a data feature.
#' @param known.mz The known m/z values to be matched to.
#' @param mass.matched A vector of indicators of whether the m/z of each EIC is
#' matched to the known m/z values. The default is NA, in which case it is
#' calculated within the function.
#' @param match.tol.ppm The tolerance level of m/z match.
#' @param do.plot Whether to produce plots of the ROCS.
#' @param pos.confidence The confidence level for the features matched to the
#' known feature list.
#' @param neg.confidence The confidence level for the features not matching to
#' the known feature list.
#' @return A matrix of four columns. The first two columns are the VUS and AUC
#' without uncertainty. The next two columns are the VUS and AUC with
#' uncertainty.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @references Bioinformatics. 30(20): 2941-2948.
#' @keywords models
eic.qual <-
function(eic.rec, known.mz, mass.matched=NA, match.tol.ppm=5, do.plot=FALSE, pos.confidence=0.99, neg.confidence=0.99)
{

    if(is.na(mass.matched[1]))
    {
        y<-mass.match(x=eic.rec[,2], known.mz=known.mz, match.tol.ppm=match.tol.ppm)
    }else{
        y<-mass.matched
    }

    rec.def<-rec.unc<-new("list")
    vus.fcauc<-matrix(nrow=ncol(eic.rec)-2, ncol=4)
    colnames(vus.fcauc)<-c("vus_1","fcauc_1","vus_uncertainty","fcauc_uncertainty")
    rownames(vus.fcauc)<-colnames(eic.rec)[-1:-2]

    for(i in 3:ncol(eic.rec))
    {
        x<-eic.rec[,i]
        r<-rocs.x(x[y==0], x[y==1], rep(1, sum(y==0)), rep(1, sum(y==1)), n.perm=1, FDR.cut=1, do.plot=FALSE)
        rec.def[[i-2]]<-r
        vus.fcauc[i-2, 1:2]<-unlist(r[1:2])

        r<-rocs.x(x[y==0], x[y==1], rep(0.99, sum(y==0)), rep(0.99, sum(y==1)), n.perm=1, FDR.cut=1, do.plot=FALSE)
        rec.unc[[i-2]]<-r
        vus.fcauc[i-2, 3:4]<-unlist(r[1:2])
    }

    if(do.plot)
    {
        par(mfrow=c(3,3))
        for(i in 1:min(9, length(rec.unc)))
        {
            plot(rec.def[[i]]$fp, rec.def[[i]]$tp, type="l", col="blue",xlab="FPR", ylab="TPR", main=rownames(vus.fcauc)[i])
            lines(rec.unc[[i]]$fp, rec.unc[[i]]$tp, col="red")
            abline(0,1,col="grey", lty=2)
        }
    }
    vus.fcauc
}
