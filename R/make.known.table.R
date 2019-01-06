#' Producing a table of known features based on a table of metabolites and a
#' table of allowable adducts.
#'
#' Given a table of known metabolites with original mass and charge
#' information, and a table of allowable adducts, this function outputs a new
#' table of potential features.
#'
#' For each allowable ion form, the function produces the m/z of every
#' metabolite given to it. The output table follows the format that is required
#' by the function semi.sup(), so that the user can directly use the table for
#' semi supervised feature detection.
#'
#' @param metabolite.table A table of known metabolites. See the description of
#' the object "metabolite.table" for details.
#' @param adduct.table A table of allowable adducts. See the description of the
#' object "adduct.table" for details.
#' @param ion.mode Character. Either "+" or "-".
#' @return A data frame containing the known metabolite ions. It contains 18
#' columns: "chemical_formula": the chemical formula if knonw; "HMDB_ID": HMDB
#' ID if known; "KEGG_compound_ID": KEGG compound ID if known; "neutral.mass":
#' the neutral mass if known: "ion.type": the ion form, such as H+, Na+, ...,
#' if known; "m.z": m/z value, either theoretical for known metabolites, or
#' mean observed value for unknown but previously found features;
#' "Number_profiles_processed": the total number of LC/MS profiles that were
#' used to build this database; "Percent_found": in what percentage was this
#' feature found historically amount all data processed in building this
#' database; "mz_min": the minimum m/z value observed for this feature;
#' "mz_max": the maximum m/z value observed for this feature; "RT_mean": the
#' mean retention time observed for this feature; "RT_sd": the standard
#' deviation of retention time observed for this feature; "RT_min": the minimum
#' retention time observed for this feature; "RT_max": the maximum retention
#' time observed for this feature; "int_mean.log.": the mean log intensity
#' observed for this feature; "int_sd.log.": the standard deviation of log
#' intensity observed for this feature; "int_min.log.": the minimum log
#' intensity observed for this feature; "int_max.log.": the maximum log
#' intensity observed for this feature;
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @seealso metabolite.table, adduct.table, semi.sup
#' @references Yu T, Park Y, Li S, Jones DP (2013) Hybrid feature detection and
#' information accumulation using high-resolution LC-MS metabolomics data. J.
#' Proteome Res. 12(3):1419-27.
#' @keywords models
#' @examples
#'
#'
#'
#' data(metabolite.table)
#' data(adduct.table)
#' known.table.example<-make.known.table(metabolite.table[1001:1020,], adduct.table[1:4,])
#'
#'
#' @export
make.known.table <-
function(metabolite.table, adduct.table, ion.mode="+")
###		the metabolite.table should be a four-column data frame. The columns should be chemical formula, HMDB_ID, KEGG ID, and monoisotopic mass
###		The two IDs can be substitute by others
###		the adduct table should be a four-column data frame. The columns are adduct type, mass divider (absolute value of the number of charges), mass addition (after dividing), and charge of the ion form.
{

    metabolite.table<-metabolite.table[order(metabolite.table[,4]),]
    for(i in 1:3) metabolite.table[,i]<-as.vector(metabolite.table[,i])

    metabolite.table<-cbind(metabolite.table, matrix(NA, nrow=nrow(metabolite.table),ncol=14))
    colnames(metabolite.table)[5:ncol(metabolite.table)]<-c("ion.type","m.z","Number_profiles_processed","Percent_found","mz_min","mz_max","RT_mean","RT_sd","RT_min","RT_max","int_mean(log)","int_sd(log)","int_min(log)","int_max(log)")

    ###		step 1:take out charged ions

    l<-nchar(metabolite.table[,1])
    last.char<-substr(metabolite.table[,1], l, l)
    sel<-which(last.char=="+" | last.char=="-")

    rm.first.row<-TRUE
    new.table<-matrix(NA, nrow=1, ncol=ncol(metabolite.table))
    colnames(new.table)<-colnames(metabolite.table)

    if(length(sel)>0)
    {
        new.table<-rbind(new.table, metabolite.table[sel,])

        l<-nchar(new.table[,1])-1
        num.charge<-substr(new.table[,1], l, l)
        num.charge<-as.numeric(num.charge)

        new.table[,5]<-"orig"
        new.table[,6]<-new.table[,4]
        new.table[which(!is.na(num.charge)),6]<-new.table[which(!is.na(num.charge)),6]/num.charge[which(!is.na(num.charge))]

        metabolite.table<-metabolite.table[-sel,]
    }


    ###		step 2: merge metaoblites with same m/z values into a single row

    n<-1
    m<-2
    to.remove<-rep(0, nrow(metabolite.table))

    while(m <= nrow(metabolite.table))
    {
        if(abs(metabolite.table[m,4]-metabolite.table[n,4])<1e-10)
        {
            if(metabolite.table[n,1] != metabolite.table[m,1]) metabolite.table[n,1]<-paste(metabolite.table[n,1],metabolite.table[m,1], sep="/")
            for(j in 2:3) metabolite.table[n,j]<-paste(metabolite.table[n,j],metabolite.table[m,j], sep="/")
            to.remove[m]<-1
            m<-m+1
        }else{
            n<-m
            m<-m+1
        }
        cat("(", n, m, ")")
    }

    if(sum(to.remove)>0) metabolite.table<-metabolite.table[-which(to.remove==1),]


    ###		generate ion adducts for uncharged metabolites, merge results to the new table

    neutral.table<-metabolite.table

    for(n in 1:nrow(adduct.table))
    {
        this<-neutral.table
        this[,5]<-adduct.table[n,1]
        this[,6]<-neutral.table[,4]/adduct.table[n,2]+adduct.table[n,3]

        new.table<-rbind(new.table, this)
    }
    if(rm.first.row) new.table<-new.table[-1,]

    ###		ion mode

    to.remove<-rep(0, nrow(new.table))
    l<-nchar(new.table[,1])
    last.char<-substr(new.table[,1], l, l)
    if(ion.mode == "+") to.remove[which(last.char == "-")]<-1
    if(ion.mode == "-") to.remove[which(last.char == "+")]<-1
    if(sum(to.remove)>0) new.table<-new.table[-which(to.remove==1),]

    ###		are there m/z overlaps?

    new.table<-new.table[order(new.table[,6]),]

    n<-1
    m<-2
    to.remove<-rep(0, nrow(new.table))

    while(m <= nrow(new.table))
    {
        if(abs(new.table[m,6]-new.table[n,6])<1e-10)
        {
            if(new.table[n,1] != new.table[m,1]) new.table[n,1]<-paste(new.table[n,1],new.table[m,1], sep="/")
            for(j in 2:5) new.table[n,j]<-paste(new.table[n,j],new.table[m,j], sep="/")
            to.remove[m]<-1
            m<-m+1
        }else{
            n<-m
            m<-m+1
        }
        cat("*(", n, m, ")")
    }

    if(sum(to.remove)>0) new.table<-new.table[-which(to.remove==1),]

    ###		output the table
    return(new.table)
}
