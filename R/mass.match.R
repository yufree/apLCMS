#' An internal function: finding matches between two vectors of m/z values.
#' 
#' Given two vectors of m/z values and the tolerance ppm level, find the
#' potential matches between the two vectors.
#' 
#' 
#' @param x m/z values from the data.
#' @param known.mz m/z values from the known feature table.
#' @param match.tol.ppm tolerance level in ppm.
#' @return A vector the same length as x. 1 indicates matched, and 0 indicates
#' unmatched.
#' @author Tianwei Yu <tyu8@@emory.edu>
#' @keywords models
#' @examples
#' 
#' 
#' 
mass.match <-
function(x, known.mz, match.tol.ppm=5)
{
    mass.matched.pos<-rep(0, length(x))
    for(i in 1:length(x))
    {
        this.d<-abs((x[i]-known.mz)/x[i])
        if(min(this.d) < match.tol.ppm/1e6) mass.matched.pos[i]<-1
    }
    return(mass.matched.pos)
}
