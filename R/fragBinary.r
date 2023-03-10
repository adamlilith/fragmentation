#' Calculate fragmentation indices for a binary raster
#'
#' The function calculates a set of fragmentation indices as per Riitters et al. (2000).
#'
#' @param rast \code{SpatRaster} with binary values (1 or 0 or \code{NA}).
#' @param w Integer, number of cells wide and high of the window used to calculate fragmentation. This must be an odd integer (default is 3).
#' @param pad Logical, if \code{TRUE} then add virtual rows and columns around the raster so that there are no edge effects. The virtual rows and columns are set to equal \code{padValue}. Default is \code{FALSE}.
#' @param padValue Value to which to set the values of the virtual cells used to pad the raster if \code{pad} is \code{TRUE}.
#' @param calcDensity Logical, if \code{TRUE} (default) then calculate density raster.
#' @param calcConnect Logical, if \code{TRUE} (default) then calculate connectivity raster.
#' @param calcClass Logical, if \code{TRUE} (default) then calculate classification raster. Note that to calculate the classification raster the density and connectivity rasters must also be calculated (\code{calcDensity} and \code{calcConnect} should both be \code{TRUE}). If they are not then the will be forced to \code{TRUE} with a warning.
#' @param na.rm Logical, if \code{FALSE} (default) then \code{NA} cells count as part of the area potentially occupied in a window (i.e., the count in the denominator when calculating density and they are counted as potential links in the connectance calculations if a neighboring cell has a value of 1). If \code{FALSE} then areas that border \code{NA} cells could still be classified as "interior" or otherwise have less apparent fragmentation if the occupied cells are fully surrounded by other occupied cells (except for the \code{NA} cells).
#' @param undet Character. When classifying this defines what is done with "undetermined" cases (when density is >= 0.6 and density == connectivity). Possible values include (partial matching of strings is used):
#' \itemize{
#' 	\item \code{'undetermined'}: Undetermined cases will be assigned a value of 5 (which is not assigned to any other case; default).
#' 	\item \code{'perforated'}: Undetermined cases will be assigned a value of 3 ("perforated").
#' 	\item \code{'edge'}: Undetermined cases will be assigned a value of 4 ("edge").
#' 	\item \code{'random'}: Undetermined cases will be assigned a value of 3 or 4 at random ("perforated" or "edge").
#' }
#' @param ... Other arguments (not used).
#'
#' @references Riitters, K., J. Wickham, R. O'Neill, B. Jones, and E. Smith. 2000. Global-scale patterns of forest fragmentation. Conservation Ecology 4:3. \href{https://www.jstor.org/stable/26271763}{URL}. Also note the \href{https://www.ecologyandsociety.org/vol4/iss2/art3/errata/january26.2001.html}{errata} to the paper on their classification scheme.
#'
#' @returns A \code{SpatRaster} "stack" with three layers: 1) a fragmentation classification (named \code{class}); 2) the density of "1" pixels in the window (named \code{density}--called "pf" in Riitter et al. 2000); 3) and a connectivity raster (conditional probability a cell with a value of 1 has a value that is also 1 named \code{connect}--called "pff" in Riitter et al. 2000). The density and connectivity rasters have values in the range [0, 1], but the classification raster has coded values (from the erratum to Ritter et al. (2000):
#' \itemize{
#' 	\item \code{NA}: \code{NA}
#' 	\item \code{0}: No forest (or whatever is being evaluated)
#'	\item \code{1}: patch (\code{pf} < 0.4)
#'	\item \code{2}: transitional (0.4 <= \code{pf} < 0.6)
#'	\item \code{3}: perforated (\code{pf} >= 0.6 & \code{pf - pff} > 0)
#'	\item \code{4}: edge (\code{pf} >= 0.6 & \code{pf - pff} < 0)
#'	\item \code{5}: undetermined (\code{pf} >= 0.6 & \code{pf == pff})
#'	\item \code{6}: interior (\code{pf} == 1)
#' }
#' Note that this differs somewhat from the numbering scheme presented by Riitters et al. (2000) and their errata.
#' @param cores Number of processor cores to use. The default is to use 1.
#'
#' @seealso \code{\link{fragCont}} for continuously-valued forest raster data
#'
#' @example man/examples/fragmentation_examples.r
#'
#' @export

fragBinary <- function(
	rast,
	w = 3,
	pad = FALSE,
	padValue = NA,
	calcDensity = TRUE,
	calcConnect = TRUE,
	calcClass = TRUE,
	na.rm = FALSE,
	undet = 'undetermined',
	cores = 1,
	...
) {

	if (w %% 2 == 0 | w < 3) stop('Argument "w" to function fragmentation() must be an odd integer >= 3.')
	if (calcClass & (!calcDensity | !calcConnect)) {
		calcDensity <- calcConnect <- TRUE
		warning('Forcing "calcDensity" and "calcConnect" in function "fragmentation()" to be TRUE since "calcClass" is TRUE.')
	}
	
	halfWindowSize <- (w - 1) / 2
	
	# add padding around raster to avoid edge effects
	if (pad) {
		origExtent <- terra::ext(rast)
		rast <- terra::extend(rast, y=halfWindowSize, fill=padValue)
	}

	w <- matrix(1, nrow=w, ncol=w)
	
	# calculate fragmentation indices
	if (calcDensity) pf <- terra::focal(rast, w=w, fun=.fragDensity, na.rm=na.rm)
	if (calcConnect) pff <- terra::focal(rast, w=w, fun=.fragConnect, na.rm=na.rm)
	
	# calculate class
	if (calcClass) {
		out <- c(pf, pff)
		names(out) <- c('density', 'connect')
		classification <- terra::app(out, fun=.fragClassify, undet=undet, cores=cores)
		names(classification) <- 'class'
		out <- c(classification, out)
	}
	
	terra::crs(out) <- terra::crs(rast)
	if (pad) out <- terra::crop(out, origExtent)	
	out
	
}
