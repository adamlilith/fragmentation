#' Calculate fragmentation indices for a continuous raster
#'
#' The function calculates a set of fragmentation indices as per Riitters et al. (2000). The index developed by those authors is based on the assumption that forest cover can be represented in binary fashion (foerst/not forest). This function accomodates situations where the raster represents non-binary data (e.g., percent forest). This is done by thresholding the forest raster using a series of thresholds, then calculating density, connectivity, and fragmentation after thresholding for each, and finally calculating aggregate indices. For density and connectivity, this is the average value across thresholds. For fragmentation class, each cell is scored by the percentage of times it is assigned to each of the 6 fragmentation classes. A composite index can then be calculated by taking the class the cell was most often assigned to.\cr
#'
#' See documentation for \code{\link{fragBinary}} and \href{https://www.jstor.org/stable/26271763}{Riitters et al. 2020} for definitions of the six fragmentation classes.
#'
#' @param rast \code{SpatRaster} with non-binary values.
#' @param thresholds Numeric vector or a single numerical value. This is either a vector of numerical values representing the actual thresholds to use, or a single value, in which case it represents the number of equally-spaced thresholds that will be used to calculate fragmentation. If the latter, then the spacing will start at the smallest value on the raster (often 0) and proceed in equally-sized bins to the maximim value (often 100).
#' @param incIndividual Logical. If \code{FAlSE} (default), just return the composite rasters calculated across all thresholds. If \code{TRUE}, the return the composite rasters plus the rasters generated for each threshold.
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
#' @returns A multi-layered \code{SpatRaster} with these rasters (depending on the values of \code{calcDensity}, \code{calcConnect}, and \code{calcClass}):
#' \itemize{
#'	\item Rasters named "\code{densityMean}" and "\code{densitySD}": Mean and standard deviation of density across thresholds.
#'	\item Rasters named "\code{connectMean}" and "\code{connectSD}": Mean and standard deviation of connectivity across thresholds.
#'	\item Rasters named "\code{mostCommonClass}" and "\code{mostCommonClassFreq}": Most common class assigned to each cell and the frequency the most common class was assigned.
#' }
#' Additionally, if \code{incIndividual} is \code{TRUE}, then there will be one to three more rasters per threshold. These will be named like \code{density_}*, \code{connect_}*, and \code{class_}* where the * is the value of the threshold used to calculate the metrics. See the help page for \code{\link{fragBinary}} for information on how to interpret these values.
#'
#' @seealso \code{\link{fragBinary}} for binary-valued forest raster data
#'
#' @example man/examples/fragmentation_examples.r
#'
#' @export

fragCont <- function(
	rast,
	thresholds,
	incIndividual = FALSE,
	w = 3,
	pad = FALSE,
	padValue = NA,
	calcDensity = TRUE,
	calcConnect = TRUE,
	calcClass = TRUE,
	na.rm = FALSE,
	undet = 'undetermined',
	cores = 1,
	verbose = FALSE,
	...
) {

	# debugging block
	if (FALSE) {

		thresholds <- 4
		incIndividual <- FALSE
		w <- 3
		pad <- FALSE
		padValue <- NA
		calcDensity <- TRUE
		calcConnect <- TRUE
		calcClass <- TRUE
		na.rm <- FALSE
		undet <- 'undetermined'
		cores <- 4
		verbose <- TRUE


	}

	# thresholds
	mm <- terra::minmax(rast)
	if (length(thresholds) == 1L) {
		start <- mm[1L, 1L]
		end <- mm[2L, 1L]
		thresholds <- seq(start, end, length.out=thresholds)
	} else {
		if (any(thresholds < mm[1L, 1L])) stop('At least one threshold is lower than the lowest value in the raster.')
		if (any(thresholds > mm[2L, 1L])) stop('At least one threshold is higher than the highest value in the raster.')
	}

	# for each threshold evaluate fragmentation
	for (i in seq_along(thresholds)) {

		threshold <- thresholds[i]

		if (verbose) omnibus::say('Evaluating fragmentation for threshold ', round(threshold, 5), '.')

		trast <- rast >= threshold
		# frag <- forestFrag(rast=trast, ...)
		frag <- fragBinary(rast=trast, 	w = 3, pad = FALSE, padValue = NA, calcDensity = TRUE, calcConnect = TRUE, calcClass = TRUE, na.rm = FALSE, undet = 'undetermined', cores = 4) # for debugging

		# remember
		if (calcDensity) {
			if (exists('density', inherits=FALSE)) {
				density <- c(density, frag[['density']])
			} else {
				density <- frag[['density']]
			}
		}

		if (calcConnect) {
			if (exists('connect', inherits=FALSE)) {
				connect <- c(connect, frag[['connect']])
			} else {
				connect <- frag[['connect']]
			}
		}

		if (calcClass) {
			if (exists('class', inherits=FALSE)) {
				class <- c(class, frag[['class']])
			} else {
				class <- frag[['class']]
			}
		}

	} # next threshold

	### calculate integrated values
	if (calcDensity) {

		densityMean <- app(density, 'mean', na.rm=TRUE)
		names(densityMean) <- 'densityMean'

		densitySD <- app(density, 'sd', na.rm=TRUE)
		names(densitySD) <- 'densitySD'

		out <- c(densityMean, densitySD)

	}

	if (calcConnect) {

		connectMean <- app(connect, 'mean', na.rm=TRUE)
		names(connectMean) <- 'connectMean'

		connectSD <- app(connect, 'sd', na.rm=TRUE)
		names(connectSD) <- 'connectSD'

		if (exists('out', inherits=FALSE)) {
			out <- c(out, connectMean, connectSD)
		} else {
			out <- c(connectMean, connectSD)
		}

	}

	if (calcClass) {

		mostCommonClass <- app(class, fun=.mostCommonClass, cores=cores)
		mostCommonClassFreq <- app(class, fun=.mostCommonClassFreq, cores=cores)

		names(mostCommonClass) <- 'mostCommonClass'
		names(mostCommonClassFreq) <- 'mostCommonClassFreq'

		out <- c(out, mostCommonClass, mostCommonClassFreq)

	}

	# include fragmentation rasters for each threshold
	if (incIndividual) {

		names(density) <- paste0('density_', thresholds)
		names(connect) <- paste0('connect_', thresholds)
		names(class) <- paste0('class_', thresholds)

		out <- c(out, density, connect, class)

	}

	out

}

#' Calculates the most common fragmentation class across a set of cells
.mostCommonClass <- function(x) {

	tab <- table(x)
	mostCommon <- names(tab[which.max(tab)])
	mostCommon <- as.numeric(mostCommon)
	mostCommon

}

#' Calculates the average frequency of the most common fragmentation class across a set of cells
.mostCommonClassFreq <- function(x) {

	tab <- table(x)
	mostCommon <- tab[which.max(tab)]
	mostCommon / sum(tab)

}
