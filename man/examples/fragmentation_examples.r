library(terra)

# forest cover
rastFile <- system.file('extdata', 'rocForest.tif', package='fragmentation')
roc <- rast(rastFile)

# fragmentation: binary forest data
rocBinary <- roc > 80 # threshold to make it bonary

binary <- fragBinary(rocBinary, cores=2)
plot(c(rocBinary, binary))

# fragmentation: continuous values with equally-spaced thresholds
contEqual <- fragCont(roc, thresholds=4, cores=2, verbose=TRUE)
plot(contEqual)

# fragmentation: continuous values with user-defined thresholds
tr <- c(70, 80, 90) # user-defined thresholds
contUser <- fragCont(roc, thresholds=tr, cores=2, verbose=TRUE)
plot(contUser)
