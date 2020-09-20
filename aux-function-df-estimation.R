library(tictoc)

Z = data$R/data$V

z = tail(Z, 1000)


cl <- makeCluster(3)
clusterExport(cl, c(".gl", "fitdist"))

tic()

gl1 <- parApply(cl = cl, z, 2, function(x)(.gl(x, "std")))

toc()

stopCluster(cl)

fBasics::basicStats(gl1)


library(PerformanceAnalytics)
stopCluster(cl)

cl <- makeCluster(3)
clusterExport(cl, c("kurtosis"))

tic()

kurt1 <- as.numeric(kurtosis(z, method = "excess"))

toc()

stopCluster(cl)


fBasics::basicStats(kurt1)

which.min(gl1)
which.max(kurt1)
