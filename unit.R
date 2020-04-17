library(sandwich)
library(survey)

source("D:/Github/transport-sim/calibrate.R")
source("D:/Github/transport-sim/tmle.R")
source("D:/Github/transport-sim/simfun.R")

iter <- 1000
n_0 <- 500
n_1 <- 500
prob <- 0.3
sig2 <- 2
scen <- "baseline"

# set.seed(06261992)

simDat <- replicate(iter, gen_data(n_1 = n_1, n_0 = n_0, prob = prob, sig2 = sig2, scenario = scen))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat, sparse = (scen == "sparse"))

tau_tmp <- do.call(rbind, estList[1,])
cp_tmp <- do.call(rbind, estList[2,])
colnames(tau_tmp) <- c("GLM", "OUT", "TMLE", "MOM", "ENT")
colnames(cp_tmp) <- c("OUTSATE", "OUTPATE", "ENTSATE", "ENTPATE", "INDPATE")

tau <- apply(tau_tmp, 2, mean, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
cp <- apply(cp_tmp, 2, mean, na.rm = TRUE)
