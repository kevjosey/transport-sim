library(sandwich)
library(survey)
library(geex)

source("D:/Github/transport-sim/calibrate.R")
source("D:/Github/transport-sim/tmle.R")
source("D:/Github/transport-sim/geex.R")
source("D:/Github/transport-sim/simfun.R")

iter <- 100
n_0 <- 500
n_1 <- 500
prob <- 0.5
sig2 <- 2
scen <- "baseline"

# set.seed(06261992)

simDat <- replicate(iter, gen_data(n_1 = n_1, n_0 = n_0, prob = prob, sig2 = sig2, scenario = scen))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat, sparse = (scen == "sparse"), missing = (scen == "missing"))

tau_tmp <- do.call(rbind, estList[1,])
var_tmp <- do.call(rbind, estList[2,])
cp_tmp <- do.call(rbind, estList[3,])
colnames(tau_tmp) <- c("GLM", "OUT", "AIPW", "TMLE", "MOM", "ENT")
colnames(cp_tmp) <- c("MOMSATE", "MOMPATE", "ENTSATE", "ENTPATE", "INDMOMPATE", "INDENTPATE")

tau <- apply(tau_tmp, 2, mean, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
cp <- apply(cp_tmp, 2, mean, na.rm = TRUE)
