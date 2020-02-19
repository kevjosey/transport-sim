library(sandwich)
library(survey)

source("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Additional Code/transport-sim/calibrate.R")
source("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Additional Code/transport-sim/tmle.R")
source("D:/Dropbox (ColoradoTeam)/JoseyDissertation/Additional Code/transport-sim/simfun.R")

iter <- 1000
n_0 <- 1000
n_1 <- 1000
prob <- 0.3
sig2 <- 5
scen <- "baseline"

# set.seed(06261992)

simDat <- replicate(iter, gen_data(n_1 = n_1, n_0 = n_0, prob = prob, sig2 = sig2, scenario = scen))

idx <- 1:iter # simulation iteration index
estList <- sapply(idx, simfit, simDat = simDat, sparse = (scen == "sparse"))

tau_tmp <- do.call(rbind, estList[1,])
cp_tmp <- do.call(rbind, estList[2,])
colnames(tau_tmp) <- c("GLM", "OUT", "TMLE", "MOM", "ENT")
colnames(cp_tmp) <- c("ENTSATE", "ENTPATE", "OUTSATE", "OUTPATE", "INDPATE")

tau <- apply(tau_tmp, 2, mean, na.rm = TRUE)
mcse <- apply(tau_tmp, 2, sd, na.rm = TRUE)
cp <- apply(cp_tmp, 2, mean, na.rm = TRUE)
