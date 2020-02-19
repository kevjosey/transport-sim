#Seth A. Berkowitz, MD MPH"
#Code adapted from Kara Rudolph https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12213#
#Code cannot be used for commericial purposes and all use must attribute original authors#

## site variable needs to be named 'site' and needs to have value 0 for the site where the outcome data is not used and value 1 for the site where the outcome data is used
## z variable needs to be named z and have values 0/1
## y variable is the outcome
## w variables in a dataframe named w and with names w1:wx

tmle <- function(data, S, Z, Y, nsitemodel,  nzmodel, noutmodel) {
  
  n.dat <- nrow(data)
  
  # Calculate components of clever covariate
  glm_cps <- glm(formula = nsitemodel, data = data, family = "binomial")
  cps <- predict(glm_cps, type = "response")
  
  glm_cpz <- glm(formula = nzmodel, data = data, family = "binomial")
  cpz <- predict(glm_cpz, type = "response")
  
  # Calculate clever covariate.
  ps0 <- mean(I(S == 0))
  
  g0w <- ((1 - cpz) * cps) / (1 - cps)
  g1w <- (cpz * cps) / (1 - cps)
  
  h0w <- ((1 - Z) * I(S == 1)) / g0w
  h1w <- (Z * I(S == 1)) / g1w 
  
  ymodel <- glm(formula = noutmodel, family = "gaussian", data = data, subset = S == 1)
  
  data_new0 <- data_new1 <- as.data.frame(data)
  data_new0$Z <- 0
  data_new1$Z <- 1
  
  # Initial prediction.
  q <- cbind(predict(ymodel, type = "link", newdata = data),
             predict(ymodel, type = "link", newdata = data_new0),
             predict(ymodel, type = "link", newdata = data_new1))
  
  epsilon <- coef(glm(Y ~ -1 + offset(q[, 1]) + h0w + h1w, family = "gaussian", subset = S == 1))
  
  # Update initial prediction.
  q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w), epsilon[1] / g0w , epsilon[2] / g1w )
  
  # Get efficient influence curve values for everyone
  tmleest <- mean(q1[, 3][S == 0]) - mean(q1[, 2][S == 0])
  
  return(list(estimate = tmleest))
  
}
