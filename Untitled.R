disney <- read.csv("disneyland.csv")

### QUESTION 1 ###

# a)

pop <- list(pop1 = disney[disney$Reviewer_Location == "Canada" & 
                            disney$Branch == "Disneyland_California", ], 
            pop2 = disney[disney$Reviewer_Location == "Australia" & 
                            disney$Branch == "Disneyland_California", ])

par(mfrow = c(1, 2))
hist(pop[[1]]$Rating, col=adjustcolor("firebrick", 0.7), freq = FALSE,
     xlab = "Rating", main = "Canada Residents", breaks = 0:5, ylim = c(0, 0.7))
hist(pop[[2]]$Rating, col=adjustcolor("blue", 0.7), freq = FALSE,
     xlab = "Rating", main = "Australia Residents", breaks = 0:5, ylim = c(0, 0.7))

# c) 

sdN <- function(y) {
  ybar = mean(y)
  sqrt(mean((y-ybar)^2))
}
y_test <- c(4,9,3,2,7)
sdN_test <- sdN(y_test)
print(sdN_test)

# d)

skew <- function(y) {
  ybar = mean(y)
  mean((y-ybar)^3) / sdN(y)^3
}
skew_test <- skew(y_test)
print(skew_test)

# e)

kurt <- function(y) {
  ybar = mean(y)
  mean((y-ybar)^4) / sdN(y)^4
}
kurt_test <- kurt(y_test)
print(kurt_test)

# f)

mixRandomly <- function(pop) {
  pop1 <- pop$pop1
  n_pop1 <- nrow(pop1)
  
  pop2 <- pop$pop2
  n_pop2 <- nrow(pop2)
  
  mix <- rbind(pop1, pop2)
  select4pop1 <- sample(1:(n_pop1 + n_pop2), n_pop1, replace = FALSE)
  
  new_pop1 <- mix[select4pop1, ]
  new_pop2 <- mix[-select4pop1, ]
  list(pop1 = new_pop1, pop2 = new_pop2)
}

calculatePVmulti <- function(pop, discrepancies, M_outer, M_inner) {
  
  ## Local function to calculate the significance levels over the discrepancies
  ## and return their minimum
  
  getPVmin <- function(basePop, discrepancies, M) {
    observedVals <- sapply(discrepancies, FUN = function(discrepancy) {
      discrepancy(basePop)
    })
    
    K <- length(discrepancies)
    
    total <- Reduce(function(counts, i) {
      # mixRandomly mixes the two populations randomly, so the new sub-populations
      # are indistinguishable
      NewPop <- mixRandomly(basePop)
      
      ## calculate the discrepancy and counts
      Map(function(k) {
        Dk <- discrepancies[[k]](NewPop)
        if (Dk >= observedVals[k]) 
          counts[k] <<- counts[k] + 1
      }, 1:K)
      counts
    }, 1:M, init = numeric(length = K))
    
    PVs <- total/M
    min(PVs)
  }
  
  PVmin <- getPVmin(pop, discrepancies, M_inner)
  
  total <- Reduce(function(count, m) {
    basePop <- mixRandomly(pop)
    if (getPVmin(basePop, discrepancies, M_inner) <= PVmin) 
      count + 1 else count
  }, 1:M_outer, init = 0)
  
  PVstar <- total/M_outer
  return(PVstar)
}

D1 <- function(pop){
  abs(mean(pop[[1]]$Rating)-mean(pop[[2]]$Rating))
}

D2 <- function(pop){
  abs(sdN(pop[[1]]$Rating)/sdN(pop[[2]]$Rating)-1)
}

D3 <- function(pop){
  abs(skew(pop[[1]]$Rating)/skew(pop[[2]]$Rating)-1)
}

D4 <- function(pop){
  abs(kurt(pop[[1]]$Rating)/kurt(pop[[2]]$Rating)-1)
}

p_discrepancies <- list(D1, D2, D3, D4)
cache = TRUE
pval <- calculatePVmulti(pop, p_discrepancies, M_outer = 100, M_inner = 100)
print(pval) # 0.02

### QUESTION 2 ###

bootstrap_t_interval <- function(S, a, confidence, B, D) {
  ## Inputs: 
  ##    S = an n element array containing the variate values in the sample 
  ##    a = a scalar-valued function that calculates the attribute a() of interest 
  ##    confidence = a value in (0,1) indicating the confidence level 
  ##    B = a numeric value representing the outer bootstrap count of
  ##    replicates (used to calculate the lower and upper limits) 
  ##    D = a numeric value representing the inner bootstrap count of replicates
  ##    (used to estimate the standard deviation of the sample attribute for
  ##    each (outer) bootstrap sample)
  
  Pstar <- S
  aPstar <- a(Pstar)
  sampleSize <- length(S)
  ## get (outer) bootstrap values
  bVals <- sapply(1:B, FUN = function(b) {
    Sstar <- sample(Pstar, sampleSize, replace = TRUE)
    aSstar <- a(Sstar)
    ## get (inner) bootstrap values to estimate the SD
    Pstarstar <- Sstar
    SD_aSstar <- sd(sapply(1:D, FUN = function(d) {
      Sstarstar <- sample(Pstarstar, sampleSize, replace = TRUE)
      ## return the attribute value
      a(Sstarstar)
    }))
    z <- (aSstar - aPstar)/SD_aSstar
    ## Return the two values
    c(aSstar = aSstar, z = z)
  })
  SDhat <- sd(bVals["aSstar", ])
  zVals <- bVals["z", ]
  ## Now use these zVals to get the lower and upper c values.
  cValues <- quantile(zVals, probs = c((1 - confidence)/2, (confidence + 
                                                              1)/2), na.rm = TRUE)
  cLower <- min(cValues)
  cUpper <- max(cValues)
  interval <- c(lower = aPstar - cUpper * SDhat, middle = aPstar, upper = aPstar - 
                  cLower * SDhat)
  return(interval)
}

# a)
canada <- subset(disney, disney$Reviewer_Location == "Canada")$Rating
sampIndex <- read.table("sampIndex.txt")$V1
S <- canada[sampIndex]
summary(S)

# b)
B <- 1000
n <- 100
Sstar <- sapply(1:B, FUN = function(b) {
  sample(S, n, replace = TRUE)
})

# c)

# i)
as <- mean(S)
ap <- mean(canada)

# ii)
as_star <- apply(X = Sstar, MARGIN = 2, FUN = mean)
hist(as_star, col = adjustcolor("darkgreen", 0.5), xlab = "Mean Rating",
     main = "1000 Bootstrap Replicates")
abline(v = ap, col = "purple", lwd = 2)

# iii)
a_ci_naive <- as + 1.96 * c(-1, 1) * sd(as_star)

# iv)
a_ci_per <- c(quantile(as_star, 0.025), quantile(as_star, 0.975))

# v)
a_ci_boot <- bootstrap_t_interval(S = S, a = mean, confidence = 0.95,
                                  B = 1000, D = 100)


# d)

# i)
# sdN function from question 1
sds <- sdN(S)
sdP <- sdN(canada)

# ii)
sds_star <- apply(X = Sstar, MARGIN = 2, FUN = sdN)
hist(sds_star, col = adjustcolor("darkgreen", 0.5), xlab = "Std Dev Rating",
     main = "1000 Bootstrap Replicates")
abline(v = sdP, col = "purple", lwd = 2)

# iii)
sd_ci_naive <- sds + 1.96 * c(-1, 1) * sd(sds_star)

# iv)
sd_ci_per <- c(quantile(sds_star, 0.025), quantile(sds_star, 0.975))

# v)
sd_ci_boot <- bootstrap_t_interval(S = S, a = sdN, confidence = 0.95,
                                   B = 1000, D = 100)


# e)

# i)
# skew function from question 1
skews <- skew(S)
skewP <- skew(canada)

# ii)
skews_star <- apply(X = Sstar, MARGIN = 2, FUN = skew)
hist(skews_star, col = adjustcolor("darkgreen", 0.5), xlab = "Skewness of Rating",
     main = "1000 Bootstrap Replicates")
abline(v = skewP, col = "purple", lwd = 2)

# iii)
skew_ci_naive <- skews + 1.96 * c(-1, 1) * sd(skews_star)

# iv)
skew_ci_per <- c(quantile(skews_star, 0.025), quantile(skews_star, 0.975))

# v)
skew_ci_boot <- bootstrap_t_interval(S = S, a = skew, confidence = 0.95,
                                   B = 1000, D = 100)


# f)

# i)
# kurt function from question 1
kurts <- kurt(S)
kurtP <- kurt(canada)

# ii)
kurts_star <- apply(X = Sstar, MARGIN = 2, FUN = kurt)
hist(kurts_star, col = adjustcolor("darkgreen", 0.5), xlab = "Kurtosis of Rating",
     main = "1000 Bootstrap Replicates")
abline(v = kurtP, col = "purple", lwd = 2)

# iii)
kurt_ci_naive <- kurts + 1.96 * c(-1, 1) * sd(kurts_star)

# iv)
kurt_ci_per <- c(quantile(kurts_star, 0.025), quantile(kurts_star, 0.975))

# v)
kurt_ci_boot <- bootstrap_t_interval(S = S, a = kurt, confidence = 0.95,
                                     B = 1000, D = 100)



