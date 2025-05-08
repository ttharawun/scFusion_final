## ChiDist_2.2.1_Combine.R (updated)
## - Skips rows with empty counts or expression to avoid NaNs in CalcZINB
## - Replaces `break` with `next` in the likelihood loop
## - Uses L-BFGS-B with rep(0,9) initial parameters (bounds unchanged)
## - Clamps and reuses fitted `lambda` for rnbinom to prevent NA warnings

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript ChiDist_2.2.1_Combine.R <numcell> <ChiDist_file> <filtered_file> <out_file> [param_file]")
}
numcell    <- as.integer(args[1])
data       <- read.table(args[2], sep = "\t", stringsAsFactors = FALSE)
message("Start Reading ChiDist File! Reading ", args[2])

# filter bad GC (if desired)
neiflag <- ifelse(pmax(data$V14, data$V13) <= 0.05 |
                  pmin(data$V14, data$V13) >= 0.95, 1, 0)
data    <- data[neiflag == 0, , drop = FALSE]
row.names(data) <- seq_len(nrow(data))

#–– skip GC correction ––#
gcprelist <- rep(1, 1001)

# load or fit parameters
out_file   <- args[4]
param_file <- if (length(args) >= 5) args[5] else file.path(tempdir(), "ChiDist_pars.RData")

# helper: split bracketed numeric strings
safe_split <- function(s) {
  s2 <- gsub("^\\[|\\]$", "", as.character(s))
  if (is.na(s2) || s2 == "") return(numeric(0))
  as.numeric(strsplit(s2, ",\\s*")[[1]])
}

## ZINB log-likelihood (unchanged)
CalcZINB <- function(par) {
  ll_sum <- 0
  for (i in seq_len(nrow(data))) {
    cc    <- safe_split(data$V4[i])
    expr1 <- safe_split(data$V15[i])
    expr2 <- safe_split(data$V16[i])

    # skip empty rows to avoid NaNs
    if (length(cc) == 0 || length(expr1) == 0 || length(expr2) == 0) {
      next
    }

    lambda   <- exp(par[5])
    gcpre    <- 0
    linpart  <- par[6] + par[8] * mean(expr1) + par[9] * mean(expr2)
    zeroprob <- pmin(exp(linpart) / (1 + exp(linpart)), 0.9999999999999)

    for (j in seq_along(cc)) {
      mu  <- exp(par[1] + par[2] * gcpre + exp(par[3]) * expr1[j] + exp(par[4]) * expr2[j])
      ll_sum <- ll_sum +
        dnbinom(cc[j], size = 1/lambda, prob = 1/(1 + mu * lambda), log = TRUE) +
        log(1 - zeroprob)
      ll_sum <- ll_sum - log(1 - dnbinom(0, size = 1/lambda, prob = 1/(1 + mu * lambda)))
    }

    ll_sum <- ll_sum + (numcell - length(cc)) * log(zeroprob)
    ll_sum <- ll_sum - log(1 - zeroprob^numcell - numcell * (1-zeroprob)*zeroprob^(numcell-1))

    if (is.nan(ll_sum)) {
      next
    }
  }
  -ll_sum
}

## P-value calculation with safe NB sampling
CalcPValueZINB_Fusion_2.0.0_3 <- function(par) {
  Pv <- rep(NA_real_, nrow(data))
  for (j in seq_len(nrow(data))) {
    cc     <- safe_split(data$V4[j])
    expr1  <- safe_split(data$V15[j])
    expr2  <- safe_split(data$V16[j])

    if (sum(cc) <= 3 || data$V3[j] <= 1 ||
        length(cc) == 0 || length(expr1) == 0 || length(expr2) == 0) {
      Pv[j] <- 1
      next
    }
    ave1    <- mean(expr1)/2
    ave2    <- mean(expr2)/2
    gcpre   <- 0
    linpart <- par[6] + par[8] * ave1 + par[9] * ave2
    zeroprob <- pmin(exp(linpart)/(1+exp(linpart)), 0.9995)

    cc[cc > 15] <- 15
    totalread <- sum(cc)

    ## compute mu0 and safe lambda
    mu0    <- exp(par[1] + par[2]*gcpre + exp(par[3])*2*ave1 + exp(par[4])*2*ave2)
    lambda <- max(0.001, (mu0 + mu0^2*exp(par[5]) - mu0) / mu0^2)

    ## safe NB sampling
    size_sim <- 1 / lambda
    prob_sim <- 1 / (1 + mu0 * lambda)
    prob_sim <- pmin(pmax(prob_sim, 1e-8), 1 - 1e-8)
    simuset  <- rnbinom(10000, size = size_sim, prob = prob_sim)
    sim      <- simuset[simuset > 0]
    
    while (length(sim) < 2000) {
      more <- rnbinom(10000, size = size_sim, prob = prob_sim)
      sim  <- c(sim, more[more > 0])
      if (length(sim) >= 2000) break
    }
    sim[sim >= 15] <- 15

    samples <- replicate(
      100,
      sum(sample(sim, size = rbinom(1, numcell, 1 - zeroprob), replace = TRUE))
    )
    Pv[j] <- pnorm(totalread, mean(samples), sd(samples), lower.tail = FALSE)
  }
  Pv
}

# fit or load parameters
if (file.exists(param_file)) {
  message("Loading parameters from ", param_file)
  load(param_file)
  par <- optres$par
} else {
  message("Fitting ZINB model (this may take a while)...")
  optres <- optim(
    par    = rep(0,9),
    fn     = CalcZINB,
    method = "L-BFGS-B",
    control = list(maxit = 10000)
  )
  par <- optres$par
  save(optres, file = param_file)
}

# apply p-value and FDR
data$Pvalue <- CalcPValueZINB_Fusion_2.0.0_3(par)
data$FDR    <- p.adjust(data$Pvalue, method = 'fdr')

# write out
write.table(data, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Done. Output written to ", out_file)
