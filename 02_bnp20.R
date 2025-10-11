
#####################################################
######## BNP model for ordinal drought data #########
#####################################################


code_start = Sys.time()
sink(file = "runtime.txt", append = TRUE)
cat(paste("code started at: ", Sys.time(), "\n"))
cat("\n")
sink()

library(nimble) ## for MCMC computation using NIMBLE.
library(coda) ## for manipulation of MCMC results.

set.seed(810)

###### run as an array job; run grid cell K in parallel
args = Sys.getenv('SLURM_ARRAY_TASK_ID')
K = as.numeric(args[1])
outfile = paste0("BNP20/out",K,".Rda")

###### load data
select <- read.csv("USDMData.csv") ## these csv files were generated when previously running ScalingValues.R
scalingvalues <- read.csv("scalingvalues.csv")

###### chosen to preserve grid cells with some missing stream28 data
II = which(select$time<20201231 & select$time>20110101)
select = select[II,]
  
###### select all weeks
len <- length(unique(select$time))

###### remove grids with no change in drought status
maxim = rep(0,len)
for (j in 1:length(unique(select$grid))){
  maxim[j] = max(table(select$drought[select$grid == unique(select$grid)[j]]))
}
unique(select$grid[which(maxim==len)])
grid.rem = c("N78","W98","GG14","WW88")
select = select[-which(select$grid %in% grid.rem),]

###### remove a smattering of grid cells with no stream28
stream.rem = c("S9", "S10", "T10", "PP49", "RR50", "SS51", "TT52", "UU53")

unique(select$grid)
hold = select[!duplicated(select$grid), ]
coords = hold[,c("grid","lon","lat")]
Q = nrow(coords);Q
Tobs = nrow(select)/Q;Tobs

sub = select[select$grid == unique(select$grid)[K],]

rm(select)

###### scale covariates
vars <- c("evp","stream28")
vars <- sort(vars)
Xcov = cbind(sub[,vars])
J = ncol(Xcov)
Xm=rep(0,J)
Xs=rep(0,J)
for (j in 1:J){
  Xm[j] <- scalingvalues$means[which(scalingvalues$X %in% vars[j])]
  Xs[j] <- scalingvalues$sds[which(scalingvalues$X %in% vars[j])]
}

X = matrix(NA,Tobs,J)
for(j in 1:J){
  X[,j] = (Xcov[,j]-Xm[j])/Xs[j]
}

X = cbind(rep(1,Tobs),X)
colnames(X) = c("int",vars)

bp = ncol(X)

###### make Y
hold = factor(sub$drought,levels=c("0","D0","D1","D2","D3","D4"))
yfull = as.numeric(hold)-1

I = length(unique(sub$grid))
crdtm = data.frame("lon"=sub$lon[1:I],"lat"=sub$lat[1:I])

###### H is a data frame with 522 rows and 20 columns, representing 10 pairs of Fourier sine and cosine basis functions.
numB <- 20
H <- read.csv("h20.csv")[,-1]

###### stuff for NIMBLE
mod_data=list(Y=yfull, X=X, H=H)

mod_constants=list(Tobs=Tobs, bp=bp, cut=c(0,1,2,3,4), L=50, numB=numB)

mod_inits=list(beta=rep(0,bp),
               tau.ztilde = rep(1,mod_constants$L),
               Z=yfull-0.5,
               alpha=1,
               cl = sample(1:10, size = mod_constants$Tobs, replace = TRUE), 
               beta0tilde = rep(0,mod_constants$L),
               gamma = rep(0,numB))

model_code=nimbleCode({
  
  ## Drought Variable
  for(t in 1:Tobs){
    Y[t] ~ dinterval(Z[t], cut[])
  }

  for (i in 1:Tobs) {
    tau.z[i] <- tau.ztilde[cl[i]]
    beta0[i] <- beta0tilde[cl[i]]
  }
  
  ## Latent Gaussian Variable
  for(t in 1:Tobs){
    mu[t] <- beta0[t] + inprod(X[t,2:bp],beta[2:bp]) + inprod(H[t,1:numB],gamma[1:numB])
    Z[t] ~ dnorm(mu[t], tau = tau.z[t])
  }
  
  ## CRP
  cl[1:Tobs] ~ dCRP(alpha, size = Tobs) ## cluster vector 
  alpha ~ dgamma(1, 1)
  
  for(i in 1:L) {
    beta0tilde[i] ~ dnorm(0, tau=0.04)
    tau.ztilde[i] ~ dgamma(shape=0.01, rate=0.01)
  }
  
  for (b in 2:bp){
    beta[b] ~ dnorm(0, tau = 0.04)
  }
  
  for (b in 1:numB){
    gamma[b] ~ dnorm(0, tau = 0.04)
  }
  
} ## closes NIMBLE code
)

nimble_model <- nimbleModel(code=model_code, constants=mod_constants, data=mod_data, inits=mod_inits)

compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)

monitors <- c('beta', 'tau.ztilde', 'cl', 'beta0tilde', 'alpha', 'gamma')
mcmc_conf <- configureMCMC(nimble_model,monitors=monitors,
                           control=list(adaptive=TRUE,scale=0.1,adaptInterval=100000,sliceMaxSteps=100000,maxContractions=5000000,sliceAdaptWidthMaxIter=0,sliceAdaptFactorMaxIter=0),
                           useConjugacy = TRUE,enableWAIC = TRUE)

nimble_mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(nimble_mcmc, project = nimble_model,resetFunctions = TRUE)

MCS=100000
nthin = 10
last <- MCS/(2*nthin)
samples_sub = runMCMC(compiled_mcmc,inits=mod_inits,
                    nchains = 1, nburnin=MCS/2,niter = MCS,samplesAsCodaMCMC = TRUE,thin=nthin,
                    summary = FALSE, WAIC = TRUE, progressBar=TRUE)

sampG <- getSamplesDPmeasure(compiled_mcmc)

samp <- samples_sub$samples
g.samp <- samp[ , grep('gamma', colnames(samp))]

weightIndex <- grep('weights', colnames(sampG[[1]]))

beta0TildeIndex <- grep('beta0tilde', colnames(sampG[[1]]))

tauTildeIndex <- grep('tau.ztilde', colnames(sampG[[1]]))


###### computation of the full-sample log score for the BNP model with basis functions.
lpp <- rep(0, Tobs)
for(j in 1:Tobs){
  
  if(yfull[j]==0){
    temp <- 0
    for (t in 1:length(sampG)) {
      temp <- temp + sum(sampG[[t]][,weightIndex]*
                           pnorm(0, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                   inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                 sd = 1/sqrt(sampG[[t]][,tauTildeIndex])))
    }
    lpp[j] <- log(temp/length(sampG))
  }
  else if(yfull[j]==1){
    temp <- 0
    for (t in 1:length(sampG)) {
      temp <- temp + sum(sampG[[t]][,weightIndex]*
                           (pnorm(1, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                    inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                  sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))-
                              pnorm(0, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                      inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                    sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))))
    }
    lpp[j] <- log(temp/length(sampG))
  }
  else if(yfull[j]==2){
    temp <- 0
    for (t in 1:length(sampG)) {
      temp <- temp + sum(sampG[[t]][,weightIndex]*
                           (pnorm(2, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                    inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                  sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))-
                              pnorm(1, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                      inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                    sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))))
    }
    lpp[j] <- log(temp/length(sampG))
  }
  else if(yfull[j]==3){
    temp <- 0
    for (t in 1:length(sampG)) {
      temp <- temp + sum(sampG[[t]][,weightIndex]*
                           (pnorm(3, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                    inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                  sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))-
                              pnorm(2, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                      inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                    sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))))
    }
    lpp[j] <- log(temp/length(sampG))
  }
  else if(yfull[j]==4){
    temp <- 0
    for (t in 1:length(sampG)) {
      temp <- temp + sum(sampG[[t]][,weightIndex]*
                           (pnorm(4, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                    inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                  sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))-
                              pnorm(3, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                      inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                    sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))))
    }
    lpp[j] <- log(temp/length(sampG))
  }
  else if(yfull[j]==5){
    temp <- 0
    for (t in 1:length(sampG)) {
      temp <- temp + sum(sampG[[t]][,weightIndex]*
                           (1-pnorm(4, mean = sampG[[t]][,beta0TildeIndex]+samp[t,"beta[2]"]*X[j,2]+samp[t,"beta[3]"]*X[j,3]+
                                      inprod(H[j,1:numB],g.samp[t,1:numB]), 
                                    sd = 1/sqrt(sampG[[t]][,tauTildeIndex]))))
    }
    lpp[j] <- log(temp/length(sampG))
  }
}

logscore = sum(lpp)/Tobs

grid = sub$grid[1]
SimMCMC = list(logscore, samples_sub$WAIC, grid, crdtm, K, MCS, 
               nthin, Tobs, bp, yfull, samp, sampG)

save(SimMCMC, file=outfile)

sink(file = "runtime.txt", append = TRUE)
cat(paste("code ended at: ", Sys.time(), "\n"))
cat("\n")
print(difftime(Sys.time(), code_start))
cat("\n")
sink()
