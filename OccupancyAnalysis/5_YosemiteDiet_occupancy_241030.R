
###########################
# integrated occupancy model using remote camera and scat detection data to estimate carnivore occupancy
###########################
# Run the model
library(rjags)
library(jagsUI)
library(R2jags)

setwd("~/OccupancyAnalysis/LongRun_240812")

set.seed(123)
sink("YosemiteOcc_All_2kmgridsnow.txt")
cat(" 
model {
    # derived parameters
    #Predicted effects of covariates on occupancy
    for(s in 1:S){
    for(i in 1:400){
    logit(snow.fx[s,i]) <- a0[s] + a1[s]*((snow.p[i]-snow.mean)/snow.sd) 
    logit(sdcancov.fx[s,i]) <- a0[s] + a3[s]*((sdcancov.p[i]-sdcancov.mean)/sdcancov.sd) 
    logit(cancov.fx[s,i]) <- a0[s] + a2[s]*((cancov.p[i]-cancov.mean)/cancov.sd) 
    }
    }

    # Specify priors
    # occupancy
    mu.occ ~ dunif(-10,10) 
    tau.occ <- 1/(sd.occ^2)
    sd.occ ~ dunif(0,10)

    #detection
    mu.p ~ dunif(-10,10)
    tau.p <-  1/(sd.p^2)
    sd.p ~ dunif(0,10)
    
    #detection
    mu.pscat ~ dunif(-10,10)
    tau.pscat <-  1/(sd.pscat^2)
    sd.pscat ~ dunif(0,10)
    
    for (s in 1:S){
    # occupancy
    a0[s] ~ dnorm(mu.occ, tau.occ)
    a1[s] ~ dnorm(0, 0.01)
    a2[s] ~ dnorm(0, 0.01)
    a3[s] ~ dnorm(0, 0.01)
    
    # detection
    b0[s] ~ dnorm(mu.p, tau.p)
    b1[s] ~ dnorm(0, 0.01)
    b2[s] ~ dnorm(0, 0.01)
    b0scat[s] ~ dnorm(mu.pscat, tau.pscat)
    b1scat[s] ~ dnorm(0, 0.01)
    b2scat[s] ~ dnorm(0, 0.01)
    } # s
    
    # Derived quantities
    for (s in 1:S) {
    occ.fs[s] <- sum(z[,s])       # Number of occupied sites among those studied
    psi.fs[s] <- occ.fs[s]/GG  
    }
    
    #### coyote
    # Observation model for camera detections
    #for (s in 1:S){
    for (r in 1:RR){                    # For each station
    for (t in 1:TT){                    # For each year
    for (m in 1:MM){                    # For each month
    ycoy[r,t,m] ~ dbern(muccoy[r,t,m]*deployed[r,t,m])        # presence/absence each year on cam
    muccoy[r,t,m] <- z[stgrid[r],1]*pcamcoy[r,t,m]  # true occupancy in the grid cell where camera was placed * probability of detection
    logit(pcamcoy[r,t,m]) <- b0[1] + b1[1]*cameffort[r,t,m] + b2[1]*elcamdetstd[r] 
    #} #s
    } #m
    } #t
    } #r
    
    # Observation model for scat locations
    #for(s in 1:S){
    for(x in 1:XX){   # for each detection grid cell
    for(t in 1:TT){   # for each year
    for(m in 1:MM){   # for each month
    yscatcoy[x,t,m] ~ dbern(mu3coy[x,t,m]) 
    mu3coy[x,t,m] <- z[gridpix[x],1]*pscatcoy[x,t,m]                        
    logit(pscatcoy[x,t,m]) <- b0scat[1] + b1scat[1]*effortseason[x,t,m] + b2scat[1]*eldetstd[x]
    #} #s
    } #m
    } #t
    } #x
    
    #### bobcat
    # Observation model for camera detections
    #for (s in 1:S){
    for (r in 1:RR){                    # For each station
    for (t in 1:TT){                    # For each year
    for (m in 1:MM){                    # For each month
    ybob[r,t,m] ~ dbern(mucbob[r,t,m]*deployed[r,t,m])        # presence/absence each year on cam
    mucbob[r,t,m] <- z[stgrid[r],2]*pcambob[r,t,m]  # true occupancy in the grid cell where camera was placed * probability of detection
    logit(pcambob[r,t,m]) <- b0[2] + b1[2]*cameffort[r,t,m] + b2[2]*elcamdetstd[r] 
    #} #s
    } #m
    } #t
    } #r
    
    # Observation model for scat locations
    #for(s in 1:S){
    for(x in 1:XX){   # for each detection grid cell
    for(t in 1:TT){   # for each year
    for(m in 1:MM){   # for each month
    yscatbob[x,t,m] ~ dbern(mu3bob[x,t,m]) 
    mu3bob[x,t,m] <- z[gridpix[x],2]*pscatbob[x,t,m]                        
    logit(pscatbob[x,t,m]) <- b0scat[2] + b1scat[2]*effortseason[x,t,m] + b2scat[2]*eldetstd[x]
    #} #s
    } #m
    } #t
    } #x
    
    #### cougar
    # Observation model for camera detections
    #for (s in 1:S){
    for (r in 1:RR){                    # For each station
    for (t in 1:TT){                    # For each year
    for (m in 1:MM){                    # For each month
    ycoug[r,t,m] ~ dbern(muccoug[r,t,m]*deployed[r,t,m])        # presence/absence each year on cam
    muccoug[r,t,m] <- z[stgrid[r],3]*pcamcoug[r,t,m]  # true occupancy in the grid cell where camera was placed * probability of detection
    logit(pcamcoug[r,t,m]) <- b0[3] + b1[3]*cameffort[r,t,m] + b2[3]*elcamdetstd[r] 
    #} #s
    } #m
    } #t
    } #r
    
    # Observation model for scat locations
    #for(s in 1:S){
    for(x in 1:XX){   # for each detection grid cell
    for(t in 1:TT){   # for each year
    for(m in 1:MM){   # for each month
    yscatcoug[x,t,m] ~ dbern(mu3coug[x,t,m]) 
    mu3coug[x,t,m] <- z[gridpix[x],3]*pscatcoug[x,t,m]                        
    logit(pscatcoug[x,t,m]) <- b0scat[3] + b1scat[3]*effortseason[x,t,m] + b2scat[3]*eldetstd[x]
    #} #s
    } #m
    } #t
    } #x
    
    #### fox
    # Observation model for camera detections
    #for (s in 1:S){
    for (r in 1:RR){                    # For each station
    for (t in 1:TT){                    # For each year
    for (m in 1:MM){                    # For each month
    yfox[r,t,m] ~ dbern(mucfox[r,t,m]*deployed[r,t,m])        # presence/absence each year on cam
    mucfox[r,t,m] <- z[stgrid[r],4]*pcamfox[r,t,m]  # true occupancy in the grid cell where camera was placed * probability of detection
    logit(pcamfox[r,t,m]) <- b0[4] + b1[4]*cameffort[r,t,m] + b2[4]*elcamdetstd[r] 
    #} #s
    } #m
    } #t
    } #r
    
    # Observation model for scat locations
    #for(s in 1:S){
    for(x in 1:XX){   # for each detection grid cell
    for(t in 1:TT){   # for each year
    for(m in 1:MM){   # for each month
    yscatfox[x,t,m] ~ dbern(mu3fox[x,t,m]) 
    mu3fox[x,t,m] <- z[gridpix[x],4]*pscatfox[x,t,m]                        
    logit(pscatfox[x,t,m]) <- b0scat[4] + b1scat[4]*effortseason[x,t,m] + b2scat[4]*eldetstd[x]
    #} #s
    } #m
    } #t
    } #x
    
    #### marten
    # Observation model for camera detections
    #for (s in 1:S){
    for (r in 1:RR){                    # For each station
    for (t in 1:TT){                    # For each year
    for (m in 1:MM){                    # For each month
    ymart[r,t,m] ~ dbern(mucmart[r,t,m]*deployed[r,t,m])        # presence/absence each year on cam
    mucmart[r,t,m] <- z[stgrid[r],5]*pcammart[r,t,m]  # true occupancy in the grid cell where camera was placed * probability of detection
    logit(pcammart[r,t,m]) <- b0[5] + b1[5]*cameffort[r,t,m] + b2[5]*elcamdetstd[r] 
    #} #s
    } #m
    } #t
    } #r
    
    # Observation model for scat locations
    #for(s in 1:S){
    for(x in 1:XX){   # for each detection grid cell
    for(t in 1:TT){   # for each year
    for(m in 1:MM){   # for each month
    yscatmart[x,t,m] ~ dbern(mu3mart[x,t,m]) 
    mu3mart[x,t,m] <- z[gridpix[x],5]*pscatmart[x,t,m]                        
    logit(pscatmart[x,t,m]) <- b0scat[5] + b1scat[5]*effortseason[x,t,m] + b2scat[5]*eldetstd[x]
    #} #s
    } #m
    } #t
    } #x

    # Ecological submodel for occupancy
    for (s in 1:S){  # species s
    for (g in 1:GG){ # grid cell g
    z[g,s] ~ dbern(psi[s,g]) # Occupancy is bernoulli distributed (0/1)
    logit(psi[s,g]) <- a0[s] + a1[s]*snow[g] + a2[s]*cancov[g] + a3[s]*sdcancov[g] 
    } #g
    } #s
    
    }
    ",fill = TRUE)
sink()

# Bundle data
#data_YNP <- list(ycoy = YNPCoyote, 
#                 ybob = YNPBob,
#                 ymart = YNPMarten,
#                 ycoug = YNPCoug,
#                 yfox = YNPFox,
#                 GG = nrow(grid), 
#                 RR = nrow(stgrid), 
#                 XX = nrow(griddet),
#                 S = dim(YNPCarn)[1], 
#                 stgrid = as.vector(stgrid[,1]),
#                 gridpix = gridingridpix$gridingridpix,
#                 TT = dim(YNPCarn)[3],
#                 MM = dim(YNPCarn)[4],
#                 yscatcoy = coyscat_season*season_sampled2, 
#                 yscatfox = foxscat_season*season_sampled2, 
#                 yscatbob = bobscat_season*season_sampled2, 
#                 yscatcoug = cougscat_season*season_sampled2,
#                 yscatmart = martscat_season*season_sampled2,
#                 effortseason = effort_seasonstd2,
#                 cameffort = cameffort_std,
#                 sampledseason = season_sampled2,
#                 sdcancov = sdcancovstd,
#                 cancov = cancovstd,
#                 snow = snpckstd, 
#                 snow.mean = snow.mean,
#                 snow.sd = snow.sd,
#                 snow.p = snow.p,
#                 sdcancov.mean = sdcancov.mean,
#                 sdcancov.sd = sdcancov.sd,
#                 sdcancov.p = sdcancov.p,
#                 cancov.mean = cancov.mean,
#                 cancov.sd = cancov.sd,
#                 cancov.p = cancov.p,
#                 elcamdetstd = elcamdetstd,
#                 eldetstd = eldetstd[,1],
#                 deployed = deployed)
#
#### save for later
data_YNP <- readRDS("YosemiteOcc_jagsdata_240814.txt")
##### Initial values
zst_YNP <- array(1, c(nrow(grid_poly), dim(YNPCarn)[1]))

# Initial values
jags_inits1 <- list(".RNG.seed" = 1, ".RNG.name" = "base::Wichmann-Hill", z = array(1, dim = c(dim(grid_poly)[1], dim(YNPCarn)[1])))
jags_inits2 <- list(".RNG.seed" = 2, ".RNG.name" = "base::Wichmann-Hill", z = array(1, dim = c(dim(grid_poly)[1], dim(YNPCarn)[1])))
jags_inits3 <- list(".RNG.seed" = 3, ".RNG.name" = "base::Wichmann-Hill", z = array(1, dim = c(dim(grid_poly)[1], dim(YNPCarn)[1])))

inits_YNP <- list(jags_inits1,
                 jags_inits2,
                 jags_inits3)
# Parameter
params_YNP <- c("a0",
                "a1",
                "a2",
                "a3",
                "b0",
                "b1",
                "b2",
                "b0scat",
                "b1scat",
                "b2scat",
                "occ.fs",
                "psi.fs",
                "cancov.fx",
                "sdcancov.fx",
                "snow.fx",
                "psi",
                "z")

# Run the model
library(rjags)
library(jagsUI)
library(R2jags)
Sys.time() 
YNP_Occ_alllargegridsnow <- jagsUI::jags(data_YNP, 
                 inits = inits_YNP, 
                 params_YNP, 
                 "YosemiteOcc_All_2kmgridsnow.txt", 
                 n.chains = 3,
                 n.thin = 5,
                 n.iter = 75000,
                 n.adapt = NULL,
                 n.burnin = 50000,
                 parallel = TRUE,
                 #save.all.iter = TRUE,
                 #Rhat.limit = 1.05,
                 n.cores = 12,
                 verbose = TRUE)
Sys.time()

