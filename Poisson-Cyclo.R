
######################################################################################################
## This R script has been created to experiment with ways to integrate geochronologic dates ("nails") and durations ("springs") into an age-depth model or a geologic time scale
## by dr. David De Vleeschouwer 
######################################################################################################

rm(list=ls()) # Empty global environment
install.packages(c("astrochron", "Hmisc", "Bchron", "devtools"))
library(astrochron) # Load some useful R libraries that must be installed
library(Hmisc)
library(Bchron)
library(devtools)

######################################################################################################
# Let's take a look at the dummy dataset.
######################################################################################################
par(mar = c(4,4,1,1))
dates = matrix(c(0,9,34,46,66,72,0,44.8,116,129,523,797.3,0.0001,10,20,40,120,140), nrow = 6, ncol = 3)
colnames(dates) <- c("Depth (m)", "Age (ka)", "Uncertainty (2sigma, kyr)")
dates = as.data.frame(dates)
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim=c(0,1000),xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}

######################################################################################################
# [1] We use a Poisson distribution to simulate break-points in-between dated levels where sedimentation rate changes. 
# Both in depth and time-domain, break-points are sampled from a uniform distribution. 
# This is undoubtedly over-pessimistic because this allows the Monte Carlo procedure to deviate far away from the straight line between two subsequent dates leves. But, let's keeps things simple for now. 
######################################################################################################
Nsimulations = 10000
ages_lie_MonteCarlo_Poisson=matrix(data = NA, nrow = 101, ncol = Nsimulations ) # This matrix will contain all individual Monte-Carlo-simulated age-depth chronologies 
depths_i=seq(0,100,1) # We'll make an age-depth model with a 1-meter resolution

for (i in 1:Nsimulations){
  age_depth_MonteCarlo=c()
  depth_Poisson=c()
  age_Poisson=c()
  
  for (j in 1:6){age_depth_MonteCarlo[j]=rnorm(1, mean = dates$`Age (ka)`[j], sd = dates$`Uncertainty (2sigma, kyr)`[j]/2)}
    if (length(which(diff(age_depth_MonteCarlo)<0))>0){ages_lie_MonteCarlo_Poisson[,i] = rep(NA, 101) # Prior knowledge: Law of superposition
    }  else{
      for (k in 1:5){
        thickness = dates$`Depth (m)`[k+1]-dates$`Depth (m)`[k]
        lambda = thickness / 20
        N_breaks = rpois(1, lambda = lambda)
        temp1=runif(N_breaks, min = dates$`Depth (m)`[k], max = dates$`Depth (m)`[k+1]) # temp1 represents a depth-level at which a break-point is introduced. temp1 is sampled from a uniform distribution between the two adjacent dated levels
        temp2=runif(N_breaks, min = age_depth_MonteCarlo[k], max = age_depth_MonteCarlo[k+1]) # temp2 represents the age of the introduced break-point. temp2 is sampled from a uniform distribution between the two adjacent dates
        depth_Poisson = c(depth_Poisson, temp1) 
        age_Poisson = c(age_Poisson, temp2)
        }
      depths_sim = c(dates$`Depth (m)`, depth_Poisson) # We now concatenate the depths of the dated levels and the introduced break-points for this specific Monte-Carlo simulation
      ages_sim = c(age_depth_MonteCarlo, age_Poisson) # We now concatenate the ages of the dated levels and the introduced break-points for this specific Monte-Carlo simulation
      temp=approxExtrap(depths_sim, ages_sim , xout = depths_i)$y # Temp is the age-depth model (through linear interpolation) of this specific Monte-Carlo simulation
      ages_lie_MonteCarlo_Poisson[,i]=temp # ages_lie_MonteCarlo_Poisson is a large matrix, storing each of the 1000 age-depth models that are generated during Monte-Carlo. 
      }
  }

ages_lie_MC_Poisson=matrix(data = NA, nrow = 101, ncol = 3) # This matrix will hold the Median Age and the Confidence Levels
colnames(ages_lie_MC_Poisson) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
  ages_lie_MC_Poisson[k,]=quantile(ages_lie_MonteCarlo_Poisson[k,], probs = c(0.025,0.5,0.975), na.rm = T)}

# Let's plot things up
dev.off()
par(mar = c(4,4,1,1))
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_MC_Poisson[,2], depths_i, col = "grey")
lines(ages_lie_MC_Poisson[,1], depths_i) # Upper confidence level
lines(ages_lie_MC_Poisson[,3], depths_i) # Lower confidence level
text(1300,30,"Uncertainties strongly depend on the adopted Poisson model", col = "red", cex = 0.75)
text(1300,34,"Here, lambda is set at 20.", col = "red", cex = 0.7)
text(1300,37,"This means that sedimentation rates", col = "red", cex = 0.7)
text(1300,40,"are expected to change every 20 meters.", col = "red", cex = 0.7)


######################################################################################################
# [2] We use a Poisson distribution to simulate break-points in-between dated levels where sedimentation rate changes. 
# In the depth domain, break-points are sampled from a uniform distribution. 
# But, instead of using a uniform distribution to determine the ages of the break-points, we use the cyclostratigraphy-derived sedimentation rates. 
######################################################################################################
Nsimulations = 10000
ages_lie_MonteCarlo_Poisson=matrix(data = NA, nrow = 101, ncol = Nsimulations )
depths_i=seq(0,100,1)

sedrate_cyclo = c(rep(4,55),rep(20,20), rep(30,40)) # These are the "dummy" cyclostratigraphy-derived sedimentation rates. 

# When we plot them up as an age-depth model, it becomes clear that they are largely congruent with the dated levels, but provide several new pieces of information. ----
cyclo_age_model = c(0)
for (i in 2:101){
  cyclo_age_model[i]=cyclo_age_model[i-1]+sedrate_cyclo[i]
}
dev.off()
par(mar = c(4,4,1,1))
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
lines(cyclo_age_model, depths_i, col = "blue", lwd = 2)

# Now let's get serious ----
# We use the same philosophy as in [1], using a uniform distribution to determine the depth-positions of the break-points.
# But, instead of using a uniform distribution to determine the ages of the break-points, we use the cyclostratigraphy-derived sedimentation rates. 
# The age of the break point is determined by applying the cyclostratigraphy-derived sedimentation rate AT the break-point to the depth-interval in-between the break-point and the dated level BELOW the break-point. 
# The "AT" and "BELOW" are user-dependent choices, and other choices could be made here (e.g. average sedimentation rates over a certain interval). But the "AT" and "BELOW" choice produce the most reasonable results in this artificial case.

for (i in 1:Nsimulations){
  age_depth_MonteCarlo=c()
  depth_Poisson=c()
  age_Poisson=c()
  
  for (j in 1:6){age_depth_MonteCarlo[j]=rnorm(1, mean = dates$`Age (ka)`[j], sd = dates$`Uncertainty (2sigma, kyr)`[j]/2)}
  if (length(which(diff(age_depth_MonteCarlo)<0))>0){ages_lie_MonteCarlo_Poisson[,i] = rep(NA, 101) # Prior knowledge: Law of superposition
  } else{
    for (k in 1:6){ 
      if (k==6){thickness = 100-dates$`Depth (m)`[k]
        } else{thickness = dates$`Depth (m)`[k+1]-dates$`Depth (m)`[k]}
      lambda = thickness / 20
      N_breaks = rpois(1, lambda = lambda)
      if (k==6){temp1=runif(N_breaks, min = dates$`Depth (m)`[k], max = 100)
      } else{temp1=runif(N_breaks, min = dates$`Depth (m)`[k], max = dates$`Depth (m)`[k+1])} # temp1 represents a depth-level at which a break-point is introduced. temp1 is sampled from a uniform distribution between the two adjacent dated levels
      temp2 = dates$`Age (ka)`[k]+(temp1-dates$`Depth (m)`[k])*sedrate_cyclo[ceiling(temp1)] # temp2 represents the age of the introduced break-point. temp2 is calculated using information on sedimentation rate from cyclostratigraphy
      depth_Poisson = c(depth_Poisson, temp1) 
      age_Poisson = c(age_Poisson, temp2) 
    }
    depths_sim = c(dates$`Depth (m)`, depth_Poisson) # We now concatenate the depths of the dated levels and the introduced break-points for this specific Monte-Carlo simulation
    ages_sim = c(age_depth_MonteCarlo, age_Poisson) # We now concatenate the ages of the dated levels and the introduced break-points for this specific Monte-Carlo simulation. 
    temp=approxExtrap(depths_sim, ages_sim , xout = depths_i)$y # Temp is the age-depth model (through linear interpolation) of this specific Monte-Carlo simulation
    ages_lie_MonteCarlo_Poisson[,i]=temp # ages_lie_MonteCarlo_Poisson is a large matrix, storing each of the 1000 age-depth models that are generated during Monte-Carlo. 
  }
}

ages_lie_MC_Poisson=matrix(data = NA, nrow = 101, ncol = 3)
colnames(ages_lie_MC_Poisson) <- c("2.5% CL", "Median Age", "97.5% CL")
for (k in 1:101){
  ages_lie_MC_Poisson[k,]=quantile(ages_lie_MonteCarlo_Poisson[k,], probs = c(0.025,0.5,0.975), na.rm = T)}

dev.off()
par(mar = c(4,4,1,1))
plot(dates$`Age (ka)`, dates$`Depth (m)`, ylim = c(100,0), xlim = c(0, 2500), xaxs = "i", yaxs = "i", col = "blue", pch = 19, xlab = "Age (ka)", ylab = "Depth (m)")
for (j in 1:6){segments(dates[j,2]-dates[j,3],dates[j,1],dates[j,2]+dates[j,3],dates[j,1], col = "blue", lwd = 2)}
points(ages_lie_MC_Poisson[,2], depths_i, col = "grey")
lines(ages_lie_MC_Poisson[,1], depths_i) # Upper confidence level
lines(ages_lie_MC_Poisson[,3], depths_i) # Lower confidence level
text(1300,30,"Uncertainties strongly depend on the adopted Poisson model", col = "red", cex = 0.75)
text(1300,34,"Here, lambda is set at 20.", col = "red", cex = 0.7)
text(1300,37,"Here, this means that sedimentation rates", col = "red", cex = 0.7)
text(1300,40,"are expected to change every 20 meters.", col = "red", cex = 0.7)
lines(cyclo_age_model, depths_i, col = "blue", lwd = 2)


