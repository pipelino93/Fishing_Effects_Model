### Edited on 10/11/2017...

library(foreign)
library(reshape2)
library(rgdal)
library(maptools)


rm(list = ls())
# Set working directory
setwd("D:\\Dropbox\\Fishing_effects_model\\to_maedhbh")


gear = "all"  ## Pelagic trawls: "PTR", Non-pelagic trawls: "NPT", Hook and line: "HAL", Jig: "JIG", Pots/traps: "POT"

habFeature = "both"  ## Keep "b":biological, "g": geological, or "both"

if(habFeature == "both"){
	habFeatToKeep = c("G", "B")  
}else{
	habFeatToKeep = toupper(habFeature)
	}


if(gear == "all"){
	gearToKeep = c("PTR", "NPT", "HAL", "JIG", "POT")
}else{
	gearToKeep = toupper(gear)
	}

	
### Output file name
year.now = substr(Sys.time(), 1,4)
month.now = substr(Sys.time(), 6,7)
day.now = substr(Sys.time(), 9,10)


outFile = paste("disturbProps_deepSteep_", gear, "Gear_",habFeature, "Struct_", year.now, month.now, day.now,  sep = "")

	
	

# Import data
grid5k = readOGR(dsn = "Fishing_effects_model.gdb", layer = "Grid5k")
grid5k.dat = grid5k@data

fe = read.csv("R_input_tables\\aggregated_fishing_effort_fake_data.csv")
#fe = subset(fe, YEAR <2016)  ## CIA only goes through June of 2015

gt = read.csv("R_input_tables\\GearWidthTable_022316.csv")

recovery_table = read.csv("R_input_tables\\Recovery_table_DeepCorals.csv")





### Convert total fishing to contact adjusted
fe = merge(fe, gt[,c("GearID", "Gear", "SASI_gear", "adjLow","adjMed", "adjHigh")], by.x = "GEARID", by.y = "GearID", all.x = T)



## Gear mods
fe[fe$GEARID %in% 45:53 & fe$YEAR < 2011, ]$adjLow = 1  # Pre 2011 gear change
fe[fe$GEARID %in% 45:53 & fe$YEAR < 2011, ]$adjMed = 1
fe[fe$GEARID %in% 45:53 & fe$YEAR < 2011, ]$adjHigh = 1


fe[fe$GEARID %in% c(6,7,9,10) & fe$YEAR < 2014, ]$adjLow = 1  # Pre Feb 2014 gear change
fe[fe$GEARID %in% c(6,7,9,10) & fe$YEAR < 2014, ]$adjMed = 1
fe[fe$GEARID %in% c(6,7,9,10) & fe$YEAR < 2014, ]$adjHigh = 1

## None of these gears in Jan 2014
#fe[fe$GEARID %in% c(6,7,9,10) & fe$YEAR == 2014 & fe$MONTH == 1, ]$adjLow = 1  # Pre Feb 2014 gear change
#fe[fe$GEARID %in% c(6,7,9,10) & fe$YEAR == 2014 & fe$MONTH == 1, ]$adjMed = 1
#fe[fe$GEARID %in% c(6,7,9,10) & fe$YEAR == 2014 & fe$MONTH == 1, ]$adjHigh = 1




fe[!(fe$Gear %in% gearToKeep), ]$SUM_Shape_ = 0 # Change area of gears not interested in keeping to zero


fe$adjArea = fe$SUM_Shape_ * runif(nrow(fe), min = fe$adjLow, max = fe$adjHigh) ## Random uniform contact adj from min and max



## aggregate fishing adjusted fishing effort
fe.agg = aggregate(adjArea ~ Grid5k_ID + YEAR + MONTH + SASI_gear, data = fe, sum)









grid_order = sort(unique(fe$Grid5k_ID))

# Sediment 
sed = read.dbf("R_input_tables\\sediment_v4_deepCorals.dbf")


# Create sediment matrix for model.  Make sure grid order is same as I_a 
# and keep only sediment areas
sedProps = as.matrix(sed[match(grid_order, sed$Grid5k_ID), 
                         c("mudProp","sandProp","grpeProp","cobProp","bouldProp", "deepProp") ])




subst_types = c("Mud", "Sand", "Gran.Peb", "Cobble", "Boulder", "DeepSteep")

SASI_gears = c("trawl", "longline", "trap")

# Set parameters
nYears = length(unique(fe$YEAR))
nSubAnnual = 12
nGrid = length(unique(fe$Grid5k_ID))
nGear = 3  # Number of SASI gears:  Trawl, Trap, Longline
nSubst = ncol(sedProps) #Number of substrates

gear_types = levels(fe$SASI_gear)

eg = expand.grid(Grid5k_ID=unique(fe$Grid5k_ID), SASI_gear=unique(fe$SASI_gear)) # create all combos of gear and grid cells

m = merge(grid5k.dat, fe.agg, by = "Grid5k_ID")

m$prop = m$adjAre/m$Shape_Area

m$MONTH = as.numeric(as.character(m$MONTH))










# Populate Fishing effort array
F_a = array(NA, dim = c(nYears, nSubAnnual, nGrid, nGear)) #create empty array

year.i = 1  # year index counter

for(year in min(m$YEAR):max(m$YEAR)){
  my = subset(m, YEAR == year)
  for(month.i in 1:12){
    mym = subset(my, MONTH==month.i)
    
    mym = merge(x = eg, y = mym, 
                by.x = c("Grid5k_ID","SASI_gear"), 
                by.y = c("Grid5k_ID", "SASI_gear"), 
                all.x=T)
    
    mym[is.na(mym$prop),]$prop = 0
    
    mym.x = dcast(Grid5k_ID ~ SASI_gear, data=mym, value.var = "prop", fun.aggregate = function(x) sum(x))
    
    mym.x = mym.x[order(mym.x$Grid5k_ID),]
    mym.x = mym.x[,c(gear_types)]
    
    F_a[year.i, month.i, ,] = as.matrix(mym.x)
    
    
  }
  year.i = year.i + 1
}





# Suceptibility 

suscept.f = function(){
	gear.q = matrix(NA, nrow = nGear, ncol = nSubst)

	i = 1
	for(gear in SASI_gears){
	  gear.m = read.csv(paste("R_input_tables\\Susceptibilty_table_deepCoral_", 
							  gear, ".csv", sep=""))
							  
	  gear.m = subset(gear.m, FeatureClass %in% habFeatToKeep) ## Choose geological or biological features
	  
	  gear.m = gear.m[,subst_types]
	    
	  
	  for(column in 1:ncol(gear.m)){
		gear.m[gear.m[,column] %in% 0, column] = 
		  runif(sum(gear.m[,column] %in% 0), min = 0, max = 0.1)
		
		gear.m[gear.m[,column] %in% 1, column] = 
		  runif(sum(gear.m[,column] %in% 1), min = 0.1, max = 0.25)
		
		gear.m[gear.m[,column] %in% 2, column] = 
		  runif(sum(gear.m[,column] %in% 2), min = 0.25, max = 0.5)
		
		gear.m[gear.m[,column] %in% 3, column] = 
		  runif(sum(gear.m[,column] %in% 3), min = 0.5, max = 1)
	  }
	  
	  gear.q[i,] = colMeans(gear.m, na.rm=T)
	  
	  i = i + 1
	  
	}

	gear.q.df = data.frame(SASI_gear = SASI_gears, gear.q)
	names(gear.q.df)[-1] = subst_types



	q_m = as.matrix(gear.q.df[,subst_types])
	return(q_m)
}


#Fishing impacts (I') 
I.prime_a = array(NA, dim = c(nYears, nSubAnnual, nGrid, nSubst))

for(y in 1:nYears){
  for(m in 1:nSubAnnual){
	q_m = suscept.f()  # Get new susceptibility table for each month
    I_m = F_a[y,m,,] %*% q_m
    I.prime_a[y,m,,] = 1-exp(-I_m)
  }
}



# Recovery (rho')


recovery_table = recovery_table[,subst_types]
tau_m = recovery_table[,subst_types] # Make sure sediments are in correct order

recovery.f = function(){
	for(column in 1:ncol(tau_m)){
	  tau_m[recovery_table[,column] %in% 0, column] = 
		runif(sum(tau_m[,column] %in% 0), min = 0, max = 1)
		
	  tau_m[recovery_table[,column] %in% 1, column] = 
		runif(sum(tau_m[,column] %in% 1), min = 1, max = 2)
		
	  tau_m[recovery_table[,column] %in% 2, column] = 
		runif(sum(tau_m[,column] %in% 2), min = 2, max = 5)
		
	  tau_m[recovery_table[,column] %in% 3, column] = 
		runif(sum(tau_m[,column] %in% 3), min = 5, max = 10)
		
	  tau_m[recovery_table[,column] %in% 4, column] = 
		runif(sum(tau_m[,column] %in% 4), min = 10, max = 50)
	}

	tau_v = colMeans(tau_m, na.rm=T) # Average recovery over all habitat features

	rho_v = 1 / (tau_v * nSubAnnual)  # Convert recovery time in years to rates per month
	
	return(rho_v)

}




rho.prime_a = array(NA, dim = c(nYears, nSubAnnual, nGrid, nSubst))




for(y in 1:nYears){
  for(m in 1:nSubAnnual){
	rho_v = recovery.f()  # Get new recovery values for each month
    for(i in 1:nGrid){
		rho.prime_a[y,m,i,] = 1-exp(-rho_v)
    }
  }
}








# Fishing Effects Model function
FishingEffectsModel = function(I.prime_a, rho.prime_a, H_prop_0){
  model_nYears = dim(I.prime_a)[1]
  model_nSubAnnual = dim(I.prime_a)[2]
  model_nGrid = dim(I.prime_a)[3]
  model_nSubst = dim(I.prime_a)[4]
  
  #Make array to hold H
  H_prop = array(dim = c(model_nYears, model_nSubAnnual, model_nGrid, model_nSubst))  
  
  for(y in 1:model_nYears){
    for(m in 1:model_nSubAnnual){
      
      if(y == 1 & m == 1){      # First time step use H_prop_0 for t-1
        prior_state = H_prop_0
      } else if (m == 1){    
        prior_state = H_prop[y-1,model_nSubAnnual,,]
      } else{                 
        prior_state = H_prop[y,m-1,,]
      }
      
      H_from_H = (1-I.prime_a[y,m,,])*prior_state  # undisturbed remaining undisturbed
      H_from_h = (1-prior_state) * (rho.prime_a[y,m,,]) # disturbed recovered to undisturbed
      H_prop[y,m,,] = H_from_H + H_from_h  # Total proportion disturbed
      
    }
  }
  
  return(H_prop)
}  # end function







##### Run model with H =burnin initial conditions


# Define five year burnin in
H_prop_0 = matrix(runif(nGrid*nSubst, 0,1), nrow = nGrid, ncol = nSubst)  # First start with random conditions

# Cycle through 2003-2005 ten times to create burnin starting conditions
for(i in 1:10){
	H_burn = FishingEffectsModel(I.prime_a[1:3,,,], rho.prime_a[1:3,,,], H_prop_0) 
	H_prop_0 = H_burn[3,12,,]  # extract only the last time step as the new intial condition
}
 

##### Run model
H_tot = FishingEffectsModel(I.prime_a, rho.prime_a, H_prop_0)




# Calculate undisturbed Areas
fished_cells = grid5k.dat[match(grid_order, grid5k.dat$Grid5k_ID), ]


undistProps = matrix(NA, ncol = nSubAnnual*nYears, nrow = length(grid_order))

i = 1
for(y in 1:nYears){
  for(m in 1:12){
    undistProps[,i] = rowSums(H_tot[y,m,,]*sedProps)
    i = i + 1
  }
}

undistProps = data.frame(Grid5k_ID = grid_order, undistProps)


disturbProps = data.frame(Grid5k_ID = undistProps[,1], 
                          apply(undistProps[,-1],2,
                                function(x) 1 - x))								
								
								



## Write results to a shapefile
grid5k@data$orderID = 1:nrow(grid5k@data)

grid5k_results = grid5k@data



grid5k_results = merge(grid5k_results, disturbProps, by = "Grid5k_ID", all.x = T)


months = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

counter = 6  # Start counter at first column that has disturbAreas
for(y in 2003:2016){
	for(m in 1:12){
		names(grid5k_results)[counter] = paste(months[m], y, sep = "")
		counter = counter + 1
		}
	}
	

dist.out = grid5k

dist.out@data = grid5k_results[order(grid5k_results$orderID),]


## Create file name and write it to output

writeSpatialShape(dist.out, paste("FE_model_output\\", outFile, ".shp", sep = ""))













