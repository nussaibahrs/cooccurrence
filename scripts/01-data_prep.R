# Load libraries and functions --------------------------------------------
library(chronosphere)
library(divDyn)

source(file.path("scripts", "functions.R"))

# Load data ---------------------------------------------------------------
dir.create(file.path("data", "chronodat"))
pbdb <- fetch("pbdb", ver="20211029", datadir = file.path("data", "chronodat"))

# Reef  environment -----------------------------------------------------------
#reefs <- getReefs(pbdb)

data(keys)
pbdb$reefs <- categorize(pbdb$environment, keys$reefs)
pbdb$bath <- categorize(pbdb$environment, keys$bath)

# Tropics -----------------------------------------------------------------
keys$lat$t <- c(0, 35) #35 degrees

pbdb$tropical <- categorize(abs(pbdb$paleolat), keys$lat)

# reefs only
reefs <- subset(pbdb, reefs=="reef")

# shallow tropical non-reef
nonreefs <- subset(pbdb, reefs=="non-reef")
nonreefs <- subset(nonreefs, bath=="shallow")
nonreefs <- subset(nonreefs, tropical=="t")

reefs <- rbind(reefs, nonreefs)
table(reefs$reefs)

# bin data
data(keys)
data(tens)


# the 'stg' entries (lookup)
stgMin <- categorize(reefs[ ,"early_interval"], keys$tenInt)
stgMax <- categorize(reefs[ ,"late_interval"], keys$tenInt)

# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

# empty container
reefs$stg <- rep(NA, nrow(reefs))

# select entries, where
stgCondition <- c(
	# the early and late interval fields indicate the same stg
	which(stgMax==stgMin),
	# or the late_intervar field is empty
	which(stgMax==-1))

# in these entries, use the stg indicated by the early_interval
reefs$ten[stgCondition] <- stgMin[stgCondition] 

save(reefs, file = file.path("data", "reefs.Rdata")) #pbdb data
