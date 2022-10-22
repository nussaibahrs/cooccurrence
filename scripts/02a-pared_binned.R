library(chronosphere)
library(divDyn)
library(tidyverse)

source(file.path("scripts", "functions.R"))
pal <- c("#5f0f40", "#9a031e", "#fb8b24", "#e36414", "#0f4c5c")

pared <- read.csv(file.path("data", "pared_data.csv"))
pared <- subset(pared, type !=3) # remove mud mound
pared <- pared[!(pared$biota_main==3 & pared$biota_sec==3),] # remove microbial mats only
pared <- pared[pared$max_ma <= 541,] # Phanerozoic only
nrow(pared)

pared %>% dplyr::select(r_number,thickness, width, extension) %>% 
  pivot_longer(cols=thickness:extension) %>% 
  na.omit() %>% 
  group_by(name) %>% 
  tally() %>% 
  mutate(prop=n/nrow(pared)*100)

pared %>% dplyr::select(r_number, m_thick, m_width, m_ext) %>% 
  pivot_longer(cols=m_thick:m_ext) %>% 
  na.omit() %>% 
  group_by(name) %>% 
  tally() %>% 
  mutate(prop=n/nrow(pared)*100)

intervals <- read.csv(file.path("data", "pareddat", "l_intervall_new.csv"))
pbdb_reef_bin <- read.csv(file.path("data", "pareddat", "Reef_bins_new.csv"))

intervals <- merge(intervals, pbdb_reef_bin[,c("Bin", "PBDB_Bin")], 
	  by.x="reef_bin", by.y="Bin")
intervals <- intervals[!is.na(intervals$PBDB_Bin),]

int_summary <- gts(2020, "ten")

pared <- merge(pared, intervals, by.x="intervall", by.y="Interval")
pared$bin <- pared$PBDB_Bin - 1

# Associations from pared -------------------------------------------------
source(file.path("scripts", "02-categorise_reef_assoc.R"))

dat <- pared[,c("r_number", "reef_bin", "biota_deta", "pal_lat_scotese", "pal_long_scotese")]
colnames(dat) <- c("reefs", "bin", "ass_num", "paleolat", "paleolong")

reefs <- merge(dat, assoc_long, by.x="ass_num", by.y="ass_num") #add group
save(reefs, file="data/reefs_group.RData")


# Thickness ---------------------------------------------------------------
# est <- tapply(pared$m_thick, pared$thickness, mean, na.rm=T)
# 
# n.thick <- rep(NA, nrow(pared))
# n.thick <- pared$thickness
# n.thick[pared$thickness == 1] <- round(est[1])
# n.thick[pared$thickness == 2] <- round(est[2])
# n.thick[pared$thickness ==3] <- round(est[3])
# n.thick[pared$thickness ==4] <- round(est[4])
# 
# pared$n.thick <- n.thick
# 
# # Mean thickness per time bin
# n.thick <- data.frame(thick=tapply(pared$n.thick, 
#                                    pared$reef_bin, mean, na.rm=T))
# n.thick$bin <- as.numeric(rownames(n.thick))

# Number of reefs ---------------------------------------------------------
n.thick <- data.frame(thick=tapply(pared$r_number, 
								   pared$reef_bin, length))
n.thick$bin <- as.numeric(rownames(n.thick))


# Volume ------------------------------------------------------------------

est <- tapply(pared$m_thick, pared$thickness, mean, na.rm=T)

n.thick <- numeric(nrow(pared))
n.thick <- pared$thickness
n.thick[pared$thickness == 1] <- round(est[1])
n.thick[pared$thickness == 2] <- round(est[2])
n.thick[pared$thickness ==3] <- round(est[3])
n.thick[pared$thickness ==4] <- round(est[4])
n.thick[!is.na(pared$m_thick)] <- pared$m_thick[!is.na(pared$m_thick)]

#pared$n.thick <- n.thick

# 2. Extent
# check best estimates for metric extent # could be more conservative
est <- tapply(pared$m_ext, pared$extension, mean, na.rm=T)
n.ext <- numeric(nrow(pared))
n.ext <- pared$extension
n.ext[pared$extension == 1] <- round(est[1])
n.ext[pared$extension == 2] <- round(est[2])
n.ext[pared$extension ==3] <- round(est[3])
n.ext[pared$extension ==4] <- round(est[4])
n.ext[!is.na(pared$m_ext)] <- pared$m_ext[!is.na(pared$m_ext)]

#pared$n.ext <- n.ext

# 3. Width
# check best estimates for metric extent # could be more conservative
est <- tapply(pared$m_width, pared$width, mean, na.rm=T)
n.width <- numeric(nrow(pared))
n.width <- pared$width
n.width[pared$width == 1] <- round(est[1])
n.width[pared$width == 2] <- round(est[2])
n.width[pared$width ==3] <- round(est[3])
n.width[pared$width ==4] <- round(est[4])
n.width[!is.na(pared$m_width)] <- pared$m_ext[!is.na(pared$m_width)]

#pared$n.width <- n.width


# Also include pared for which dimensions are completely unknown
pared$ext_belt[is.na(pared$ext_belt)] <- 0

# get the three dimensions separtely
# Assume 10 m when thickness is not
# thick <- pared$n.thick
n.thick[is.na(n.thick)] <- 10

# extent <- pared$n.ext
n.ext[is.na(n.ext)] <- 100   # was 300    

# width <- pared$n.width
n.width[is.na(n.width)] <- n.thick[is.na(n.width)]  # take thickness estimate for width   


# now modify extent by extent of reef tract
# need to find a way to handle the "belonging" entries
extent <- n.ext
# For short belts take half of the estimate directly
extent[pared$ext_belt>0 & pared$ext_belt<=10] <- pared$ext_belt[pared$ext_belt>0 & pared$ext_belt<=10]*1000/2 # 
# Longer belts
extent[pared$ext_belt>10 & pared$ext_belt<=25] <- 12000
extent[pared$ext_belt>25 & pared$belonging!=""] <- 18000
extent[pared$ext_belt>25 & pared$belonging=="" & pared$ext_belt<=100] <- 25000
extent[pared$ext_belt>100 & pared$belonging==""] <- 100000

# Consider cases when actually observed extent is greater than inferred from tract length
extent[extent<n.ext] <- n.ext[extent<n.ext]

vol1 <- n.thick*n.width*n.ext
vol <- n.thick*n.width*extent

pared <- cbind(pared, vol, vol1, n.thick)


