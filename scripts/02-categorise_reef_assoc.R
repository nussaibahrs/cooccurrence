library(data.table)
library(dplyr)
library(gt)
library(gtsummary)

# Associations from pared -------------------------------------------------
biota_associations <- readxl::read_excel(file.path("data", "pareddat", "Legend_Biota_Detailed.xlsx"))

#Some cleaning
biota_associations$ASSOCIATION <- gsub("only", "", biota_associations$ASSOCIATION, ignore.case = T)
biota_associations$ASSOCIATION <- gsub("tabulate/rugose corals|rugose/tabulate corals", "tabulate corals - rugose corals ", biota_associations$ASSOCIATION,
									   ignore.case = T)
biota_associations$ASSOCIATION <- gsub("tabulate/rugose corals|rugose/tabulate corals", "tabulate corals - rugose corals ", biota_associations$ASSOCIATION,
									   ignore.case = T)
biota_associations$ASSOCIATION <- gsub("Non-", "non", biota_associations$ASSOCIATION)
biota_associations$ASSOCIATION <- gsub("/", "-", biota_associations$ASSOCIATION)
biota_associations$ASSOCIATION <- gsub("stromatoporoids or chaetetids", "stromatoporoids - chaetetids", biota_associations$ASSOCIATION, ignore.case = T)
biota_associations$ASSOCIATION <- tolower(trimws(gsub("\\(.+\\)", "", biota_associations$ASSOCIATION))) #remove everything between brackets
biota_associations$ASSOCIATION <- gsub("-$", "", biota_associations$ASSOCIATION)
biota_associations$ASSOCIATION <- gsub(", ", "-", biota_associations$ASSOCIATION)

# View(biota_associations)

# get list of organisms per association
assoc <- data.table::rbindlist(
	lapply(strsplit(biota_associations$ASSOCIATION, "-"), function(x) as.data.frame(matrix(x, nrow=1))),
	fill=T)
assoc$ass_num <- biota_associations$ASS_NUM
assoc_long <- na.omit(reshape2::melt(assoc, id.vars="ass_num", value.name="group"))[,c(1,3)]
assoc_long$group <- trimws(assoc_long$group)
assoc_long <- assoc_long[-which(assoc_long$group==""),]
# View(assoc_long)

table(assoc_long$group)

# recategorisations
assoc_long$group[grep("foram", assoc_long$group, ignore.case = TRUE)] <- "forams"

#Algae
assoc_long$group[assoc_long$group %in% c("solenoporaceans", "palaeoaplysina", "coralline algae")] <- "red algae"

alg <- c("algae", "phylloid and other platy algae", "noncoralline algae",       
		 "green algae" ,"phylloid algae","spygoria",         
		 "receptaculid algae" )
assoc_long$group[assoc_long$group %in% alg] <- "other algae"

# Corals
assoc_long$group[grep("corals", assoc_long$group)] <- "corals"

#Calcareous sponges
calc_sp <- c("stromatoporoids", "coralline sponges", "sponges", "pharetronid sponges",
			 "pharetronids", "spongiomorphids", "archaeocyaths", "chaetetids",   "calcisponges" , "calcareous sponges", "sclerosponges")

# "lithistid sponges" -> silicieous

assoc_long$group[assoc_long$group %in% calc_sp] <- "calcareous sponges"
assoc_long$group[grep("sphinctozoan", assoc_long$group)] <- "calcareous sponges"


# Microbes
micr<- c("cyanoyphyceans","cyanophyceans", "tubiphytes", "stromatolites", "girvanella",
		 "microbial micrite", "calathium", "stromatolites")
assoc_long$group[assoc_long$group %in% micr] <- "microbes"

# Bivalves
biv <- c("oysters", "nonrudist bivalves", "bivalves", "massive bivalves")
assoc_long$group[assoc_long$group %in% biv] <- "other bivalves"

#Brachiopods
assoc_long$group[grep("brachiopod", assoc_long$group)] <- "brachiopods"

# Vermetids
assoc_long$group[assoc_long$group == "vermetid gastropods"] <- "vermetids"

#View(assoc_long)

# Create table
original <-readxl::read_excel(file.path("data", "pareddat", "Legend_Biota_Detailed.xlsx"))
original <- original[,1:2]
colnames(original) <- tolower(colnames(original))


groups <- assoc_long %>% group_by(ass_num) %>% 
	mutate(group = stringr::str_to_sentence(group)) %>% 
	summarise(group = paste(group, collapse=", "))

gtsummary::theme_gtsummary_compact()

tab1 <- merge(original, groups) %>% 
	gt %>% 
	fmt_markdown(columns = everything()) %>% 
	cols_label(ass_num=md("**Code**"),
			   association=md("**Association details**"),
			   group = md("**Standardised Groups**")) %>% 
	cols_width(
		ass_num ~ px(50),
		association ~ px(350),
		group ~ px(200)
	)

tab1
gtsave(tab1, "output/table_s1.rtf")
