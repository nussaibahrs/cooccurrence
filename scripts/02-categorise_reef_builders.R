
# Load data ---------------------------------------------------------------
load(file.path("data","reefs.Rdata")) #pbdb data
microbe <- readxl::read_excel(file.path("data","calcimicrobes.xls"))

table(reefs$phylum)
# Categories --------------------------------------------------------------
reefs$group <- NA

# Corals
# Scleractinian

reefs$group[reefs$order=="Scleractinia"] <- "scleractinian corals"

t.gen <- c("Turnacipora", "Pachypseudoplasmopora", "Procatenipora", "Moorowipora", "Columnopora",
           "Tetraporella", "Ligulodictyum")
t.ord <- c("Tetradiida","Sarcinulida","Heliolitida", "Favositida", "Auloporida", "Lichenariida")


dat <- data.table::rbindlist(list(dat, cbind.data.frame(group="corals", subgroup="scleractinian corals", 
				 genus=paste(t.gen, collapse=", "), order=paste(t.ord, collapse=", "))),
				 fill=T)
#Tabulate

reefs$group[reefs$genus %in% t.gen | 
              reefs$order %in% t.ord | reefs$accepted_name == "Tabulata"] <- "tabulate corals"

r.gen <- c("Stauromatidium", "Yokoyamaella (Maoriphyllum)", "Stereodepasophyllum", "Dinostrophinx",
           "Pseudoamplexus", "Omphyma", "Prodarwinia", "Pseudomicroplasma", "Radiastrea", "Botrophyllum",
           "Sinospongophyllum", "Timanophyllum", "Ferganophyllum", "Duplocarinia", "Cystihexagonaria",
           "Enallophrentis", "Puanophyllum", "Frechastrea", "Nephelophyllum", "Diversiphyllum",
           "Parvaxon", "Leptelasma", "Orthophyllum", "Ketophyllum", "Xiangzhouphyllum", "Yishanophyllum",
           "Hamarophyllum", "Nardophyllum", "Mixogonaria", "Pseudohexagonaria", "Cuctienophyllum",
           "Smithicyathus", "Pantophyllum", "Piceaphyllum", "Agastophyllum", "Parasunophyllum",
           "Scissoplasma", "Thryptophyllum", "Molophyllum", "Parasiphonophyllia", "Solominella",
           "Cystodactylon", "Gudmania")
r.fam <- c("Protozaphrentidae", "Chonophyllidae") 
r.ord <- c("Stauriida", "Cystiphyllida")
reefs$group[reefs$genus %in% r.gen | reefs$family %in% r.fam |
              reefs$order %in% r.ord | reefs$accepted_name == "Rugosa"] <- "rugose corals"

reefs$group[reefs$order == "Heterocorallia"] <- "heterocorals" # could be grouped with rugosans but not yet

# Unknown corals, Hydrozoans and scyphozoans (no structural integrity)
h.cla <- c("Scyphozoa", "Hydrozoa", "Conulata", "Hydroconozoa")
h.ord <- c("Gorgonacea", "Helioporacea", "Alcyonacea", "Pennatulacea") # octocorals
h.ord <- c(h.ord, "Tabulaconida")
h.fam <- c("Helioporidae")
reefs$group[reefs$class %in% h.cla | reefs$order %in% h.ord | reefs$accepted_name == "Anthozoa" |reefs$family %in% h.fam] <- "other cnidarians"

# What about these?
not_classified <- reefs[is.na(reefs$group) & reefs$phylum=="Cnidaria",] #how to deal with these?
unique(not_classified$order)

# Sponges
ker.or <- c("Ancorinida", "Clavulina", "Dictyoceratida", "Haplosclerida", "Poecilosclerida", "Dystactospongia")

reefs$group[reefs$class=="Irregulares" | reefs$order %in% c("Archaeocyathida", "Ajacicyathida", "Monocyathida", "Capsulocyathida", 
							   "Kazakhstanicyathida", "Hetairacyathida", "Acanthinocyathida") | 
				reefs$family %in% c("Archaeocyathidae")] <- "archaeocyaths"

reefs$group[reefs$class == "Stromatoporoidea"] <- "stromatoporoids s.str."

reefs$group[reefs$order == "Axinellida" | reefs$family=="Burgundiidae"] <- "Mesozoic stromatoporoids"

reefs$group[reefs$genus=="Lovcenipora" | reefs$order=="Chaetetida" | reefs$order=="Tabulospongida"] <- "chaetetids"

# Needs cleaning in the PBDB
reefs$group[reefs$class=="Heteractinida"] <- "pharetronids"
reefs$group[reefs$order%in%c("Aspiculata", "Sphinctozoa", "Pharetronida", "Sphaerocoeliida", "Vaceletida",
							 "Permosphincta", "Permosphincta", "Stellispongiida", "Agelasida", "Porata")] <- "pharetronids"
reefs$group[reefs$family%in%c("Disjectoporidae", "Minchinellidae", "Tabasiidae")] <- "pharetronids"


reefs$group[reefs$order%in%c("Lithistida", "Hexatinellida", "Hexactinosa", "Hadromerida", "Lyssacinosa", 
							 "Lychniscosa", "Lyssakida", "Epipolasida",
							 "Innaecoeliida","Reticulosa","Protomonaxonida", "Spirosclerophorida",
							 "Streptosclerophorida", "Amphidiscosa", 
							 "Megalithistida", "Monalithistida","Tetralithistida")] <- "siliceous sponges"

# Echidoderms
reefs[reefs$class %in% c("Crinoidea", "Cystoidea", "Parablastoidea", 
						 "Paracrinoidea", "Eocrinoidea", "Rhombifera", "Diploporita"),]$group <- "pelmatozoans"


# Microbes
reefs$group[reefs$genus %in% microbe$Genus | reefs$genus=="Allonema"] <- "microbes"
reefs$group[reefs$phylum == "Cyanobacteria"] <- "microbes"
reefs$group[reefs$genus %in% c("Shamovella", "Tubiphytes", "Crescentiella")] <- "tubiphytes"

#Bryozoans
reefs$group[reefs$phylum=="Bryozoa"] <- "other bryozoans"
reefs$group[reefs$class == "Stenolaemata" | reefs$identified_name == "Fasciculipora janinae"] <- "reef-building bryozoans"

# Forams
reefs$group[reefs$phylum=="Foraminifera"] <- "other foraminifera"
reefs$group[reefs$family %in% c("Acervulinidae", "Homotrematidae",
								"Planorbulinidae")] <- "encrusting foraminifera"

# Molluscs
reefs$group[reefs$class=="Gastropoda"] <- "other gastropods"
reefs$group[reefs$family=="Vermetidae"] <- "vermetid gastropods"

reefs$group[reefs$class=="Bivalvia"] <- "other bivalves"
reefs$group[reefs$order=="Hippuritida"] <- "rudists"
reefs$group[reefs$order=="Ostreida" | reefs$genus == "Tridacna"] <- "other reef-building bivalves"

# Brachiopods
reefs$group[reefs$phylum=="Brachiopoda"] <- "other brachiopods"
reefs$group[reefs$family == "Richthofeniidae"] <- "richthofeniid brachiopods"

# Worms
reefs$group[reefs$class=="Polychaeta"] <- "worms"

# Algae
reefs$group[reefs$phylum=="Rhodophyta"] <- "red algae"
reefs$group[reefs$phylum=="Chlorophyta"] <- "green algae"

x <-table(reefs$group)
x
sum(x)
table(reefs[is.na(reefs$group),]$phylum)
#View(reefs[is.na(reefs$group),])

## Assign NA to all non-reef builders
not.build <- c("green algae", "keratose sponges", "other bivalves", "other brachiopods",
               "other bryozoans", "other cnidarians", "other foraminifera", "other gastropods")
reefs$group[reefs$group%in%not.build] <- NA

reefs <- reefs[!is.na(reefs$group),]
table(reefs$reefs)

# Store the results ...
save(reefs, file = file.path("data", "reefs.Rdata"))

