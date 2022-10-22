library(tidyverse)
library(ggstream)
library(cooccur)

source("scripts/functions.R")

pal <- c("#5f0f40", "#9a031e", "#fb8b24", "#e36414", "#0f4c5c")

load("data/reefs_group.RData")

reefs$group[reefs$group %in% c("rudists", "other bivalves")] <- "bivalves"
reefs$group[grep("algae", reefs$group)] <- "algae"
reefs$group[grep("sponge", reefs$group)] <- "sponges"

imp <- c("corals", "sponges", "algae", "microbes", "bivalves")
reefs$group[!reefs$group %in% imp] <- "other"
reefs$group <- factor(reefs$group, levels=c(imp, "other"))

intervals <- read.csv(file.path("data", "pareddat", "l_intervall_new.csv"))

int_summary <- data.frame(bin=1:50,
						  min_ma=tapply(intervals$min.ma, intervals$reef_bin, min),
						  max_ma=tapply(intervals$max.ma, intervals$reef_bin, max))
int_summary$mid <- (int_summary$max_ma + int_summary$min_ma)/2

nreefs <- reefs %>% group_by(bin, group) %>% 
	tally() %>% 
	left_join(int_summary)

nn <- 500

p <- ggplot(nreefs, aes(mid, n, fill = group, label=stringr::str_to_sentence(group))) +
	geom_stream(extra_span = 0.02, type="mirror", n_grid = nn) +
	geom_stream_label(size = 4, type="mirror", n_grid = nn, col="white") +
	scale_fill_manual(values=c(pal, "darkgrey")) +
	scale_x_continuous(trans="reverse", expand=expansion(0),
					   limits=c(550,0)) +
	labs(x = "Age (ma)",
		 y = "Number of reef sites") +
	cowplot::theme_minimal_vgrid(font_size = 14) +
	theme(legend.position ="none",
		  axis.title = element_text(face=2))

x1 <- c(550, 541, 252, 66,0)
p <- p + annotate("rect", xmin=x1[1:4], xmax=x1[2:5], ymin=-200, ymax=-230,
				  col="black", fill="white") +
	annotate("text", x=(x1[1:4] + x1[2:5])/2, y=-215, 
			 label=c("", "Paleozoic", "Mesozoic", "Cenozoic"))
p
ggsave("figs/nreefs.svg", w=8, h=5)

nreefs2 <- cbind(int_summary, 
				 nr= as.numeric(table(reefs$bin[duplicated(reefs$reefs)==FALSE])))
nreefs2$x <- "a"

p <- ggplot(nreefs2, aes(x=mid, y=nr, fill=x)) +
	geom_stream(type="ridge", bw=0.4)+
	scale_x_continuous(trans="reverse", expand=expansion(0),
					   limits = c(550,0)) +
	labs(x="Age (ma)", y="Number of reef sites") +
	scale_fill_manual(values="#EBB55A") +
	theme_minimal(base_size = 14) +
	theme(axis.text = element_text(face=2),
		  legend.position = "none")

p


t <- nreefs2 %>% 
	slice_max(order_by=nr, n=6)
t <- t[-3,]
t$nr[4] <- 140
t$nr[2] <- 225
t$nr[1] <- 390


p <- p + annotate("rect", xmin=x1[1:4], xmax=x1[2:5], ymin=-25, ymax=0,
				col="black", fill="white") +
	annotate("text", x=(x1[1:4] + x1[2:5])/2, y=-12.5, 
			 label=c("", "Paleozoic", "Mesozoic", "Cenozoic"))
p
ggsave("figs/nreefs.png", p, w=8, h=5)

p + annotate("segment", x=t$mid, y=t$nr, yend=t$nr+20, xend=t$mid, 
			 arrow=arrow(ends="first", length = unit(0.05, "inches"), 
			 			type="closed"), col=pal[1]) 

ggsave("figs/nreefs2_1.png", w=8, h=5)

x2 <- c(355, 250, 201, 174, 56)
p + geom_vline(xintercept = x2, lty=2)

ggsave("figs/nreefs2_crises.svg", w=8, h=5)


#legend
p <- ggplot(nreefs[nreefs$group %in% c("corals", "bivalves"),], 
			aes(mid, n, fill = group, label=stringr::str_to_sentence(group))) +
	geom_stream(extra_span = 0.5, type="mirror", n_grid = nn) +
	scale_fill_manual(values=c("darkgrey", "lightgrey")) +
	scale_x_continuous(trans="reverse", expand=expansion(0),
					   limits=c(200, -100)) +
	theme_void() +
	theme(legend.position = "none")

p + annotate("text", x=0, y=0, 
			 label="Number of reefs in\nspecific reef-building group",
			 hjust=0)

ggsave("figs/nreefs_legend.svg", w=6, h=3)

# reefs groups
load("data/reefs_group.RData")

reefs$group[grep("algae", reefs$group)] <- "algae"

reefs_builders <- reefs %>%  distinct(reefs, group) %>% 
	group_by(reefs) %>% 
	tally() %>% 
	ungroup()
	

reefs_builders$n[reefs_builders$n > 3] <- "3+"	

reefs_summary <- reefs_builders %>% 
	group_by(n) %>% 
	summarise(nn = n()) %>%
	mutate(freq = nn / sum(nn))

reefs_summary$fll <- 0
reefs_summary$fll[reefs_summary$n==2] <- 1

reefs_summary$x <- 1:nrow(reefs_summary)

ggplot(reefs_summary, aes(x=x, y=nn, fill=factor(fll))) +
	geom_bar(stat="identity", width = 0.6) +
	scale_fill_manual(values=pal[c(1,3)], guide="none") +
	labs(x="Number of reef-builders per reef", 
		 y="Number of reef sites") +
	ggthemes::theme_hc() +
	geom_text(aes(label=nn), hjust=1.2, col="white", fontface="bold", size=3.5) +
	scale_x_continuous(trans="reverse", breaks=1:4, labels = reefs_summary$n) +
	theme(axis.title=element_text(face="bold")) +
	coord_flip()

ggsave("figs/reef_builders.png", w=6, h=4)

# two categories
nn <- reefs_builders$reefs
xx <- reefs[reefs$reefs %in% nn,]

spp <- table(xx$group, xx$reefs)
spp[spp>0] <- 1
spp <- spp[,!colSums(spp)==0] # remove any column that doesn't have any observations

xx_cc <- do_cooccur(spp)$results
xx_cc <- xx_cc[order(xx_cc$obs_cooccur, decreasing = T),]
xx_cc$prop <- xx_cc$obs_cooccur/sum(xx_cc$obs_cooccur)

assoc <- xx_cc[1:5,c("assoc", "prop")]
assoc <- rbind(assoc, 
			   cbind(assoc="other", prop=1-sum(assoc$prop)))
assoc$prop <- as.numeric(assoc$prop)
assoc$assoc <- factor(assoc$assoc, levels = rev(assoc$assoc))

ggplot(assoc, aes(x = 2, y = prop, fill = assoc)) +
	geom_col(position = "stack", width = 1) +
	xlim(0.5, 2.5) +
	coord_polar("y", start=-0) +
	scale_fill_manual(values=c("lightgrey", rev(pal))) +
	theme_void() +
	theme(legend.position="bottom")

ggsave("figs/bipartite_reefs,.svg", w=6, h=6)
