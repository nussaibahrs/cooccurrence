# Libraries ---------------------------------------------------------------

library(chronosphere)
library(divDyn)
library(cooccur)
library(tidyverse)
library(ggraph)
library(igraph)
library(patchwork)
library(scales)

source(file.path("scripts", "functions.R"))

source("scripts/02a-pared_binned.R")

#vol=without reef tract, vol1=with reef tract, n.thick=thickness
var <- "vol" 

n.thick <- data.frame(thick=tapply(pared[,var], 
                                   pared$reef_bin, mean, na.rm=T)) %>% 
  tibble::rownames_to_column(var="bin") %>% 
  mutate(bin=as.numeric(bin))

#palette
pal <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

res <- list()

# across time
for(i in 2:50){
	temp <- subset(reefs, bin==i)
	
	# create site x sp matrix
	spp <- table(temp$group, temp$reefs)
	spp[spp>0] <- 1
	
	spp <- spp[,!colSums(spp)==0] # remove any column that doesn't have any observations
	
	res[[i]]  <- tryCatch(do_cooccur(spp),
						  error=function(e)return(NULL))
	
	
}


reef <- lapply(res, function(x) x$results)
df <- do.call(rbind,add_stg(reef))
df <- merge(df, int_summary, by.x="stg", by.y="ten")

# Categorisation

# number of groups according to var
cutoff <- quantile(n.thick$thick, na.rm=T, probs=c(0.25, 0.75))

df.reefs <- reefs %>% 
	group_by(bin, reefs) %>% 
	tally() %>% 
	group_by(bin) %>% 
	summarise(n=median(n)) %>% 
	ungroup() %>% 
	left_join(n.thick) %>% 
	mutate(cat = case_when(thick > cutoff[2] ~ "boom",
						   thick < cutoff[1] ~ "bust",
						   TRUE ~ "normal"))


df <- merge(df, n.thick, by.x="stg", by.y="bin")

mainreef <- c("corals", "sponge", "rudist", "microbe")
df.main <- df[grep(paste(mainreef, collapse="|"), df$assoc),]

n.thick$cat <- "background"
n.thick$cat[n.thick$thick >= cutoff[2]] <- "boom"
n.thick$cat[n.thick$thick <= cutoff[1]] <- "depression"

bins_boom <- n.thick$bin[n.thick$cat=="boom"]

table(n.thick$cat)
sum(table(n.thick$cat))

n.thick <- merge(n.thick, int_summary, by.x="bin", by.y="ten")

p2 <- ggplot(n.thick, aes(x=mid, y=thick)) +
	geom_hline(yintercept = cutoff, col="grey50", linetype="dashed") +
	geom_line(col="grey60", size=1) +
	geom_point(aes(col=cat), size=2.5, shape=21, fill="white", stroke=1) +
	xlim(541,0) +
	scale_color_manual(values=pal[c(1,4,5)], 
					   breaks = c("depression", "boom", "background"))+
	scale_shape_manual(values=c(18, 17, 16), breaks = c("depression", "boom", "background")) +
	coord_cartesian(clip="off")+
	labs(x="Age (Ma)", y="Reef volume (m3)", 
		 col="Reef interval", shape="Reef interval") +
	annotate("rect", xmin=540, xmax=0, ymin=0, ymax=cutoff[1], alpha=0.2)+
	annotate("text", x=0, y=(0+cutoff[1])/2, label="Reef\ndepression", col=pal[1], fontface=2, hjust=1, cex=3) +
	annotate("text", x=0, y=(cutoff[2]+cutoff[1])/2, label="Background\nintervals", col=pal[5], fontface=2, hjust=1, cex=3) +
	annotate("text", x=0, y=(cutoff[2]+max(n.thick$thick))/2, label="Reef\nbooms", col=pal[4], fontface=2, hjust=1, cex=3) +
	ggthemes::theme_hc() +
	theme(legend.position = "none",
		  axis.title=element_text(face=2)) +
	scale_y_continuous(trans="log10",
					   labels = trans_format("log10", math_format(10^.x)))


p3 <- ggplot(n.thick, aes(y=thick))+
	geom_boxplot(col="grey60") +
	scale_y_continuous(trans="log10",
					   labels = trans_format("log10", math_format(10^.x))) +
	ggthemes::theme_hc() +
	theme(axis.text=element_blank(),
		  axis.title = element_blank()
	)

svg("figs/fig_01_reef_vol.svg", w=8, h=4)
cowplot::ggdraw(deeptime::gggeo_scale(p2, size=3, height = unit(1, "lines"))) +
	p3 +
	plot_layout(width=c(1,0.1))
dev.off()

# check cooccurrences in bin with high thickness
boom <- n.thick$bin[n.thick$thick > cutoff[2]]
bust <- n.thick$bin[n.thick$thick < cutoff[1]]
length(boom)
length(bust)

# network
df.net <- df[df$value == 1 & !df$stg %in% c(boom, bust),] #only sig occurrences
df.boom2 <- df[df$stg %in% boom & df$value==1, ]

clust.algo <- cluster_optimal

# Boom ----
edges <- df.boom2 %>% 
	dplyr::select(to=sp1_name.x, from=sp2_name.x, obs=obs_cooccur) %>% 
	group_by(to, from) %>% 
	summarise(weight=sum(obs)) 

n.main <- reefs %>% 
	filter(bin %in% boom) %>% 
	group_by(group) %>% 
	tally() 

g1 <- graph_from_data_frame(edges, directed = F)
plot(g1)

V(g1)$n  <- n.main$n[n.main$group %in% V(g1)$name] #adding1 tally

set.seed(42)
wtc1 <- clust.algo(g1)
modularity(wtc1)

E(g1)$cross <- ifelse(crossing(wtc1, g1),1,0)

V(g1)$community <- wtc1$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(g1, vertex.color=colrs[V(g1)$community])

p1 <- ggraph(g1, "with_gem") +
	geom_edge_link(
		aes(edge_width=weight, edge_colour=cross),
		
	) +
	#highlight cross
	geom_edge_link(
		aes(edge_width=weight, alpha=cross),
		colour="lightgrey",
	) +
	
	scale_edge_width(range=c(0.5,1)) +
	geom_node_point(aes(col=factor(community), size=n), show.legend = FALSE) +
	scale_size(range=c(3,11), limits=c(0,1800))+
	geom_node_text(aes(label=name), size=3, vjust=-2, hjust=0.5, fontface="bold") +
	coord_cartesian(clip="off")+
	theme_graph() +
	scale_color_manual(values=pal[c(4,3, 1)]) +
	theme(legend.position = "none") 
p1

# Background ----
edges <- df.net %>% 
	#filter(value ==1) %>% 
	dplyr::select(to=sp1_name.x, from=sp2_name.x, obs=obs_cooccur) %>% 
	group_by(to, from) %>% 
	summarise(weight=sum(obs))

g3 <- graph_from_data_frame(edges, directed = F)

V(g3)$n  <- n.main$n[n.main$group %in% V(g3)$name] #adding1 tally

set.seed(42)
wtc3 <- clust.algo(g3)
modularity(wtc3)

V(g3)$community <- wtc3$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(g3, vertex.color=colrs[V(g3)$community])

E(g3)$cross <- ifelse(crossing(wtc3, g3),1,0)

p3 <- ggraph(g3, "with_gem") +
	geom_edge_link(
		aes(edge_width=weight, edge_colour=cross),
	) +
	#highlight cross
	geom_edge_link(
		aes(edge_width=weight, alpha=cross),
		colour="lightgrey",
	) +
	
	scale_edge_width(range=c(0.5,1)) +
	geom_node_point(aes(col=factor(community), size=n), show.legend = FALSE) +
	scale_size(range=c(3,11), limits=c(0,1800))+
	geom_node_text(aes(label=name), size=3, vjust=-2, hjust=0.5, fontface="bold") +
	coord_cartesian(clip="off")+
	theme_graph() +
	scale_color_manual(values=pal[c(3,4)]) +
	theme(legend.position = "none")
p3


# Post permian -----------------------------------------------------------

# * Boom ----
edges <- df.boom2 %>% 
  filter(stg > 25)  %>% # 25 == Permian 4
  dplyr::select(to=sp1_name.x, from=sp2_name.x, obs=obs_cooccur) %>% 
  group_by(to, from) %>% 
  summarise(weight=sum(obs)) 

n.main <- reefs %>% 
  filter(bin %in% boom) %>% 
  group_by(group) %>% 
  tally() 

g4 <- graph_from_data_frame(edges, directed = F)
plot(g4)

V(g4)$n  <- n.main$n[n.main$group %in% V(g4)$name] #adding1 tally

set.seed(42)
wtc4 <- clust.algo(g4)
modularity(wtc4)

E(g4)$cross <- ifelse(crossing(wtc4, g4),1,0)

V(g4)$community <- wtc4$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(g4, vertex.color=colrs[V(g4)$community])

p4 <- ggraph(g4, "with_gem") +
  geom_edge_link(
    aes(edge_width=weight, edge_colour=cross),
    
  ) +
  #highlight cross
  geom_edge_link(
    aes(edge_width=weight, alpha=cross),
    colour="lightgrey",
  ) +
  
  scale_edge_width(range=c(0.5,1)) +
  geom_node_point(aes(col=factor(community), size=n), show.legend = FALSE) +
  scale_size(range=c(3,11), limits=c(0,1800))+
  geom_node_text(aes(label=name), size=3, vjust=-2, hjust=0.5, fontface="bold") +
  coord_cartesian(clip="off")+
  theme_graph() +
  scale_color_manual(values=pal[c(4,3, 1)]) +
  theme(legend.position = "none") 

p4

# * Boom ----
edges <- df.boom2 %>% 
  filter(stg <= 25)  %>% # 25 == Permian 4
  dplyr::select(to=sp1_name.x, from=sp2_name.x, obs=obs_cooccur) %>% 
  group_by(to, from) %>% 
  summarise(weight=sum(obs)) 

n.main <- reefs %>% 
  filter(bin %in% boom) %>% 
  group_by(group) %>% 
  tally() 

g5 <- graph_from_data_frame(edges, directed = F)
plot(g5)

V(g5)$n  <- n.main$n[n.main$group %in% V(g5)$name] #adding1 tally

set.seed(42)
wtc5 <- clust.algo(g5)
modularity(wtc5)

E(g5)$cross <- ifelse(crossing(wtc5, g5),1,0)

V(g5)$community <- wtc5$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(g5, vertex.color=colrs[V(g5)$community])

p5 <- ggraph(g5, "with_gem") +
  geom_edge_link(
    aes(edge_width=weight, edge_colour=cross),
    
  ) +
  #highlight cross
  geom_edge_link(
    aes(edge_width=weight, alpha=cross),
    colour="lightgrey",
  ) +
  
  scale_edge_width(range=c(0.5,1)) +
  geom_node_point(aes(col=factor(community), size=n), show.legend = FALSE) +
  scale_size(range=c(3,11), limits=c(0,1800))+
  geom_node_text(aes(label=name), size=3, vjust=-2, hjust=0.5, fontface="bold") +
  coord_cartesian(clip="off")+
  theme_graph() +
  scale_color_manual(values=pal[c(4,3, 1)]) +
  theme(legend.position = "none") 

p5



# Summarising -------------------------------------------------------------
mat <- c(
	modularity(wtc1),
	modularity(wtc3),
	modularity(wtc4),
	modularity(wtc5))

names(mat) <- c("boom", "normal", "post-PT", "pre-PT")

# rownames(mat) <- c("post-T", "pre-T") 
mat
mat[1]/mat[2]
mat[3]/mat[4]

# size
rbind(c(vcount(g1), gsize(g1)),
	  c(vcount(g3), gsize(g3))
)


# different from random
it=2000

r <- lapply(list(g1, g3, g4, g5), random.mod, Nperm=it)
r <- do.call(rbind, r)
r


# add annotation
p_boom <- p1 + 
	annotate("text", x=Inf, y=-Inf, label=paste("m=", signif(modularity(wtc1), 3), ", p=", signif(r[1,2],3)), 
			 fontface=3, hjust=1)
p_boom

p_background <- p3 + 
	annotate("text", x=Inf, y=-Inf, label=paste("m=", signif(modularity(wtc3), 3), ", p=", signif(r[2,2],3)), 
			 fontface=3, hjust=1)

p_background

p_postPT <- p4 + 
  annotate("text", x=Inf, y=-Inf, label=paste("m=", signif(modularity(wtc4), 3), ", p=", signif(r[3,2],3)), 
           fontface=3, hjust=1)#

p_prePT <- p5 + 
  annotate("text", x=Inf, y=-Inf, label=paste("m=", signif(modularity(wtc5), 3), ", p=", signif(r[4,2],3)), 
           fontface=3, hjust=1)


svg(paste0("figs/fig_03_network_", var, ".svg"), w=8, h=8)
p_background + p_boom +
  p_prePT + p_postPT +
	plot_layout(nrow=2) +
	plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

dev.off()

# Per time ----------------------------------------------------------------
head(df)

bns <- unique(df$stg)

res <- list()
it <- 500

for(i in 1:length(bns)){ 
	
	cat("\r", i, " out of ", length(bns))
	
	edges <- df %>% 
		filter(stg == bns[i] & value==1) %>% 
		dplyr::select(to="sp1_name.x", from="sp2_name.x", weight=obs_cooccur ) %>% 
		group_by(to, from) %>% 
		summarise(weight=sum(weight))
	
	if(nrow(edges) > 0){
		g3 <- graph_from_data_frame(edges, directed = F)
		
		V(g3)$n  <- n.main$n[n.main$group %in% V(g3)$name] #adding1 tally
		
		set.seed(42)
		wtc3 <- clust.algo(g3)
		modularity(wtc3)
		
		r <- random.mod(g3, Nperm=it, p=0.5, clust.algo=clust.algo, plot=F)
		
		res[[i]] <- list(stg = bns[i],
						 edges=edges,
						 graph=g3,
						 modularity=modularity(wtc3),
						 random = r)
	}
	
}


res <- Filter(Negate(is.null), res) # remove empties

modt <- data.frame(stg=sapply(res, function(x) x$stg),
				   modularity=sapply(res, function(x) x$modularity),
				   p = sapply(res, function(x){
				   	r <- x$random
				   	ifelse(r[1] > 0.5, r[2], r[1])
				   	
				   }),
				   n = sapply(res, function(x) nrow(x$edges))
)

modt$boom <- 0
modt$boom[modt$stg %in% bins_boom] <- 1

modt <- merge(modt, int_summary, by.x="stg", by.y="ten")

p5 <- ggplot(modt, aes(x=mid, y=modularity)) +
	geom_vline(xintercept=252, col="grey80", linetype="dashed") +
	geom_line(col="grey60", size=1, linetype="dashed") +
	geom_point(aes(col=as.factor(boom), shape=as.factor(boom)),
			   size=2.5, stroke=1, fill="white") +
	scale_x_continuous(trans="reverse") +
	ggthemes::theme_hc() +
	scale_color_manual(values=pal[c(3,4)], labels=c("0"="no","1"="yes")) +
	scale_shape_manual(values=c(21, 16), labels=c("0"="no","1"="yes")) +
	labs(x="Age (Ma)", y="Modularity", shape="Reef boom", 
		 color="Reef boom") +
	theme(axis.title=element_text(face="bold"),
		  legend.title=element_text(face="bold"),
		  legend.position = "none")

p5 <- cowplot::ggdraw(deeptime::gggeo_scale(p5, size=3, height = unit(1, "lines")))
p5

ggsave("figs/fig04_modularity_time_1.svg", p5, w=6, h=4)

# reproduce evenness vs time

modt$era <- "Paleozoic intervals"
modt$era[modt$bottom < 250] <- "Mesozoic & Cenozoic intervals"

modt$era <- factor(modt$era, levels=c("Paleozoic intervals", "Mesozoic & Cenozoic intervals"))

# boom vs background
wilcox.test(
  modt %>% filter(boom == 1) %>% pull(modularity),
  modt %>% filter(boom == 0) %>% pull(modularity),
  alternative = "two.sided", , exact=FALSE
)

# booms vs background in the Mesozoic and Cenozoic
wilcox.test(
  modt %>% filter(boom == 1 & era != "Paleozoic intervals") %>% pull(modularity),
  modt %>% filter(boom == 0 & era != "Paleozoic intervals") %>% pull(modularity), 
  alternative = "two.sided", , exact=FALSE
)

# Palezoic vs Mesozoic and Cenozoic
wilcox.test(
  modt %>% filter(era == "Paleozoic intervals") %>% pull(modularity),
  modt %>% filter(era != "Paleozoic intervals") %>% pull(modularity),
  lternative = "two.sided", , exact=FALSE
)


# add volume
modt <- modt %>%  left_join(n.thick %>%  select(stg=bin, vol=thick))

p6 <- ggplot(modt, aes(x=vol, y=modularity, col=era)) +
  geom_point(aes(shape=as.factor(boom)), size=2.5) +
  scale_shape_manual(values=c(21, 16), labels=c("0"="no","1"="yes")) +
  scale_color_manual(values=pal[c(2,5)]) +
  scale_x_log10(labels = label_log()) +
  labs(x="Reef volume (m3, logged-scale)", y="Modularity", shape="Reef boom", 
       color="Era") +
  theme(axis.title=element_text(face="bold"),
        legend.title=element_text(face="bold"))

ggsave("figs/fig04_modularity_time_2.svg", p6, w=6, h=4)  
