library(cooccur)
library(ggplot2)
library(dplyr)
library(patchwork)
library(divDyn)

source("scripts/functions.R")
pal <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]

# Load data ---------------------------------------------------------------
load(file.path("data", "reefs.Rdata"))

reefs <- bin_occ(reefs, "tens")
nrow(reefs)

reefs <- reefs[!is.na(reefs$stg),]
nrow(reefs)

# recategorisations

# Corals
reefs$group[grep("corals", reefs$group)] <- "corals"

#Calcareous sponges
calc_sp <- c("stromatoporoids", "coralline sponges", "sponges", "pharetronid sponges",
             "pharetronids", "spongiomorphids", "archaeocyaths", "chaetetids",
             "lithistid sponges",   "calcisponges" , "calcareous sponges", "sclerosponges", "stromatoporoids s.str.")

reefs$group[reefs$group %in% calc_sp] <- "calcareous sponges"
reefs$group[grep("stromatoporoids", reefs$group)] <- "calcareous sponges"

#microbes
reefs$group[grep("tubiphytes", reefs$group)] <- "microbes"

# other
reefs$group[grep("bivalves", reefs$group)] <- "other bivalves"
reefs$group[grep("bryozoans", reefs$group)] <- "bryozoans"
reefs$group[grep("vermetid", reefs$group)] <- "vermetids"
reefs$group[grep("foram", reefs$group)] <- "foraminifera"

# make grid
reefs <- reefs[!is.na(reefs$paleolat)|!is.na(reefs$paleolng),]

comp <- list()

resolutions <- c(0.01, 0.05, 0.1, 0.5, 1, 5)

for(re in resolutions){
  
  # sz <- findClosest(re, tessguide$meanEdgeLength_km)
  # gr <- hexagrid(sz, sp=TRUE) 
  # reefs$cell <- locate(gr, reefs[, c("paleolng","paleolat")])
  #reefs$cell2 <- paste(reefs$cell, reefs$stg) #quite some differences
  
  gr <- assign_grid_points(reefs$paleolng, reefs$paleolat, cellsize=c(re,re))
  reefs$cell <- as.numeric(factor(paste(gr$x, gr$y)))
  reefs$cell2 <- paste(reefs$cell, reefs$stg) #temporal
  
  cooc <- purrr::map(c("reef", "non-reef"), function(r){
    res <- list()
    
    for (i in 1:95){
      tempr <- subset(reefs, stg==i & reefs==r)
      
      if(nrow(tempr)>0){
        # create site x sp matrix for reefs
        spp <- table(tempr$group, 
                     tempr$cell2 # in each cell per time
        )
        spp[spp>0] <- 1
        
        sppr <- spp[,!colSums(spp)==0] # remove any column that doesn't have any observations
        
        rr <- do_cooccur(sppr, 
                         thresh=F)
        
        if(!is.null(rr)){
          rr$results$stg <- i
          
          res[[i]]<- rr
        }
        
        
        
        
      }
      
    }
    
    return(res)
    
  }
  )
  
  
  
  reef <- Filter(Negate(is.null), cooc[[1]])
  nonreef <- Filter(Negate(is.null), cooc[[2]])
  
  
  nr<- do.call(rbind, 
               lapply(reef, function (x) data.frame(table(assoc_type(x$cooccur)[[1]]$value),
                                                    stg=unique(x$results$stg)))
  )
  
  nnr<- do.call(rbind, 
                lapply(nonreef, function (x) data.frame(table(assoc_type(x$cooccur)[[1]]$value),
                                                        stg=unique(x$results$stg)))
  )
  
  
  nr$cat <- "r"
  nnr$cat <- "nr"
  
  nr <- nr %>% 
  group_by(cat, stg) %>% 
    mutate(Freq=Freq/sum(Freq))
  
  nnr <- nnr %>% 
    group_by(cat, stg) %>% 
    mutate(Freq=Freq/sum(Freq))
  
  tt <- matrix(nrow = 3, ncol=3)
  tt[,1] <- c(0, 1, -1)
  
  for(assoc in 1:3){
    
    temp<-  t.test(nr[nr$Var1==tt[assoc, 1],]$Freq, nnr[nnr$Var1==tt[assoc, 1],]$Freq)
    tt[assoc, 2] <- temp$statistic
    tt[assoc, 3] <- temp$p.value
  }
  
  colnames(tt) <- c("association", "statistic", "p-value")
  
  write.csv(tt, file.path("output", 
                          sprintf("assoc_comparison_%s.csv", gsub("\\.", "_", re))),
            row.names=F)
  
  tot <- rbind(nr, nnr) %>% 
    group_by(cat, Var1) %>% 
    summarise(mean=mean(Freq), ci=qnorm(.975)*(sd(Freq)/sqrt(length(Freq))))
  
  tot$cat <- factor(tot$cat, levels=c("r", "nr"))
  
  comp[[which(resolutions == re)]] <- ggplot(tot, aes(x=Var1, y=mean, ymin=mean-ci, ymax=mean+ci, col=cat)) +
    geom_point(position=position_dodge(0.5)) +
    geom_errorbar(width=0.2, position=position_dodge(0.5)) +
    scale_color_manual(values=pal[c(4,2)],
                       label=c(r="Reef", nr="Non-Reef")) +
    scale_x_discrete(labels=c("0"="random", "1"="positive", "-1"="negative")) +
    #scale_y_continuous(breaks=seq(-2,15, 2)) +
    labs(x="Association", y="Mean proportion of pairs", col="Environment",	
         caption=paste("Size of cell =",
                       re ,"degrees"))+
    ggthemes::theme_hc() +
    theme(axis.title = element_text(face="bold"),
          legend.title=element_text(face="bold"),
          plot.caption=element_text(face="italic")) +
    annotate("text",x=as.character(tt[,1]), y=Inf, 
             label=paste("italic(W)==", format(tt[,2], digits=2)), parse=T,
             vjust=1, size=3) +
    annotate("text",x=as.character(tt[,1]), y=Inf, 
             label=paste("italic(p)==", format(tt[,3], digits=2)), parse=T,
             vjust=2, size=3)
  
}

# ggsave("figs/fig_02_reef_nonreef.svg", comp[[1]], w=4, h=4)

svg("figs/supplement/fig_02_reef_nonreef.svg", w=12, h=8)

dev.off()

p <- wrap_plots(comp, 2,3) +
  plot_annotation(tag_levels="a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect") &
  theme(plot.tag=element_text(size=10),
        legend.position = "bottom")

ggsave("figs/fig_02_reef_nonreef.svg", p, w=12, h=8)
