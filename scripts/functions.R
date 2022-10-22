getReefs <- function(pbdb){
  reefs <- pbdb[grep("Reef", pbdb$collection_aka),]
  
  #remove Old Reef Limestone
  n <- grep("Old Reef", reefs$collection_aka)
  
  reefs <- reefs[-n,]
  reefs$reef_number <- reefs$collection_aka
  
  #triple
  n <- grep("Reefs [0-9]{2,4}, [0-9]+, [0-9]+$", reefs$reef_number)
  reefsn <- reefs$reef_number[n]
  reefsn <- gsub("Reefs ", "", reefsn)
  reefsn <- do.call(rbind, (strsplit(reefsn, ",")))
  
  reefstri <- reefs[n,]
  reefstri <- rbind(reefstri, reefstri, reefstri)
  
  reefstri$reef_number <- paste("Reef ", c(reefsn[,1], reefsn[,2], reefsn[,3]))
  reefs <- rbind(reefs[-n,], reefstri)  
  
  #double
  n <- grep("Reefs [0-9]{2,4}, [0-9]+$", reefs$reef_number)
  reefsn <- reefs$reef_number[n]
  reefsn <- gsub("Reefs ", "", reefsn)
  reefsn <- do.call(rbind, (strsplit(reefsn, ",")))
  
  reefsdou <- reefs[n,]
  reefsdou <- rbind(reefsdou, reefsdou)
  
  reefsdou$reef_number <- paste("Reef ", c(reefsn[,1], reefsn[,2]))
  reefs <- rbind(reefs[-n,], reefsdou)  
  
  #single
  n <- grep("Reef", reefs$reef_number)
  reefs <- reefs[n,]
  
  # Other details
  reefs$reef_number <- gsub(".*(Reef [0-9]+).*", "\\1", reefs$reef_number)
  reefs$reefs <- as.numeric(trimws(gsub("Reef ", "", reefs$reef_number)))
  reefs$reef_number <- NULL
  
  # reefs[is.na(reefs$reefs),]$collection_aka
  
  #with all reef numbers
  reefs <- reefs[!is.na(reefs$reefs),]
  
  return(reefs)
} 

bin_occ <- function (pbdb, level="stage"){
  data(keys, package = "divDyn")
  
  
  ff <- switch(level,
               stage="stgInt", 
               tens = "tenInt")
  
  # the 'stg' entries (lookup)
  stgMin <- divDyn::categorize(pbdb[, "early_interval"], keys[[ff]])
  stgMax <- divDyn::categorize(pbdb[, "late_interval"], keys[[ff]])
  
  
  # convert to numeric
  stgMin <- as.numeric(stgMin)
  stgMax <- as.numeric(stgMax)
  
  # empty container
  pbdb$stg <- rep(NA, nrow(pbdb))
  
  # select entries, where
  stgCondition <- c(which(stgMax == stgMin), which(stgMax == -1))
  
  # in these entries, use the stg indicated by the early_interval
  pbdb$stg[stgCondition] <- stgMin[stgCondition]
  
  # convert to numeric
  stgMin <- as.numeric(stgMin)
  stgMax <- as.numeric(stgMax)
  
  return(pbdb)
}

assoc_type <- function(x){
  dim <- x$species
  comat_pos <- comat_neg <- matrix(nrow=dim,ncol=dim)
  
  co_tab <- x$result
  
  for (i in 1:nrow(co_tab)){
    comat_pos[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_gt"]
    comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_gt"]
    
    row.names(comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]])
    
  }
  
  for (i in 1:nrow(co_tab)){
    comat_neg[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_lt"]
    comat_neg[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_lt"]
  }
  
  comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
  colnames(comat) <- 1:dim
  row.names(comat) <- 1:dim
  
  if ("spp_key" %in% names(x)){
    
    sp1_name <- merge(x=data.frame(order=1:length(colnames(comat)),sp1=colnames(comat)),y=x$spp_key,by.x="sp1",by.y="num",all.x=T)
    sp2_name <- merge(x=data.frame(order=1:length(row.names(comat)),sp2=row.names(comat)),y=x$spp_key,by.x="sp2",by.y="num",all.x=T)
    
    colnames(comat) <- sp1_name[with(sp1_name,order(order)),"spp"]  
    row.names(comat) <- sp2_name[with(sp2_name,order(order)),"spp"]
    
  }  
  
  #ind <- apply(comat, 1, function(x) all(is.na(x)))
  #comat <- comat[!ind,]
  #ind <- apply(comat, 2, function(x) all(is.na(x)))
  #comat <- comat[,!ind]
  
  comat[is.na(comat)] <- 0
  
  # SECTION TO REMOVE SPECIES INTERACTION WITH NO OTHERS
  
  rmrandomspp <- function(orimat,plotrand = FALSE){
    if(plotrand == FALSE){
      ind <- apply(orimat, 1, function(x) all(x==0))
      orimat <- orimat[!ind,]    
      ind <- apply(orimat, 2, function(x) all(x==0))
      orimat <- orimat[,!ind]
    }
    return(orimat)
  }
  
  comat <- rmrandomspp(orimat = comat)
  ####################################################### 
  
  comat <- comat[order(rowSums(comat)),]
  comat <- comat[,order(colSums(comat))]
  
  #comat <- rmrandomspp(orimat = comat, ...)
  
  ind <- apply(comat, 1, function(x) all(x==0))
  comat <- comat[!ind,]
  ind <- apply(comat, 2, function(x) all(x==0))
  comat <- comat[,!ind]
  
  # ind <- apply(comat, 1, function(x) all(x==0))
  # comat <- comat[names(sort(ind)),]
  # ind <- apply(comat, 2, function(x) all(x==0))
  # comat <- comat[,names(sort(ind))]
  
  #comat
  data.m = reshape2::melt(comat)
  colnames(data.m) <- c("X1","X2","value")
  data.m$X1 <- as.character(data.m$X1)
  data.m$X2 <- as.character(data.m$X2)
  
  
  
  dfids <- subset(data.m, X1 == X2)
  
  X1 <- data.m$X1
  X2 <- data.m$X2
  
  df.lower = subset(data.m[lower.tri(comat),],X1 != X2)
  
  return(list(df.lower, dfids))
  
}

### 
plot.cooccur <-
  function(x, title="Species Co-occurrence Matrix", 
           cols= c("#FFCC66","dark gray","light blue"), ...){
    require(ggplot2)
    
    xx <- assoc_type(x)
    df.lower <- xx[[1]]
    df.lower <- df.lower[order(df.lower$X1),]
    meas <- as.character(unique(c(df.lower$X2, df.lower$X1)))
    
    X1 <- df.lower$X1
    X2 <- df.lower$X2
    value <- df.lower$value
    
    p <- ggplot(df.lower, aes(X1, X2)) + 
      geom_tile(aes(fill = factor(value,levels=c(-1,0,1))), 
                colour ="white") 
    p <- p + scale_fill_manual(values = cols, name = "", labels = c("negative","random","positive"),drop=FALSE) + 
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.text=element_text(size=12)
      ) + 
      ggtitle(title) + 
      xlab("") + ylab("") + 
      scale_x_discrete(limits=meas, expand = c(0.3,0),drop=FALSE) + 
      scale_y_discrete(limits=meas, expand = c(0.3,0),drop=FALSE) 
    p + geom_text(data=xx[[2]],aes(label=X1),hjust=1,vjust=0,angle = -22.5)#, color="dark gray")
  }

sort_assoc <- function(sp1_name, sp2_name){
  tmp <- cbind(sp1_name, sp2_name)
  
  for(i in 1:nrow(tmp)){
    tmp[i,] <- sort(tmp[i,])
  }
  
  apply(tmp, 1, function(x) paste(x, collapse=" x "))
}

# do cooccurrence, get effect size
do_cooccur <- function(spp, thresh=T){
  if(length(nrow(spp)) > 0){
    mat <- matrix(spp, nrow=dim(spp)[1])
    rownames(mat) <- row.names(spp)
    colnames(mat) <- colnames(spp)
    
    # run cooccur
    cooccur.reefs <- tryCatch(cooccur::cooccur(mat,
                                               type = "spp_site",
                                               thresh = thresh,
                                               spp_names = TRUE), error=function(e) return (NULL))
    
    if(!is.null(cooccur.reefs)){
      if(cooccur.reefs$co_occurrences !=0){
        #save results
        temp2 <- cooccur::prob.table(cooccur.reefs)
        temp2$assoc <- sort_assoc(temp2$sp1_name, temp2$sp2_name)
        
        assoc_temp <- assoc_type(cooccur.reefs)[[1]]
        colnames(assoc_temp) <- c("sp1_name", "sp2_name", "value")
        
        assoc_temp$assoc <- sort_assoc(assoc_temp$sp1_name, assoc_temp$sp2_name)
        
        eff_temp <- effect.sizes(cooccur.reefs)
        colnames(eff_temp) <- c("sp1_name", "sp2_name", "effects")
        eff_temp$assoc <- sort_assoc(eff_temp$sp1_name, eff_temp$sp2_name)
        
        return(list(cooccur = cooccur.reefs, 
                    results=merge(merge(temp2, assoc_temp, by="assoc"), 
                                  eff_temp, by="assoc")))
      } else return(NULL)
    } else return(NULL)
  }else return(NULL)
}

add_stg <- function(reef){
  nreef <- cbind(stg=1:length(reef), 
                 times=unlist(lapply(reef, function(x){
                   if(is.null(x)) {return(0)} else return(nrow(x))
                 })))
  nreef <- apply(nreef, 1, list)
  
  nreef <- lapply(nreef, function(x) rep(x[[1]]["stg"], x[[1]]["times"]))
  
  mapply(function(x, y) "[<-"(x, "stg", value = y) ,
         reef, nreef, SIMPLIFY = FALSE)
  
}

typ <- c(sort_assoc("corals", "calcareous sponges"),
         sort_assoc("red algae", "corals"),
         sort_assoc("corals", "microbes")
)


# Random network modularity
random.mod <- function(g, Nperm=1000, plot=F, p=0.5, 
                       clust.algo=cluster_spinglass){
  
  randomized.modularity <- lapply(seq_len(Nperm), function(x){  
    
    tryCatch({
      randomnet <- rewire(g, with=each_edge(p)) #rewire vertices with constant probability
      E(randomnet)$weight <- sample(E(g)$weight) #shuffle initial weights and assign them randomly to edges
      
      return(clust.algo(randomnet)$modularity)
    }, 
    error=function(e) return(NA)
    )
    
  })
  
  obs <- clust.algo(g)$modularity
  randomized.modularity <- unlist(randomized.modularity)
  randomized.modularity <- randomized.modularity[!is.na(randomized.modularity)]
  
  if(plot){
    plot(density(randomized.modularity), main="Observed vs. randomized",)
    abline(v=obs, col="red", lwd=2, xlab="")
    
  }
  
  c(mean(randomized.modularity), 
    length(randomized.modularity[randomized.modularity > obs])/length(randomized.modularity))
}

# function to plot pdp
pdp_plot <- function (pdp, olddata=NULL, levels=NULL, unscale=FALSE, col="darkred"){
  lab <- levels(pdp$`_vname_`)
  pdp <- as.data.frame(pdp)
  
  if (unscale == TRUE) {pdp$`_x_` <- unscale(pdp$`_x_`, olddata)}
  
  if(is.factor(pdp$`_x_`)){
    pdp$`_x_` <- factor(as.character(pdp$`_x_`), levels=levels)
    
    p <- ggplot(pdp, aes(x=`_x_`, y=`_yhat_`)) + geom_bar(stat="identity", aes(fill=`_x_`), width = 0.6) +
      labs(x=lab, y="Predicted Value") +
      scale_y_continuous(expand=expand_scale(mult = c(0, .1))) +
      theme_light(base_size = 15) +
      theme(axis.title = element_text(size = 12, face="bold"),
            axis.text.x = element_text(family = "Roboto Mono", size = 10),
            panel.grid = element_blank(),
            legend.position = "none")
  } else {
    r <- range(pdp$`_yhat_`)
    
    p <- ggplot(pdp, aes(x=`_x_`, y=`_yhat_`)) + geom_line(col=col) +
      labs(x=lab, y="Predicted Value") +
      scale_y_continuous(expand=expand_scale(c(0, 0.1)))+
      geom_linerange(aes(x=`_x_`, ymin=r[1], ymax=r[2]/20+r[1]), col="darkgrey") +
      theme_light(base_size = 15) +
      theme(axis.title = element_text(size = 12, face="bold"),
            axis.text = element_text(size = 10),
            panel.grid = element_blank(),
            legend.position = "none")
  }
  
  return(p)
}

gts <- function(timescale, level="stage"){
  valid <-c(2004, 2008, 2012, 2016, 2020)
  
  if(!timescale %in% valid) errorCondition("Please specify a valid timescale")
  if(!level %in% c("stage", "tens")) errorCondition("Please specify a valid level")
  
  timescale <- paste0("gts", timescale)
  
  scl <-read.csv("data/timescale.csv")
  
  
  if(level=="stage"){
    data("stages", package="divDyn")
    
    stages$bottom <- scl[,timescale]
    stages$top <- c(scl[2:95, timescale],0)
    
  } else{
    data("tens", package = "divDyn")
    scl <- scl[4:95,]
    scl <- sort(
      tapply(scl[,timescale], scl$tens, max), decreasing=T)
    stages <- tens
    
    stages$bottom <- scl
    stages$top <- c(scl[-1], 0)
  }
  
  stages$mid <- (stages$bottom+stages$top)/2
  stages$dur <- stages$bottom - stages$top
  
  
  return(stages)
}
# Find stage
find_stg <- function(age, timescale){
  scl <-read.csv("data/timescale.csv")
  
  if(is.null(from)) errorCondition("Please specify a valid timescale")
  
  valid <-c(2004, 2008, 2012, 2016, 2020)
  
  if(!timescale %in% valid) errorCondition("Please specify a valid timescale")
  
  timescale <- paste0("gts", timescale)
  
  # check which bin it corresponds to
  base <- c(scl[,timescale], 0)
  
  bin <- rep(NA, length(age))
  
  for(i in 1:95){
    n <- which(age <= base[i] & age > base[i+1])
    
    if (length(n)>0) bin[n] <- i
  }
  
  bin
}

### assign grid points based on origin and cell size -----

assign_grid_points <- function(x,y, origin=c(0,0), cellsize=c(5,5)) {
  xy <- cbind(x,y)
  temp <- t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
  return(as.data.frame(temp))
}
