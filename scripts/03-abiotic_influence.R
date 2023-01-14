library(rgdal)
library(rgeos)
library(raster)
library(forecast)
library(corrplot)
library(gtsummary)
library(broom)

# Load data ---------------------------------------------------------------
source("scripts/02a-pared_binned.R")

#vol=without reef tract, vol1=with reef tract, n.thick=thickness
var <- "vol" 

n.thick <- data.frame(thick=tapply(pared[,var], 
                                   pared$reef_bin, mean, na.rm=T)) %>% 
  tibble::rownames_to_column(var="bin") %>% 
  mutate(bin=as.numeric(bin))


#palette
u_col <- ggsci::pal_uchicago("default")(9)[c(2,5,4,3,1)]


# Continental shelf

paleomap <- fetch("paleomap", "paleocoastlines", datadir = "data/chronodat")

# Tropical only, between 35 degrees N & S
available <- rep(NA, nrow(paleomap))

for(i in 1:nrow(paleomap)){
  t1 <- raster::crop(paleomap[rownames(paleomap)[i], "margin"], 
                     extent(c(-180, 180, -35, 35)))
  t2 <- raster::crop(paleomap[rownames(paleomap)[i], "coast"], 
                     extent(c(-180, 180, -35, 35)))
  
  margin <- raster::area(t1) / 1000000
  coast <- raster::area(t2) / 1000000
  
  available[i] <- sum(margin) - sum(coast)
}

names(available) <- rownames(paleomap)

n.thick$shelf <- available[matchtime(rownames(paleomap), int_summary$mid[n.thick$bin])]

# Temperature -------------------------------------------------------------
temperature <- read.csv("data/Song_2019_low_lat.csv") #Uses GTS2012 according to Xu Dai

modt <- mgcv::gam(temp~s(age), data=temperature, method="REML")
temperature <- data.frame(age=541:0)
temperature$temp <- predict(modt, newdata=temperature)

tms2012 <- gts("2012", "tens")

temperature$bin <- NA

for(i in 1:nrow(tms2012)){
  n <- which(temperature$age <= tms2012$bottom[i] & temperature$age > tms2012$top[i])
  temperature$bin[n] <-tms2012$ten[i]
}

temperature <- temperature[temperature$bin %in% n.thick$bin,]

n.thick$temperature <- NA

n.thick$temperature[n.thick$bin %in% sort(unique(temperature$bin))] <- tapply(temperature$temp, temperature$bin, mean)

# Isotopes ----------------------------------------------------------------
iso <- read.csv("data/prokoph2013.csv") # uses GTS2004 according to paper

tms2004 <- gts("2004", "tens")
iso$bin <- NA

for(i in 1:nrow(tms2004)){
  n <- which(iso$age <= tms2004$bottom[i] & iso$age > tms2004$top[i])
  iso$bin[n] <-tms2004$ten[i]
}

iso$bin[iso$age==0] <- 95

strontium <- data.frame(Sr=tapply(iso$Sr, iso$bin, mean))

strontium$bin <- row.names(strontium)

n.thick <- merge(n.thick, strontium)

sulphur <- data.frame(sulphur=tapply(iso$S, iso$bin, mean))
sulphur$bin <- row.names(sulphur)

n.thick <- merge(n.thick, sulphur)

carbon <- read.csv("data/prokoph2008.csv") # Also uses #GTS2004
colnames(carbon)[c(1,3)] <- c("age", "deltaC") 

modc <- mgcv::gam(deltaC~s(age), data=carbon, method="REML")
carbon <- data.frame(age=541:0)
carbon$deltaC <- predict(modc, newdata=carbon)

carbon$bin <- NA

for(i in 1:nrow(tms2004)){
  n <- which(carbon$age <= tms2004$bottom[i] & carbon$age > tms2004$top[i])
  carbon$bin[n] <-tms2004$ten[i]
}

carbon <- carbon[carbon$bin %in% n.thick$bin,]

n.thick$deltaC <- NA

n.thick$deltaC[n.thick$bin %in% sort(unique(carbon$bin))] <- tapply(carbon$deltaC, carbon$bin, mean)

# Combined ----------------------------------------------------------------
vars <- c("shelf", "temperature", "Sr", "sulphur", "deltaC")
labels <- c(shelf="Available Tropical Shelf Area",
            temperature="Tropical Sea Surface Temperature",
            Sr="<sup>87Sr</sup>/<sup>86</sup>Sr",
            sulphur="&delta;34S",
            deltaC="&delta;13C")

# Autocorrelation
ac <- sapply(vars, function(x){
  ac <- acf(n.thick[,x], plot=F)
  ac$acf[2:6]
})

series <- n.thick[,vars[1]]
significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(series)))

ac <- setNames(data.frame(t(ac)), paste("lag", 1:5)) %>% 
  tibble::rownames_to_column("variable")

gtsummary::theme_gtsummary_compact()

builder <- function(x, Limit){cells_body(columns = !!sym(x), rows = !!sym(x) > Limit)}
ac$Limit <- significance_level

tab2 <- ac %>% 
  gt %>% 
  fmt_number(columns=2:6, decimals=3) %>% 
  tab_style(style = list(cell_fill(color = "grey80"), cell_text(style = "italic")),
            locations = lapply(paste("lag", 1:5), builder, Limit = sym(Limit))) %>% 
  tab_style(style=gt::cell_text(weight = "bold"),
            locations= cells_column_labels()) %>% 
  cols_hide(Limit)

tab2 

gtsave(tab2, "output/table_s_acf.rtf")
file.show("output/table_s_acf.rtf")

# Collinearity -----------------------------------------------------------
labs <- c("Reef thickness", "Shelf area", "Temperature",
          expression("87Sr/86Sr"), "δ34S", "δ13C")

vars.resid <- purrr::map(c("thick", vars), function (x) as.numeric(residuals(auto.arima(n.thick[,x])))) # deal with correlation
vars.resid <- setNames(as.data.frame(vars.resid), c("thick", vars))

# Lag 0
vars.resid <- na.omit(vars.resid)
M = cor(vars.resid, method="spearman")

M.p <- corrplot::cor.mtest(vars.resid, conf.level = 0.95, method="spearman")
colnames(M) <- row.names(M) <- colnames(M.p$p) <- row.names(M.p$p) <- labs

# Lag 1
vars.lag0 <- vars.resid # save for later
vars.resid$thick <- c(NA, vars.resid$thick[-nrow(vars.resid)]) # lag1

vars.resid <- na.omit(vars.resid)
M2 = cor(vars.resid, method="spearman")

M.p2 <- corrplot::cor.mtest(vars.resid, conf.level = 0.95, method="spearman")
colnames(M2) <- row.names(M) <- colnames(M.p2$p) <- row.names(M.p2$p) <- labs

#svg("figs/supplement/fig_s05_correlation.svg", w=8, h=4)
par(mfrow=c(1,2))

corrplot::corrplot(M, p.mat = M.p$p, method = 'circle', type = 'lower', insig='pch', 
                   addCoef.col ='black', number.cex = 0.8, diag=FALSE, tl.col="grey40",
                   tl.cex = 0.8, cl.pos = "n")
mtext("(a)",
      side = 3, adj = 0, line = 0, cex=0.8)

corrplot::corrplot(M2, p.mat = M.p2$p, method = 'circle', type = 'lower', insig='pch',
                   addCoef.col ='black', number.cex = 0.8, diag=FALSE, tl.col="grey40",
                   tl.cex = 0.8, cl.pos = "n")

mtext("(b)",
      side = 3, adj = 0, line = 0, cex = 0.8)

#dev.off()

# Linear regression ----------------------------------------------------------

max <- length(vars)

res <- list()

#linear regression with raw data
mod <- lm(thick~., data=n.thick[,-1])
summary(mod)

final <- step(mod)
summary(final)

#linear regression with autocorrelation removed
mod1 <- lm(thick~., data=vars.resid)
summary(mod1)

final1 <- step(mod1, 
               keep = function(model, aic){
                 list(model=model, 
                      terms = paste(names(model$coefficients)[-1], collapse="+"), 
                      aic = aic)}
               )
summary(final1)

summary(final1)$adj.r.squared
summary(final1)$r.squared

str(summary(final1))

final2 <- lm(thick~temperature +Sr, data=vars.resid)
summary(final2)

modsumm <- broom::tidy(final2) %>% 
  add_row(term="R2", estimate=summary(final1)$adj.r.squared, p.value=anova(final1)$`Pr(>F)`[1]) 

# %>% 
#     bind_rows(
#       broom::tidy(final) %>% 
#         add_row(term="R2", estimate=summary(final)$adj.r.squared, p.value=anova(final)$`Pr(>F)`[1])
#     )


# save results
gtsummary::theme_gtsummary_compact()

tab3 <- modsumm %>% 
  gt() %>% 
  fmt_number(
    columns = 2:3,
    decimals = 3,
    use_seps = TRUE
  ) %>% 
  fmt_number(
    columns = 4:5,
    decimals = 2,
    use_seps = FALSE
  ) %>% 
  cols_label(term=md("**Term**"),
             estimate=md("**Estimate**"),
             std.error = md("**Standard error**"),
             statistic = md("**Statistic**"),
             p.value= md("**P-value**")) %>% 
  cols_width(
    term ~ px(100)
  )

tab3	

gtsave(tab3, "output/table_s_model_final.rtf")
file.show("output/table_s_model_final.rtf")

tab4 <- t(final1$keep[-1,]) %>% 
  as.data.frame() %>% 
  add_column(r2=sapply(final1$keep[1,], function(x){
    summary(x)$adj.r.squared
  })
  ) %>% 
  gt() %>% 
  fmt_number(columns=r2, decimals=3)

gtsave(tab4, "output/table_s_model_final_aic.rtf")
file.show("output/table_s_model_final_aic.rtf")


