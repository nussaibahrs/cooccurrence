library(tabulizer)

doc <- "data/prokoph2013.pdf"
pngfile <- pdftools::pdf_convert(doc, dpi = 600)

# Page 3 ------------------------------------------------------------------
pg1 <- tesseract::ocr(pngfile[3])
pg1 <- do.call(rbind.data.frame, strsplit(
strsplit(pg1, "\n")[[1]][7:77], " "))

View(pg1)
pg1 <- rbind(setNames(pg1[,1:4], 1:4),
             setNames(pg1[,8:11], 1:4))
colnames(pg1) <- c("age", "a10", "Sr", "S")

# Page 4 ------------------------------------------------------------------

pg2 <- tesseract::ocr(pngfile[4])
pg2 <- do.call(rbind.data.frame, strsplit(
  strsplit(pg2, "\n")[[1]][7:76], " "))

pg2 <- rbind(setNames(pg2[,1:4], 1:4),
             setNames(pg2[,8:11], 1:4))
View(pg2)
colnames(pg2) <- c("age", "a10", "Sr", "S")

# Page 5 ------------------------------------------------------------------

pg3 <- tesseract::ocr(pngfile[5])
pg3 <- do.call(rbind.data.frame, strsplit(
  strsplit(pg3, "\n")[[1]][6:76], " "))

View(pg3)
pg3 <- rbind(setNames(pg3[,1:4], 1:4),
             setNames(pg3[,6:9], 1:4))

colnames(pg3) <- c("age", "a10", "Sr", "S")
View(pg3)
pg3[14,"Sr"] <- "0.70793"

# Page 4 ------------------------------------------------------------------

pg4 <- tesseract::ocr(pngfile[6])
pg4 <- strsplit(
  strsplit(pg4, "\n")[[1]][7:77], " ")

pg4 <- lapply(pg4, function(x) setNames(rbind.data.frame(x), 1:length(x)))

pg4 <- as.data.frame(data.table::rbindlist(pg4, fill=T))

pg4 <- rbind(setNames(pg4[,1:4], 1:4),
             setNames(pg4[,6:9],1:4))

colnames(pg4) <- c("age", "a10", "Sr", "S")

dat <- rbind(pg1,pg2,pg3,pg4)
dat <- na.omit(dat)

dat$age[78] <- 77
dat$age <- as.numeric(dat$age)
range(dat$age)

dat[252,3:4] <- c(0.707513, 17.7684)
dat[367,] <- c(366, 0.05199, 0.708397, 20.0302)
dat$Sr <- as.numeric(dat$Sr)
range(dat$Sr)

dat$S <- as.numeric(dat$S)
dat[296,"S"] <- 12.5937
range(dat$S)

write.csv(dat[,c(1,3:4)], "data/prokoph2013.csv", row.names = F)

system("rm *.png")
