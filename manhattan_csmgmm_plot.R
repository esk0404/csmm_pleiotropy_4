library(data.table)


getwd()

setwd("/Users/ekim14/Downloads/Pleiotropy_Project/Final")
pleio_t1 <- fread("~/Downloads/Pleiotropy_Project/Final/Pleiotropy_K4_ssh_t1.csv")
pleio_t2 <- fread("~/Downloads/Pleiotropy_Project/Final/Pleiotropy_K4_ssh_t2.csv")

library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)

manData_t1 <- pleio_t1 %>%
  select(rsid, lfdrResults, Chr, BP, EA, OA)  %>%
  mutate(
    Chr = as.numeric(Chr),
    BP = as.numeric(BP),
    newLfdr = lfdrResults
  )

manData_t2  <- pleio_t2 %>%
  select(rsid, lfdrResults, Chr, BP, EA, OA)  %>%
  mutate(
    Chr = as.numeric(Chr),
    BP = as.numeric(BP),
    newLfdr = lfdrResults
  )


## for colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## plot manhattan function
plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  # arrange data by chromosome
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))
  
  chrCounts <- plotRes %>%
    group_by(Chr) %>%
    summarise(chrLength = max(BP, na.rm = TRUE)) %>%
    arrange(Chr) 
  
  chrOffsets <- cumsum(chrCounts$chrLength)
  names(chrOffsets) <- chrCounts$Chr
  
  truePos <- numeric(nrow(plotRes))
  
  for (chr_it in seq_along(uniqueChrs)) {
    
    tempChr <- uniqueChrs[chr_it]
    idx <- plotRes$Chr == tempChr
    
    offsetVal <- if (chr_it == 1) 0 else chrOffsets[chr_it - 1]
    
    truePos[idx] <- plotRes$BP[idx] + offsetVal
  }
  
  # x-axis ticks at the end of each chromosome
  xBreaks <- chrOffsets 
  xBreaksLabs <- ifelse(uniqueChrs %% 2 == 0, "", uniqueChrs)
  
  plotDat <- plotRes %>% mutate(truePos = truePos)
  
  returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr))) +
    geom_point() +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    ylim(ylimits) +
    theme_cowplot() 
  
  return(returnPlot)
}



manData_t1 <- manData_t1[order(manData_t1$newLfdr), ]
manData_t1 <- manData_t1[!duplicated(manData_t1$rsid), ]
manData_t1<- manData_t1[, c("rsid","Chr","BP","newLfdr","cat")]



manPlot_t1 <- plotManhattan(plotRes = manData_t1, chrCounts,
                            colValues = gg_color_hue(1), shapeValues=c(16), ylimits=c(0, 15), legName="")




manData_t2 <- manData_t2[order(manData_t2$newLfdr), ]
manData_t2 <- manData_t2[!duplicated(manData_t2$rsid), ]
manData_t2<- manData_t2[, c("rsid","Chr","BP","newLfdr")]

manPlot_t2 <- plotManhattan(plotRes = manDataTwo_t2, chrCounts,
                               colValues = gg_color_hue(1), shapeValues=c(16), ylimits=c(0, 15), legName="")





manPlot_t1
manPlot_t2


getwd()
setwd("/Users/ekim14/Downloads/Pleiotropy_Project/Final")

## Save plot
ggsave(
  filename = "Manhattan_Pleio4_t1.png",
  plot = manPlot_t1,
  width = 17,
  height = 7,
  dpi = 300
)


ggsave(
  filename = "Manhattan_Pleio4_t2.png",
  plot = manPlot_t2,
  width = 17,
  height = 7,
  dpi = 300
)
