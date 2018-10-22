
#' ** Analysis of processed LeafJ data
#' mailto:j.s.hoyer@wustl.edu
#' Started 2016-10-07.
#'
#' The main question of interest for this experiment
#' (and the two before it)
#' is whether transgenes driven by truncated versions
#' versions of the AGO7 promoter
#' are sufficient to largely complement ago7 (zip-1) mutants.
#' Based on the results it is clear
#' that the 298 bp promoter and shorter constructs
#' were not sufficient for complementation
#' whereas longer constructs do suffice.
#' The secondary question is whether
#' complementation with the 453 bp and 422 bp constructs
#' (lacking SPL and TCP binding sites)
#' differs from complementation with longer constructs.
#' Such small differences could indicate a fine-tuning role
#' for the TF binding sites of interest,
#' though more work would be necessary to prove this.
#'
#' 33-dps_leaf-series-scans/readme.txt

#' library(dplyr)
#' library(ggplot2); theme_set(theme_classic())
source("../2016-08-26_zip/30-dps_trichome-and-leaf-6-analysis.R")
#' This in turn runs
#' source("../R/axis-label-strings-and-helper-function.R")
#'
#' I think the following may be equivalent
#' and yield the one data frame we need,
#' but I have not tested this.
incsvfile <- "../2016-08-26_zip/working-data/workedUp923.csv"
#' zworkedUp923 <- read.csv(incsvfile)
#' identical(zworkedUp923, workedUp923) # FALSE
## Factor levels differ.  Possibly numeric precision too.

leafj33 <- read.delim(
    "../2016-08-26_zip/33-dps_leafj-measurements.txt")
leafj35 <- read.delim(
    "../2016-08-26_zip/35-dps_leafj-measurements.txt")

leafj33$dps <- 33
#' dim(leafj33)
leafj35$dps <- 35
#' dim(leafj35)

leafj <- rbind(leafj33, leafj35)
#' Recall that leafJ spits out a reminder of the image file name
#' (which we do not need).
leafj <- select(leafj, -matches("file"))
#' dim(leafj)

#' Scans were done at 600 dpi, and ImageJ automatically recognized this scale.
#' This can be be confirmed by seeing that the 15 cm scale
#' is roughly 15 / 2.54 * 600 pixels wide.
#' Therefore convert LeafJ measurements from inches to centimeters.
str(leafj)

leafj <- mutate(leafj,
                petioleLength  = petioleLength * 2.54,
                bladeLength    = bladeLength   * 2.54,
                bladeWidth     = bladeWidth    * 2.54,
                bladeArea      = bladeArea * 2.54 * 2.54,
                bladePerimeter = bladePerimeter * 2.54)

leafj <- mutate(leafj,
                bpRatio = bladeLength/petioleLength,
                lwRatio = bladeLength/bladeWidth)

str(leafj)

#' Somewhat clumsy replication of the leaf 6 caliper data here:
leafj <- left_join(leafj, workedUp923, by = c("transplantNum" = "plantNum"))
## Reassigning dataframe name here was quite error-prone
## when I was using semi_join!
## Multi-rosette plant 660 not measured with calipers.
filter(leafj, promoter2 == "WT") %>% select(transplantNum)
filter(leafj, promoter2 == "zip-1") %>% select(transplantNum)

#' Could manually set to NA:
#notLeaf6 <- leafj$leaf != 6
#leafj$[notLeaf6] <- NA
#' Maybe filtering is easier.
leaf6j <- filter(leafj, leaf == 6)
str(leafj)

sans645 <- filter(leafj, transplantNum != 645)
sans645_660 <- filter(sans645, transplantNum != 660)
str(sans645_660)

## Axis label strings:
leafStr           <- "Leaf number"
petioleLengthStr  <- "Petiole length (cm)"
bladeLengthStr    <- "Leaf blade length (cm)"
bladeWidthStr     <- "Leaf blade width (cm)"
areaStr           <- "Leaf blade area (cm^2)"
perimeterStr      <- "Leaf blade perimeter (cm)"
circularityStr    <- "Leaf blade circularity"

oneToTen <- scale_x_discrete(limits = 1:10, labels = format(1:10, digits = 0))
oneToEight <- scale_x_discrete(limits = 1:8, labels = format(1:8, digits = 0))

range(leafj$bpRatio)
range(leafj$lwRatio)
filter(leafj, lwRatio > 3)

bpRange <- scale_y_continuous(limits = c(0.7, 3.6))
lwRange <- scale_y_continuous(limits = c(1  , 3.45))

#' ** Colors a hybrid between two palettes:
#' http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
#' http://colorbrewer2.org/#type=sequential&scheme=OrRd&n=9

orangeRed <- c("#fdbb84",
               "#fc8d59",
               "#ef6548",
               "#d7301f",
               "#b30000",
               "#7f0000")

cols <- c(
    "WT" = "#1b9e77",      # green
    "zip-1" = "#7570b3",   # purple
    "1934 bp" = orangeRed[6],
    "495 bp"  = orangeRed[5],
    "453 bp"  = orangeRed[2],
    "422 bp"  = orangeRed[4],
    "298 bp"  = orangeRed[1],
    "0 bp"    = orangeRed[3])

scm <-  scale_colour_manual(values = cols)

first8  <- filter(leafj, leaf < 9)
leafj33b <- filter(leafj, dps == 33)
leafj35b <- filter(leafj, dps == 35)

#' ** Manual jittering and grouping
#' Promoter not annotated for plant 660:
#' tail(leafj, n = 39)
is660 <- leafj$transplantNum == 660
leafj$transgene[is660] <- "pMDC99 GUS"
leafj$promoter[is660] <- "GUS (Col-0)"
leafj$promoter2[is660] <- "WT"
## See 28-to-30-dps_abaxial-trichome-positions-and-true-leaf-6-caliper-measurements.txt
##  and ../2016-05-13_zip/parent-vial-genotypes.txt

set1 <- leafj$promoter2 %in% c("453 bp", "298 bp", "0 bp", "zip-1")
set2 <- leafj$promoter2 %in% c("WT", "1934 bp", "495 bp", "422 bp")
which(set1 == set2)   ## No plants should be in both groups!
leafj$leafjit[set1] <- leafj$leaf[set1] + 0.05
leafj$leafjit[set2] <- leafj$leaf[set2] - 0.05
set.seed(-0.5843456)  # Picked by running rnorm(1)

leafj$leafjit <- leafj$leafjit + runif(479, max = 0.01)

leafj$set <- "a"
leafj$set[leafj$promoter2 %in% c("1934 bp",  "0 bp")] <- "f"
leafj$set[leafj$promoter2 %in% c("495 bp", "453 bp")] <- "g"
leafj$set[leafj$promoter2 %in% c("422 bp", "298 bp")] <- "h"

refLines <- filter(leafj, promoter2 %in% c("WT", "zip-1"))
setA     <- filter(leafj, promoter2 %in% c("WT", "zip-1", codesInOrder[3:5]))
setB     <- filter(leafj, promoter2 %in% c("WT", "zip-1", codesInOrder[6:8]))
setC     <- filter(leafj, promoter2 %in% c("WT", "zip-1", codesInOrder[3:4]))
setD     <- filter(leafj, promoter2 %in% c("WT", "zip-1", codesInOrder[5:6]))
setE     <- filter(leafj, promoter2 %in% c("WT", "zip-1", codesInOrder[7:8]))
setF     <- filter(leafj, promoter2 %in% c("WT", "zip-1", "1934 bp", "0 bp"))
setG     <- filter(leafj, promoter2 %in% c("WT", "zip-1", "495 bp", "453 bp"))
setH     <- filter(leafj, promoter2 %in% c("WT", "zip-1", "422 bp", "298 bp"))

setI     <- filter(leafj, promoter2 %in% c("WT"))
setJ     <- filter(leafj, promoter2 %in% c("WT", "zip-1", "1934 bp"))
setK     <- filter(leafj, promoter2 %in% c("WT", "zip-1",  "495 bp"))
setL     <- filter(leafj, promoter2 %in% c("WT", "zip-1",  "453 bp"))
setM     <- filter(leafj, promoter2 %in% c("WT", "zip-1",  "422 bp"))
setN     <- filter(leafj, promoter2 %in% c("WT", "zip-1",  "298 bp"))
setO     <- filter(leafj, promoter2 %in% c("WT", "zip-1",    "0 bp"))
setP     <- filter(leafj, promoter2 %in% c("zip-1"))

setA2    <- filter(leafj, promoter2 %in% codesInOrder[3:5])
setB2    <- filter(leafj, promoter2 %in% codesInOrder[6:8])
setC2    <- filter(leafj, promoter2 %in% codesInOrder[3:4])
setD2    <- filter(leafj, promoter2 %in% codesInOrder[5:6])
setE2    <- filter(leafj, promoter2 %in% codesInOrder[7:8])
setF2    <- filter(leafj, promoter2 %in% c("1934 bp", "0 bp"))
setG2    <- filter(leafj, promoter2 %in% c("495 bp", "453 bp"))
setH2    <- filter(leafj, promoter2 %in% c("422 bp", "298 bp"))

#' ** Mean calculations for subsets
grouped     <- group_by(leafj,   leaf, promoter2)
grouped2    <- group_by(leafj,   leaf, promoter2, dps)
grouped1to8 <- group_by(first8,leaf, promoter2)
grouped33   <- group_by(leafj33b, leaf, promoter2)
grouped35   <- group_by(leafj35b, leaf, promoter2)
groupedRef  <- group_by(refLines, leaf, promoter2)

groupedA    <- group_by(setA,   leaf, promoter2)
groupedB    <- group_by(setB,   leaf, promoter2)
groupedC    <- group_by(setC,   leaf, promoter2)
groupedD    <- group_by(setD,   leaf, promoter2)
groupedE    <- group_by(setE,   leaf, promoter2)
groupedF    <- group_by(setF,   leaf, promoter2)
groupedG    <- group_by(setG,   leaf, promoter2)
groupedH    <- group_by(setH,   leaf, promoter2)

groupedI    <- group_by(setI,   leaf, promoter2)
groupedJ    <- group_by(setJ,   leaf, promoter2)
groupedK    <- group_by(setK,   leaf, promoter2)
groupedL    <- group_by(setL,   leaf, promoter2)
groupedM    <- group_by(setM,   leaf, promoter2)
groupedN    <- group_by(setN,   leaf, promoter2)
groupedO    <- group_by(setO,   leaf, promoter2)
groupedP    <- group_by(setP,   leaf, promoter2)

groupedA2    <- group_by(setA2,   leaf, promoter2)
groupedB2    <- group_by(setB2,   leaf, promoter2)
groupedC2    <- group_by(setC2,   leaf, promoter2)
groupedD2    <- group_by(setD2,   leaf, promoter2)
groupedE2    <- group_by(setE2,   leaf, promoter2)
groupedF2    <- group_by(setF2,   leaf, promoter2)
groupedG2    <- group_by(setG2,   leaf, promoter2)
groupedH2    <- group_by(setH2,   leaf, promoter2)


summarizeIt <- function(groupedDF)
{
    summarize(groupedDF,
              meanPetioleLength    = mean(petioleLength),
              meanBladeLength      = mean(bladeLength),
              meanBladeWidth       = mean(bladeWidth),
              meanBpRatio          = mean(bpRatio),
              meanLwRatio          = mean(lwRatio),
              meanBladeArea        = mean(bladeArea),
              meanBladePerimeter   = mean(bladePerimeter),
              meanBladeCircularity = mean(bladeCircularity),
              #
              sdPetioleLength    = sd(petioleLength),
              sdBladeLength      = sd(bladeLength),
              sdBladeWidth       = sd(bladeWidth),
              sdBpRatio          = sd(bpRatio),
              sdLwRatio          = sd(lwRatio),
              sdBladeArea        = sd(bladeArea),
              sdBladePerimeter   = sd(bladePerimeter),
              sdBladeCircularity = sd(bladeCircularity),
              #
              semPetioleLength    = sd(petioleLength)/sqrt(6),
              semBladeLength      = sd(bladeLength)/sqrt(6),
              semBladeWidth       = sd(bladeWidth)/sqrt(6),
              semBpRatio          = sd(bpRatio)/sqrt(6),
              semLwRatio          = sd(lwRatio)/sqrt(6),
              semBladeArea        = sd(bladeArea)/sqrt(6),
              semBladePerimeter   = sd(bladePerimeter)/sqrt(6),
              semBladeCircularity = sd(bladeCircularity)/sqrt(6),
              #
              meanPlusSDPetioleLength    = meanPetioleLength    + sdPetioleLength,
              meanPlusSDBladeLength      = meanBladeLength      + sdBladeLength,
              meanPlusSDBladeWidth       = meanBladeWidth       + sdBladeWidth,
              meanPlusSDBpRatio          = meanBpRatio          + sdBpRatio,
              meanPlusSDLwRatio          = meanLwRatio          + sdLwRatio,
              meanPlusSDBladeArea        = meanBladeArea        + sdBladeArea,
              meanPlusSDBladePerimeter   = meanBladePerimeter   + sdBladePerimeter,
              meanPlusSDBladeCircularity = meanBladeCircularity + sdBladeCircularity,
              #
              meanMinusSDPetioleLength    = meanPetioleLength    - sdPetioleLength,
              meanMinusSDBladeLength      = meanBladeLength      - sdBladeLength,
              meanMinusSDBladeWidth       = meanBladeWidth       - sdBladeWidth,
              meanMinusSDBpRatio          = meanBpRatio          - sdBpRatio,
              meanMinusSDLwRatio          = meanLwRatio          - sdLwRatio,
              meanMinusSDBladeArea        = meanBladeArea        - sdBladeArea,
              meanMinusSDBladePerimeter   = meanBladePerimeter   - sdBladePerimeter,
              meanMinusSDBladeCircularity = meanBladeCircularity - sdBladeCircularity,
              #
              meanPlusSEMPetioleLength    = meanPetioleLength    + semPetioleLength,
              meanPlusSEMBladeLength      = meanBladeLength      + semBladeLength,
              meanPlusSEMBladeWidth       = meanBladeWidth       + semBladeWidth,
              meanPlusSEMBpRatio          = meanBpRatio          + semBpRatio,
              meanPlusSEMLwRatio          = meanLwRatio          + semLwRatio,
              meanPlusSEMBladeArea        = meanBladeArea        + semBladeArea,
              meanPlusSEMBladePerimeter   = meanBladePerimeter   + semBladePerimeter,
              meanPlusSEMBladeCircularity = meanBladeCircularity + semBladeCircularity,
              #
              meanMinusSEMPetioleLength    = meanPetioleLength    - semPetioleLength,
              meanMinusSEMBladeLength      = meanBladeLength      - semBladeLength,
              meanMinusSEMBladeWidth       = meanBladeWidth       - semBladeWidth,
              meanMinusSEMBpRatio          = meanBpRatio          - semBpRatio,
              meanMinusSEMLwRatio          = meanLwRatio          - semLwRatio,
              meanMinusSEMBladeArea        = meanBladeArea        - semBladeArea,
              meanMinusSEMBladePerimeter   = meanBladePerimeter   - semBladePerimeter,
              meanMinusSEMBladeCircularity = meanBladeCircularity - semBladeCircularity)
}

means   <- summarizeIt(grouped)
means2  <- summarizeIt(grouped2)
means33 <- summarizeIt(grouped33)
means35 <- summarizeIt(grouped35)
meansRef <- summarizeIt(groupedRef)
means1to8 <- summarizeIt(grouped1to8)

meansA   <- summarizeIt(groupedA)
meansB   <- summarizeIt(groupedB)
meansC   <- summarizeIt(groupedC)
meansD   <- summarizeIt(groupedD)
meansE   <- summarizeIt(groupedE)
meansF   <- summarizeIt(groupedF)
meansG   <- summarizeIt(groupedG)
meansH   <- summarizeIt(groupedH)

meansI   <- summarizeIt(groupedI)
meansJ   <- summarizeIt(groupedJ)
meansK   <- summarizeIt(groupedK)
meansL   <- summarizeIt(groupedL)
meansM   <- summarizeIt(groupedM)
meansN   <- summarizeIt(groupedN)
meansO   <- summarizeIt(groupedO)
meansP   <- summarizeIt(groupedP)

meansA2   <- summarizeIt(groupedA2)
meansB2   <- summarizeIt(groupedB2)
meansC2   <- summarizeIt(groupedC2)
meansD2   <- summarizeIt(groupedD2)
meansE2   <- summarizeIt(groupedE2)
meansF2   <- summarizeIt(groupedF2)
meansG2   <- summarizeIt(groupedG2)
meansH2   <- summarizeIt(groupedH2)

meansRef$set <- "a"
meansF$set <- "f"
meansG$set <- "g"
meansH$set <- "h"

meansI$set <- "i"
meansJ$set <- "j"
meansK$set <- "k"
meansL$set <- "l"
meansM$set <- "m"
meansN$set <- "n"
meansO$set <- "o"
meansP$set <- "p"

#' ** Leaf series plots: raw measurements
#' dev.new(height = 10, width = 10 * 16/9)
#'
#' *** Petiole
#' Later leaves get longer, until you get to the not-yet-fully-expanded ones.
aesPetL <- aes(leafjit,  petioleLength, fill = promoter2, col = promoter2)
ggplot(refLines, aesPetL) +
    xlab(leafStr) + ylab(petioleLengthStr) +
    geom_point(alpha = 0.5) +
    ## Used    alpha = 0.25 in the past.
    ## Also try alpha = 0.125, which is dimmer than I like.
    ## Difference between the more red and more yellow dots
    ## becomes hard to see when the points are that faint..
    geom_line(aes(group = transplantNum), alpha = 0.05) +
    geom_line(aes(leaf, meanPetioleLength), meansRef) + oneToTen + scm

ggplot(meansRef, aes(leaf, meanPetioleLength, col = promoter2)) +
    oneToTen + scm + geom_line() +
    geom_linerange(
        aes(ymin = meanMinusSEMPetioleLength,
            ymax = meanPlusSEMPetioleLength),
        meansRef, alpha = 0.1) +
    xlab(leafStr) + ylab(petioleLengthStr)


ggplot(leafj, aes(leaf, petioleLength, fill = promoter2, col = promoter2)) +
    xlab(leafStr) + ylab(petioleLengthStr) +
    geom_point(position = position_jitterdodge()) +
    geom_line(aes(leaf, meanPetioleLength), means) + oneToTen


#' *** Blade length
ggplot(refLines, aes(leafjit,  bladeLength, fill = promoter2, col = promoter2)) +
    xlab(leafStr) + ylab(bladeLengthStr) +
    geom_line(aes(group = transplantNum), alpha = 0.05) +
    geom_line(aes(leaf, meanBladeLength), meansRef) + scm +
    geom_point(alpha = 0.25) + oneToTen

ggplot(meansRef, aes(leaf, meanBladeLength, col = promoter2)) +
    oneToTen + scm + geom_line() +
    geom_errorbar(
        aes(ymin = meanMinusSEMBladeLength,
            ymax =  meanPlusSEMBladeLength),
        meansRef, #alpha = 0.1,
        width = 0.25) +
    xlab(leafStr) + ylab(bladeLengthStr)

ggplot(leafj, aes(leaf,  bladeLength, fill = promoter2, col = promoter2)) +
    xlab(leafStr) + ylab(bladeLengthStr) +
    geom_line(aes(leaf, meanBladeLength), means) +
    geom_point(position = position_jitterdodge()) + oneToTen


#' *** Blade width
ggplot(refLines, aes(leafjit,  bladeWidth, fill = promoter2, col = promoter2)) +
    xlab(leafStr) + ylab(bladeWidthStr) +
    geom_line(aes(leaf, meanBladeWidth), meansRef) +
    geom_line(aes(group = transplantNum), alpha = 0.05) + scm +
    geom_point(alpha = 0.25) + oneToTen

ggplot(meansRef, aes(leaf, meanBladeWidth, col = promoter2)) +
    oneToTen + scm + geom_line() +
    geom_errorbar(
        aes(ymin = meanMinusSEMBladeWidth,
            ymax =  meanPlusSEMBladeWidth),
        meansRef, width = 0.25) +
    xlab(leafStr) + ylab(bladeWidthStr)

ggplot(leafj, aes(leaf,  bladeWidth, fill = promoter2, col = promoter2)) +
    xlab(leafStr) + ylab(bladeWidthStr) +
    geom_line(aes(leaf, meanBladeWidth), means) +
    geom_point(position = position_jitterdodge()) + oneToTen


#' ** Leaf series plots: ratios
#' *** Blade length:Petiole length ratio
#'
#' Supplemental Figure S6:
#'
#' dev.new(height = 10, width = 10 * 16/9)

library(grid)

mnsblpl <-
    ggplot(meansI, aes(leaf, meanBpRatio, fill = promoter2, col = promoter2)) +
    bpRange +
    xlab(leafStr) + ylab(blpl) +
    geom_line() + oneToTen + scm  +
    geom_line(aes(leaf, meanBpRatio), meansJ) +
    geom_line(aes(leaf, meanBpRatio), meansK) +
    geom_line(aes(leaf, meanBpRatio), meansL) +
    geom_line(aes(leaf, meanBpRatio), meansM) +
    geom_line(aes(leaf, meanBpRatio), meansN) +
    geom_line(aes(leaf, meanBpRatio), meansO) +
    geom_line(aes(leaf, meanBpRatio), meansP) +
    facet_wrap( ~ set, ncol = 1)

allblpl <-
    ggplot(leafj, aes(leafjit,  bpRatio, fill = promoter2, col = promoter2)) +
    bpRange +
    xlab(leafStr) + ylab(blpl) +
    geom_point() + oneToTen + scm +
    geom_line(aes(group = transplantNum), alpha = 0.5) +
    facet_wrap( ~ promoter2, ncol = 1)

grid.newpage()
pushViewport(viewport(width = 1/2, x = 1/4))
print(mnsblpl, newpage = FALSE)
upViewport()
pushViewport(viewport(width = 1/2, x = 3/4))
print(allblpl, newpage = FALSE)

filter(leafj33b, bpRatio > 3)  # See notes in 33-dps_leaf-series-scans/readme.txt
# file:33-dps_leaf-series-scans/readme.txt::p-588


#' **** Still more: Simpler line graphs, with error bars
ggplot(meansRef, aes(leaf, meanBpRatio, col = promoter2)) +
    geom_line() +
    xlab(leafStr) + ylab(blpl) +
    geom_linerange(aes(leaf, ymin = meanMinusSEMBpRatio, ymax = meanPlusSEMBpRatio), meansRef) +
    oneToTen + scm


#' *** Blade length:Blade width ratio
#'
#' Figure 7:

mnslw <-
    ggplot(meansI, aes(leaf, meanLwRatio, fill = promoter2, col = promoter2)) +
    lwRange +
    xlab(leafStr) + ylab(blw) +
    geom_line() + oneToTen + scm  +
    geom_line(aes(leaf, meanLwRatio), meansJ) +
    geom_line(aes(leaf, meanLwRatio), meansK) +
    geom_line(aes(leaf, meanLwRatio), meansL) +
    geom_line(aes(leaf, meanLwRatio), meansM) +
    geom_line(aes(leaf, meanLwRatio), meansN) +
    geom_line(aes(leaf, meanLwRatio), meansO) +
    geom_line(aes(leaf, meanLwRatio), meansP) +
    facet_wrap( ~ set, ncol = 1)

alllw <-
    ggplot(leafj, aes(leafjit, lwRatio, fill = promoter2, col = promoter2)) +
    lwRange +
    xlab(leafStr) + ylab(blw) +
    geom_point() + oneToTen + scm +
    geom_line(aes(group = transplantNum), alpha = 0.5) +
    facet_wrap( ~ promoter2, ncol = 1)

grid.newpage()
pushViewport(viewport(width = 1/2, x = 1/4))
print(mnslw, newpage = FALSE)
upViewport()
pushViewport(viewport(width = 1/2, x = 3/4))
print(alllw, newpage = FALSE)


#' **** Similar: Simpler line graphs, with error bars
ggplot(meansRef, aes(leaf, meanLwRatio, col = promoter2)) +
    geom_line() +
    xlab(leafStr) + ylab(blw) +
    geom_linerange(aes(leaf, ymin = meanMinusSEMLwRatio, ymax = meanPlusSEMLwRatio), meansRef) +
    oneToTen + scm


#' ** Some plants of interest
filter(leafj, transplantNum == 528) %>% select(leaf, petioleLength, bladeLength,
                                               bladeWidth, bpRatio, lwRatio)
filter(leafj, transplantNum == 637) %>% select(leaf, petioleLength, bladeLength,
                                               bladeWidth, bpRatio, lwRatio)
#' filter(leafj, promoter == "422 bp", leaf == 10)

filter(leafj, transplantNum == 508) %>% select(leaf, petioleLength, bladeLength,
                                               bladeWidth, bpRatio, lwRatio)

#' ** How well do measurements from top-down RasPi photos correlate?
#' Started 2017-10-03.

#' library(dplyr)
library(tidyr)

coordTable <-
    read.delim(
        "../2016-08-26_zip/leaf-blade-base-coordinates.txt",
        na.strings = "")
## blank.lines.skip = TRUE) # by default

coordTable <- fill(coordTable, plantNum, type)

#' ** Smush things into a table by leaf.
centers <- filter(coordTable, type == "center")
centers <- select(centers, -type, -leafNum)
centers <- transmute(centers, plantNum = plantNum, x_center = x, y_center = y)

#' str(centers)
tips <- filter(coordTable, type == "tip")
tips <- select(tips, -type)
bases <- filter(coordTable, type == "base")
bases <- select(bases, -type)


tipbase <-
    left_join(tips, bases, by = c("plantNum", "leafNum", "image", "slice"),
              suffix = c("_tip", "_base"))

#' Duplication of center coordinates here.
tipbasecenter <-
    inner_join(tipbase, centers,
              by = c("plantNum")) # "image", "slice"))

tipbasecenter <- mutate(tipbasecenter,
                        bladeLpi =
                            sqrt((x_tip - x_base)^2 +
                                 (y_tip - y_base)^2),
                        petioleLpi =
                            sqrt((x_center - x_base)^2 +
                                 (y_center - y_base)^2))

relset <- filter(tipbasecenter, leafNum < 11)

bothtypes <-
    inner_join(leafj, relset, by = c("transplantNum" = "plantNum",
                                     "leaf" = "leafNum"))

########################################

#' Figure 3.9 in thesis:

petlab <- "Projected petiole length (pixels)"
blalab <- "Projected leaf blade length (pixels)"

#' dev.new(width = 7, height = 3)
grid.newpage()
pushViewport(viewport(width = 0.5, x = 0.3125))

print(
ggplot(bothtypes, aes(bladeLpi, bladeLength, col = leaf)) +
    geom_point(size = 1) + # alpha = 0.5 does not work great.
                                        # lightness matters.
    scale_color_gradient(guide = "legend") +
    ylab(paste0(bladeLengthStr, " from flatbed scans")) +
    xlab(blalab)
, newpage = FALSE)

upViewport()

pushViewport(viewport(width = 0.5, x = 0.75))

print(
ggplot(bothtypes, aes(petioleLpi, petioleLength, col = leaf)) +
    geom_point(size = 1) +
    scale_color_gradient(guide = "legend") +
    ylab(paste0(petioleLengthStr, " from flatbed scans")) +
    xlab(petlab)
, newpage = FALSE)
upViewport()

filter(relset, bladeLpi > 200)
filter(bothtypes, bladeLength < 0.51)

#' ** (end)
