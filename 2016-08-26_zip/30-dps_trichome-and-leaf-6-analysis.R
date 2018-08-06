
#' * zip experiment 12
#' mailto:j.s.hoyer@wustl.edu
#' Started 2016-09-25.
#'
#' Cf. ../2016-05-13_zip/abaxial-trichome-and-leaf-6-analysis.R
#' That directory has the table defining what the parent genotypes are.

#' ** INPUTS and dependencies:
library(digest)
library(ggplot2); theme_set(theme_classic())
library(grid)
library(lattice)
source("../R/axis-label-strings-and-helper-function.R")
                                        # Includes workupTrichCaliper()
                                        # Loads dplyr

gen <- read.delim("../2016-08-26_zip/00_parent-vial-genotypes.txt")

#' Table copied from 2016-09-05 excel file:
data923 <- read.delim(
    "28-to-30-dps_abaxial-trichome-positions-and-true-leaf-6-caliper-measurements.txt",
    na.strings = "")
codesInOrder <- readLines("00_tested-genotype-codes.txt")
#codesAlt <- c("WT", "zip-1", codesInOrder[3:8])

#' ** Workup etc.
data923 <- mutate(data923, parentVial = as.factor(parentVial))
data923 <- inner_join(data923, gen, by = c("parentVial" = "vial"))
data923$promoter <- factor(data923$promoter, levels = codesInOrder)
                                       # get rid of unused levels
data923$promoter <- factor(data923$promoter, levels = codesInOrder)
data923$promoter2 <- data923$promoter
levels(data923$promoter2)[1] <- "WT"
levels(data923$promoter2)[2] <- "zip-1"
data923$promoter2 <- factor(data923$promoter2,  # Alternate order/notation
                            levels = c("WT", codesInOrder[3:8], "zip-1"))

#' For convenience, when we do not want to use vial numbers:
dir.create("working-data")
writeLines(as.character(data923$promoter),
           "working-data/genotype-coded-plants-in-order.txt")
digest("working-data/genotype-coded-plants-in-order.txt",
       algo = "sha1")
# 5ada5dc3a4cc6d451329a8c9a54e57f70fe29555
writeLines(as.character(data923$promoter2),
           "working-data/genotype-coded-plants-in-order-2.txt")
digest("working-data/genotype-coded-plants-in-order-2.txt",
       algo = "sha1")
# c0c989a0e15cf49fd17cda05afeeb7d0e9112df0

#' I want to blind myself better during analysis,
#' because I remember what several of the vial numbers are now.
set.seed(117)
## I picked the number 117 by doing sample(1:1000, 1)
AtoH <- LETTERS[1:8]

data923blind <- select(data923,
                       -transgene, -background,
                       -parentVial, -sourcePlate
                       -notes)
data923blind$promoter <- factor(data923blind$promoter, labels = sample(AtoH))
data923blind$promoter <- factor(data923blind$promoter,
                                levels = AtoH)         # Scramble order in graphs
data923blind$flatPos <- rep(1:12, each = 15) # from below/old script

# To break the blind, I think one could do this:
#' set.seed(117)
#' sam <- sample(1:8)
#' cbind(LETTERS[sam], codesInOrder)
#'
#' Could also just look at (e.g.)
#' working-data/worked-up-data-simililar-to-2016-09-05-spreadsheet.csv

#' dir.create("working-data")
#' writeLines(as.character(data923blind$promoter),
#'            "working-data/letter-codes-in-order.txt")

notesOnly923 <- select(data923, plantNum, parentVial, promoter, notes)
data923 <- select(data923, -notes)  # Too wide for R session....

# Copied from 10-dps_randomize-positions-for-transplanting.R:
data923$potPos  <- c(rep(1:15, 12))
data923$flatPos <- rep(1:12, each = 15)

data923$dayScored <- rep(c("28 dps", "29 dps", "30 dps"),
                 each = 15 * 4)
                 ## I measured four flats on each day.

#' ** How were the non-measured plants distributed across the flats?
codedOfInterest <- factor(rep(NA, 180), levels = AtoH)
rowsToInclude <- is.na(data923$exclude)
codedOfInterest[rowsToInclude] <- data923blind$promoter[rowsToInclude]
excluded <- filter(data923, !rowsToInclude)
data923main <- filter(data923, rowsToInclude)
with(data923main, table(parentVial))
with(excluded,    table(parentVial))
#' with(data923main, table(promoter))
#' with(excluded,    table(promoter))

# What I wanted:
with(data923, table(promoter, flatPos))
# What I got:
with(data923main, table(promoter, flatPos))
with(data923main, table(flatPos))
# Pattern of missingness/failure to thrive:
# reasonably even, except for Flat 8,
# in which many plants failed to thrive.
with(excluded,    table(promoter, flatPos))
with(excluded,    table(flatPos))

data923$codedOfInterest <- codedOfInterest
outcsvfile1 <- "../2016-08-26_zip/working-data/better-labeled-923.csv"
write.csv(data923, outcsvfile1, row.names = FALSE)
digest(algo = "sha1", file = outcsvfile1)
# e69166dfff8d29e5ae7dad56d7c20269f9cdcaaf
#' Read in by analysis-7.R, but not really used.


#' ** How well did my quick binary classification as zippy/nonzippy work?
#' Not perfect, but pretty good:
with(data923, table(promoter, zippy))
with(data923, table(promoter, zippy,     useNA = "ifany"))
with(data923main, table(promoter, zippy)) # NAs already removed
with(excluded, table(promoter, zippy))

#' What about the plants
#' that did not match the general pattern for their genotype?
#' Were they borderline and/or difficult-to-classify cases?
#' options(width = 300)
filter(data923, promoter ==     "GUS", zippy == "no")
filter(data923, promoter == "1934 bp", zippy == "yes")
#' Hmm. Not sure how I could have marked 517 as nonzippy.
#' It had a very strong phenotype,
#' which became increasingly evident as it aged.
#' Likewise, it is very clear
#' that plant 545 displayed a mutant phenotype.

filter(data923, promoter ==  "495 bp", zippy == "yes")
filter(data923, promoter ==  "453 bp", zippy == "yes")
#' Plant 528 had a complicated phyllotaxy.
filter(data923, promoter ==  "422 bp", zippy == "yes")
filter(data923, promoter ==  "298 bp", zippy == "no")
#' Plant 507 had a very broad leaf 3, but very narrow leaves later on.
filter(data923, promoter ==    "0 bp", zippy == "no")

#' * Abaxial trichomes analysis
#' *** Graphs
workedUp923 <- workupTrichCaliper(data923main)
                                        # Script sourced above

outcsvfile <- "../2016-08-26_zip/working-data/workedUp923.csv"
write.csv(workedUp923, outcsvfile, row.names = FALSE)
digest(algo = "sha1", file = outcsvfile)
# 35fb87512213a58a59686055369a817bddb5af7e
#' Read in by analysis-7.R

with(workedUp923, table(promoter, firstLeaf))

grouped923 <- group_by(workedUp923, promoter)
trichmeansd <-
    summarize(grouped923,
              meanFirst = mean(firstLeaf),
              sdFirst   = sd(firstLeaf),
              meanBp    = mean(bpratio),
              sdBp      = sd(bpratio),
              meanLw    = mean(lwratio),
              sdLw      = sd(lwratio))

## Redundant: kept to avoid having to edit code for scatterplots below,
## which apparently broke anyway?
bpmeansd <-
    summarize(grouped923, mean(bpratio), sd(bpratio))
lwmeansd <-
    summarize(grouped923, mean(lwratio), sd(lwratio))

#' Cf. ../2015-05_zip/analysis.R
#' (on branch).
#' Plant 626 first trichome position was initially not written down:
#' filter(workedUp923, is.na(firstLeaf))
#' No longer needed after rescoring 9/26: workedUp923 <- filter(workedUp923, !is.na(firstLeaf))

#' Method 1: clumsy
ab <- ggplot(workedUp923, aes(as.factor(firstLeaf)))
#' Method 2: extends x-axis unnecessarily. (See below)
ab2 <- ggplot(workedUp923, aes(firstLeaf))

zft <- c(5, 10)

#' dev.new(width = 10 * 16/9, height = 10)
grid.newpage()
pushViewport(viewport(width = 0.95, height = 0.95))
pushViewport(viewport(x = 1, just = "right", width = 2/3))

print(newpage = FALSE,
ab + geom_histogram(fill = "gray", stat = "count") +
    facet_grid(promoter ~ .) + xlab(indx) +
    scale_y_discrete(limits = zft) +
    geom_vline(data = trichmeansd, aes(xintercept = meanFirst - 3),
                                        # Manually subtract    /\
                linetype="dashed", size = 0.5)
)

#' I cannot figure out how to set the x-limits
#' for the other version of the plot, below.
#' I currently have to add the promoter diagrams manually.

upViewport()
pushViewport(viewport(x = 0, just = "left", width = 1/3))
pushViewport(viewport(y = 0.525, height = 0.95, width = 0.9))
## A truncated representation of the 1934 bp promoter
## would be a more efficient use of space.
lengths <- c(1934, 495, 453, 422, 298, 0)
for (i in 1:6) {
    grid.lines(x = c(1 - lengths[i]/lengths[1], 1),
               y = (6.5 - i)/8)
}

# todo:
#       properly round mean/SD and draw it as a table directly next to graph

ab2 + geom_bar(fill = "gray") +
    facet_grid(promoter ~ .) + xlab(indx) +
    scale_y_discrete(limits = zft) +
    scale_x_discrete() +
    geom_vline(data = trichmeansd, aes(xintercept = meanFirst),
                linetype="dashed", size = 0.5)


#' Adapted from ../2016-05-13_zip/abaxial-trichome-and-leaf-6-analysis.R
#' ** Formal statistical analysis
summary(fm0 <- aov(firstLeaf ~ promoter, data = workedUp923))
(hsd0 <- TukeyHSD(fm0))
head(hsd0$promoter)
summaryStatFindings <- "
Abaxial trichomes appeared 1.7 leaf positions earlier on average
for mutant vs. wild-type
empty-vector-transformed plants
(95% confidence interval 0.5 to 2.9,
 p = 4 * 10^-4,
 Tukey's honest significant difference method),
consistent with previous reports.
No 3xHA-AGO7 transgenic line
showed a detectable increase
in earliest abaxial trichome position
(relative to empty-vector-transformed mutant plants;
 p > 0.3),
indicating that none of the promoter lengths tested
were able to drive full complementation of this defect.
"

#' plot(hsd0)

#' *** Try Benjamini-Hochberg adjustment
#' Revisit after discussing with Jeff on 2017-05-05.
#' He pointed out that the HSD method is conservative.
lvl <- levels(workedUp923$promoter2)

pvals <- numeric()
pvalsNegC <- numeric()
for (i in 1:7) {
    for (j in (i+1):8) {
        tresult <-
            with(workedUp923,
                 t.test(firstLeaf[promoter2 == lvl[i]],
                        firstLeaf[promoter2 == lvl[j]],
                        var.equal = TRUE))
        p <- tresult$p.value
        print(paste0(
            lvl[i], " vs. ", lvl[j], ": unadjusted p = ", p))
        pvals <- c(pvals, p)
        if (j == 8) {
            pvalsNegC <- c(pvalsNegC, p)
        }
    }
}

#' check:
#' grp1 <- filter(workedUp923, promoter2 == "WT")
#' grp2 <- filter(workedUp923, promoter2 == "1934 bp")
#' t.test(grp1$firstLeaf, grp2$firstLeaf, var.equal = TRUE)
#' pvals[1]

adjusted <- p.adjust(pvals, "BH")
# See grid: file:readme.txt::#trichome-BH-adjusted-p-values

#' I was curious about the numbers
#' for two slightly more conservative adjustments:
#' p.adjust(c(pvals[1:6], pvalsNegC), "holm")
#' p.adjust(c(pvals[1:6], pvalsNegC), "hommel")

#' * Leaf caliper measurements
#' ** General checks
#' These do not work well at all because there are *eight* genotypes.
#' Lattice starts reusing colors for multiple groups.
#' Therefore for now I just drop one of the groups:
noZero923 <- filter(workedUp923, promoter != "0 bp")
noZero923$promoter <- factor(noZero923$promoter, levels = codesInOrder[-8])
#' head(workedUp923)
#' head(noZero923)

#' How well do shape measurements correlate with trichome positions?
xyplot(bladeL   ~ jitter(firstLeaf),
       data = noZero923, groups = promoter, auto.key = TRUE)
xyplot(bladeW   ~ jitter(firstLeaf),
       data = noZero923, groups = promoter, auto.key = TRUE)
xyplot(petioleL ~ jitter(firstLeaf),
       data = noZero923, groups = promoter, auto.key = TRUE)

## See below
xyplot(bpratio   ~ jitter(firstLeaf),
       data = noZero923, groups = promoter, auto.key = TRUE)
xyplot(lwratio   ~ jitter(firstLeaf),
       data = noZero923, groups = promoter, auto.key = TRUE)

xyplot(bladeW   ~ petioleL,
       data = noZero923, groups = promoter, auto.key = TRUE)
xyplot(bladeW   ~ bladeL,
       data = workedUp923, groups = promoter)
## Tightest correlation:
xyplot(petioleL ~ bladeL,
       data = noZero923, groups = promoter, auto.key = TRUE)


#' ** Correlation between phenotype means
plot(trichmeansd$meanLw, trichmeansd$meanBp, type = "n",
     xlab = blw,    ylab = blpl)
text(trichmeansd$meanLw, trichmeansd$meanBp, codesInOrder)

#' In some sense the 422 bp promoter transgene
#' can complement the leaf shape phenotype
#' *better* than it can complement the trichome position phenotype.
#' (Discussed with Jim 2016-09-26)
plot(trichmeansd$meanFirst, lwmeansd$mean, type = "n",
     xlab = indx,      ylab = blpl)
text(trichmeansd$meanFirst, lwmeansd$mean, codesInOrder)

plot(trichmeansd$meanFirst, bpmeansd$mean, type = "n",
     xlab = indx,      ylab = blw)
text(trichmeansd$meanFirst, bpmeansd$mean, codesInOrder)


#' ** Univariate scatterplot checks
workedUp923rev <- workedUp923
workedUp923rev$promoter <- factor(workedUp923$promoter, levels = rev(codesInOrder))

dotplot(promoter ~ bladeW, data = workedUp923rev,
       xlab = bw)
dotplot(promoter ~ bladeL, data = workedUp923rev,
       xlab = bl)
filter(workedUp923, bladeL > 25)
dotplot(promoter ~ bladeL, group = zippy, data = workedUp923rev,
       xlab = pl, pch = 19, alpha = 0.5)
dotplot(promoter ~ petioleL, data = workedUp923rev,
       xlab = pl)
dotplot(promoter ~ petioleL, group = zippy, data = workedUp923rev,
       xlab = pl, pch = 19, alpha = 0.5)
filter(workedUp923, petioleL < 5.75)

#' Actual graphs of interest:
dotplot(promoter ~ bpratio, data = workedUp923rev,
       xlab = blpl)
#' Positions for severally of the plants scored as 'zippy'
#' is somewhat puzzling:
dotplot(promoter ~ bpratio, group = zippy, data = workedUp923rev,
       xlab = blpl, pch = 19, alpha = 0.5)
filter(workedUp923, bpratio > 2)
workedUp923minusOutliers <- filter(workedUp923rev, bpratio < 2)

dotplot(promoter ~ lwratio, data = workedUp923rev,
       xlab = blw)
dotplot(promoter ~ lwratio, group = zippy, data = workedUp923rev,
       xlab = blw, pch = 19, alpha = 0.5)

filter(workedUp923, lwratio > 3)
filter(workedUp923, lwratio < 1.25)

######################################################################
# The same thing, but nicer looking,
# and with grid graphics stuff
# for ~/r/thesis-co-1/reports/captions/fig-leaf-6-caliper-lw.tex

#' dev.new(width = 6.5, height = 6.5 * 2/3)

grid.newpage()
pushViewport(viewport(width = 0.45, x = 0.75))

print(
ggplot(workedUp923rev, aes(bpratio, promoter)) +
    ylab("") + xlab(blpl) +
    geom_point(alpha = 0.25) +
    geom_point(aes(meanBp, promoter), trichmeansd,
               shape = "|", col = "red", size = 4)
, newpage = FALSE)

upViewport()

# (skip this one)
#' ggplot(workedUp923minusOutliers, aes(bpratio, promoter)) +
#'     ylab("") + xlab(blpl) +
#'     geom_point(alpha = 0.25) +
#'     geom_point(aes(meanBp, promoter), trichmeansd,
#'                shape = "|", col = "red", size = 4)

pushViewport(viewport(width = 0.45, x = 0.35))

print(
ggplot(workedUp923rev, aes(lwratio, promoter)) +
    ylab("") + xlab(blw) +
    geom_point(alpha = 0.25) +
    geom_point(aes(meanLw, promoter), trichmeansd,
               shape = "|", col = "red", size = 4)
, newpage = FALSE)

upViewport()
######################################################################

dotplot(promoter ~ bpratio, group = dayScored, data = workedUp923rev,
       xlab = blpl, auto.key = TRUE)
#' filter(workedUp923, bpratio > 2)
dotplot(promoter ~ lwratio,  group = dayScored, data = workedUp923rev,
       xlab = blw)

#' Most of the plants
#' with an apparently very high blade length:petiole length ratio
#' were scored on the first of the three days.
#' I will examine those plants in detail soon.
dotplot(promoter ~ bpratio | dayScored, data = workedUp923rev,
       xlab = blpl, auto.key = TRUE)
dotplot(promoter ~ lwratio | dayScored, data = workedUp923rev,
       xlab = blw)
#' (No other patterns jumped out at me.)


#' ** Formal statistical analysis
summary(fm1 <- aov(bpratio ~ promoter, data = workedUp923))
(hsd1 <- TukeyHSD(fm1))
plot(hsd1)

summary(fm2 <- aov(lwratio ~ promoter, data = workedUp923))
(hsd2 <- TukeyHSD(fm2))
plot(hsd2)


#' ** How well do trichome positions scored on 33 and 35 dps line up with the 28 to 30 dps ones?
day33 <- read.delim(
    "../2016-08-26_zip/33-dps_leaf-series-scans/abaxial-trichome-positions.txt")

day33 <- rename(day33, firstLeaf33 = firstLeaf)
firstSetSampled <- inner_join(workedUp923, day33)
xyplot(firstLeaf ~ jitter(firstLeaf33), data = firstSetSampled)
firstSetSampled <- mutate(firstSetSampled,
                          firstLeafDiff = firstLeaf33 - firstLeaf)
select(firstSetSampled, starts_with("firstLeaf"))
length(which(firstSetSampled$firstLeaf == firstSetSampled$firstLeaf33))

#' ** More subsets
sampledAlready <- data923$plantNum %in% firstSetSampled$plantNum
codedOfInterest2 <- codedOfInterest
codedOfInterest2[sampledAlready] <- NA

data923blind <- mutate(data923blind,
                       sampledAlready = plantNum %in% firstSetSampled$plantNum)

secondRosetteLate <- scan(text = "
524
525
552
603
641
658")

data923blind <- mutate(data923blind,
                       exclude2 = plantNum %in% secondRosetteLate)

codedOfInterest3 <- codedOfInterest2
codedOfInterest3[data923blind$exclude2] <- NA


data923blindMain <- filter(data923blind, is.na(exclude))
with(data923blindMain, table(promoter))
with(data923blindMain, table(promoter, flatPos))


todo <- filter(data923blindMain, !sampledAlready)
todoStringent <- filter(todo, !exclude2)

with(todo, table(promoter))
with(todoStringent, table(promoter))

with(todo, table(flatPos))
with(todoStringent, table(flatPos))

#' Flat 3 is perfectly balanced: one plant of each genotype.
#' Flat 2 could be perfectly balanced
#' if we do not exclude multirosette plants 524 and 525.
#' (I opted to sample those two plants,
#'  as noted on p. 68.)
with(todo, table(promoter, flatPos))
with(todoStringent, table(promoter, flatPos))

write.table(todo, "working-data/remaining-plants_2016-10-03.tsv")
write.table(todoStringent, "working-data/remaining-good-looking-plants_2016-10-03.tsv")

secondSetSampled <- read.delim(
    "../2016-08-26_zip/38-dps_sample-leaves/list-of-plants-sampled.txt",
    header = FALSE)
#' Who is left after sampling those plants?
data923blind <- mutate(data923blind,
                       sampledAlready2 = plantNum %in%
                           secondSetSampled$V1)
codedOfInterest4 <- codedOfInterest3
codedOfInterest4[data923blind$sampledAlready2] <- NA
table(codedOfInterest4)

#' The following includes 23 plants that I photographed on 2016-10-10,
#' but misses six more which were too small to score at 28 dps.
todo2 <- filter(data923blind,
                is.na(exclude),
                !(exclude2),
                !sampledAlready,
                !sampledAlready2)
#' I did not score trichomes for these six plants,
#' but I decided they were worth photographing anyway:
smallRecovered <- c(505, 517, 531,
                    605, # sort of
                    677,
                    688) # sort of

day45photoFiles <-
    dir("../2016-08-26_zip/45-dps_photos/", pattern = "plant-")

thirdSet <- substring(day45photoFiles, 7, last = 9)

misfits <- filter(data923, !(plantNum %in%
                             c(firstSetSampled$plantNum,
                               secondSetSampled$V1,
                               thirdSet)))
writeLines(as.character(misfits$plantNum),
    "../2016-08-26_zip/working-data/misfits.txt")

#' 48 + 64 + 29 + 39
#' 48 + 29      # < these two sets fit together in a single freezer box,
                                        # so continue numbering: 49 to 77

#' Plant 504 and 507 looked fine up to the point of trichome scoring
#' (though 507 had a rather large leaf 3),
#' but developed a second rosette
#' (or perhaps several more,
#'  and/or many axillary leaves),
#' so I removed them in the morning, 44 dps.
#' These two plants were marked outliers for blade:petiole ratio.
#' Excluding those plants seems justifiable,
#' but I am not sure what is best.
#' See also the list above of six plants (two sampled)
#' that developed a second rosette late.
filter(data923, plantNum %in% misfits$plantNum) %>%
    select(plantNum, exclude)

#' ** (end)
