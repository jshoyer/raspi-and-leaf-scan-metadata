
#' * Transplanting for another zip-1 complementation test
#'   with raspi imaging
#' mailto:j.s.hoyer@wustl.edu
#' Started 2016-09-04,
#' based on earlier plans in [[nb:zip-12-random-plan]].

#' Inputs: none.
#'   Relevant transplant/pot counts and identifiers
#'   are hardcoded in this script.
#'
#' Outputs:
#'   A table of the parent vial
#'   for each of the 180 transplants.
#'   This is displayed according to pot arrangement
#'   and also saved as a "tall" text file.
#'
#' The best way to view output, code, and comments may be to run
#' knitr::spin("10-dps_randomize-positions-for-transplanting.R")
#+ echo=FALSE
library(knitr)
opts_chunk$set(comment = '', echo = FALSE)
#library(dplyr)
library(grid)

#' Cf. other scripts,
#' listed in ../readme-23-raspi-imaging.txt.
#' In those scripts we imagined pot positions
#' from the perspective of the camera looking down.
#' This was error-prone,
#' so let us instead consider positions
#' from the perspective of a person looking into the chamber.
#' This is just like rotating our axes by 90°,
#' so that they fit R's matrix indexing order more naturally.

#' * Define groupings
# Picked a seed based on
# sample(1:1000, 1)
set.seed(447)  # for randomization, below.

#' Half of the genotypes might be considered
#' 'reference' genotyes,
#' and the other half are not.
#' See [[sop:3xHA-AGO7-strains]].
ref <-  c(204,  768, 2833, 2834)
nonref <- c(2204, 2731, 2743, 2758)

#' There are eight plant groupings of interest:
#' each has one of the eight genotypes represented only once
#' and the other genotypes represented by two plants,
#' for a total of fifteen plants in the flat.
#' Ignore the arrangement of the genotypes/pots within the flat --
#' we deal with that below.
#'
#' The first four plant groupings (#1 to 4)
#' have a different reference genotype present only once.
grouping1 <- c(ref, ref[-1], nonref, nonref)
grouping2 <- c(ref, ref[-2], nonref, nonref)
grouping3 <- c(ref, ref[-3], nonref, nonref)
grouping4 <- c(ref, ref[-4], nonref, nonref)
grouping5 <- c(ref, ref, nonref, nonref[-1])
grouping6 <- c(ref, ref, nonref, nonref[-2])
grouping7 <- c(ref, ref, nonref, nonref[-3])
grouping8 <- c(ref, ref, nonref, nonref[-4])

groupingList <-
    list(grouping1, grouping2, grouping3, grouping4,
         grouping5, grouping6, grouping7, grouping8)


#' * Start randomization
#' There is only one flat for each of these groupings,
#' because we want a slightly higher sample size
#' for the reference genotypes:
#' 22 plants, versus 21 plants for each of the nonreference genotypes.
#' We want two of these flats on the top half of the growth chamber
#' and two on the lower half.
missing1ref <- sample(1:4)

#' We also want to balance the other four plant groupings
#' (call them groupings #5 to 8)
#' so that for a given genotype
#' there is one flat with that genotype respresented once
#' on the top half of the chamber
#' and one flat on the bottom half.
#+ results="asis"
orderOfPlantGroupings <- c(
    sample(c(5:8, missing1ref[1:2])),
    sample(c(5:8, missing1ref[3:4])))


parentVial <- numeric()
flatCounter <- 0

for (i in orderOfPlantGroupings) {
    flatCounter <- flatCounter + 1
    cat(paste0("Flat ", flatCounter,
                ", plant grouping ", i))
    pots <- sample(groupingList[[i]])
    print(kable(matrix(pots, nrow = 5), col.names = 1:3),
          format = "html")
    parentVial <- c(parentVial, pots)
}

# print(parentVial)

plantNumber <- c(
    501:590,
    601:690)


potPos  <- c(rep(1:15, 12))
flatPos <- rep(1:12, each = 15)

plants <- data.frame(plantNumber, parentVial, potPos, flatPos)

#' ** Save two copies of table, for convenience:
dir.create("working-data", showWarnings = FALSE)
write.table(plants,
            sep = "\t", quote = FALSE, row.names = FALSE,
    "working-data/transplants-in-order-by-pot-position_2016-09-05.txt")

plants <- dplyr::arrange(plants, parentVial)
write.table(plants,
            sep = "\t", quote = FALSE, row.names = FALSE,
    "working-data/transplants-in-order-by-parent-vial_2016-09-05.txt")
# Put these files into 9/5 excel file as sheets.

#' Printed out tables above,
#' wrote down source plate numbers as I transplanted.
#' (Just printed the whole thing --
#'  did not mess with getting the breaks between tables right,
#'  in contrast to previous efforts.)
#' Stapled in to notebook 10 p. 60.
#' Cf. p. 47, p. 23, and p. 16.

#' ** (I later made many nicer diagrams
#' with library(grid)
#' and used them to guide selection of plants
#' for dissection at 33-dps
#' and sampling   at 38-dps
#' See pot-indexing-drawing-helper-functions.R
#' and pot-indices-and-diagrams.R

#' * (end)
