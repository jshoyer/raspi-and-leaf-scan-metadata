
#' mailto:j.s.hoyer@wustl.edu
#' Started 2016-06-19.
#' Adapted 2016-08-29.

#' ** Load packages; Define axis label strings
library(dplyr)
#library(lattice)

indx <- "Index of first true leaf with at least one abaxial trichome"

bw <- "Blade width (mm)"
pl <- "Petiole length (mm)"
bl <- "Blade length (mm)"
blw  <- "Blade length:Blade width ratio"
blpl <- "Blade length:Petiole length ratio"

#' ** Example INPUTS:
#' Dataset 1:
#' includes phyllotaxy directions (hence 'pd').
## pd <- read.delim("../2016-05-13_zip/33-dps-abaxial-trichome-positions-and-subsequent-leaf-6-caliper-measurements.txt")

gen <- read.delim("../2016-05-13_zip/parent-vial-genotypes.txt")


#' ** Workup
workupTrichCaliper <- function(
                               pd,
                               ...)
{
#pd <- inner_join(pd, gen, by = c("parentVial" = "vial"))

#notesOnly <- select(pd, plantNum, parentVial, promoter, notes)
#pd <- select(pd, -notes)  # Too wide for R session....
pd <- mutate(pd, lwratio = bladeL / bladeW)
pd <- mutate(pd, bpratio = bladeL / petioleL)

    return(pd)
}
