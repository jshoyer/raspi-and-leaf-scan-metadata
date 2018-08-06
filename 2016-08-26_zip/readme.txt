
** Overview
As described in our manuscript,
the main question of interest for this experiment
was how sequence upstream of the AGO7 transcription start site
might be sufficient (when driving the AGO7 coding sequence)
for complementation of ago7 (zip-1) mutants.

Transplants measured in this experiment
were designated #401 to 490 (upper level of chamber)
and #501 to 590 (lower level of chamber).

- Top-down photos: https://zenodo.org/record/1322154
  - Approximate pi photography configuration for this experiment:
    https://zenodo.org/record/439652
- Leaf scans and Canon camera photos:
  https://zenodo.org/record/1322799


** Figures 2 and 3 in the PlantCV 2 paper used a photo from this experiment
https://doi.org/10.7717/peerj.4088

fov-01/2016-09-27_1025_ch129-pos01.jpg
(fov-01/2016-09-27_1030_ch129-pos01.jpg
 in the Zenodo upload
 is essentially the same.)

See under image-analysis-scripts/figure.2and3.multiplant.watershed/
at https://github.com/danforthcenter/plantcv-v2-paper/

** file:10-dps_randomize-positions-for-transplanting.R
Transplant [2016-09-05 Mon]
(see timetable below)

'dps' = days post-stratification


** file:28-to-30-dps_abaxial-trichome-positions-and-true-leaf-6-caliper-measurements.txt
Started 2016-09-23.
That sheet was adapted from a table
output by the position randomization script linked to above.
I grayed out the 'source plate' and 'parent vial' columns
in the excel sheet to try to avoid subconscious bias during data collection,
and then pasted in letter codes instead, randomly assigned with
file:30-dps_trichome-and-leaf-6-analysis.R

** Sets of photos of individual plants taken with Canon EOS Digital Rebel XT 350D camera

I used fully automatic shooting mode for all of these sets of photos,
so unfortunately the exposure times (and possibly other settings) were quite variable.

Photos enable checking that leaves were removed
for flat-bed scanning
in the correct phyllotactic order.


* Timeline
** Phases

The experiment can thought of as having six phases:
1. Initial imaging of all 180 plants/pots
2. Imaging continues after slight position disruption for each flat
   in the range of 28 to 30 days post-stratification.
3. Imaging of 156 pots,
   after 24 were dissected
4. Imaging of 132 pots,
   after 24 more were dissected
5. Imaging of 68 pots,
   after 64 plants were removed for tissue sampling.
6. Imaging of the final 29 plants,
   after 39 empty pots/multi-rosette plants were pulled
The phases started at slightly different times for different plants.

Two files are kept around for historical reasons
(until I find time to potentially streamline
 30-dps_trichome-and-leaf-6-analysis.R)
but are not important to the main analysis:
33-dps_second-set-of-abaxial-trichome-measurements.txt and
38-dps_plants-for-leaf-sampling.txt


** Table by day

| [2016-08-24 Wed] |    | Sterilize seed                             |    |
| [2016-08-25 Thu] |    |                                            |    |
| [2016-08-26 Fri] |    | Plate seed                                 |    |
| [2016-08-27 Sat] |  1 |                                            |    |
| [2016-08-28 Sun] |  2 |                                            |    |
| [2016-08-29 Mon] |  3 |                                            |    |
| [2016-08-30 Tue] |  4 |                                            |    |
| [2016-08-31 Wed] |  5 |                                            |    |
| [2016-09-01 Thu] |  6 |                                            |    |
| [2016-09-02 Fri] |  7 |                                            |    |
| [2016-09-03 Sat] |  8 |                                            |    |
| [2016-09-04 Sun] |  9 |                                            |    |
| [2016-09-05 Mon] | 10 | Transplant                                 |    |
| [2016-09-06 Tue] | 11 |                                            |    |
| [2016-09-07 Wed] | 12 | Start of time-lapse photography            |  1 |
| [2016-09-08 Thu] | 13 |                                            |  2 |
| [2016-09-09 Fri] | 14 |                                            |  3 |
| [2016-09-10 Sat] | 15 |                                            |  4 |
| [2016-09-11 Sun] | 16 |                                            |  5 |
| [2016-09-12 Mon] | 17 |                                            |  6 |
| [2016-09-13 Tue] | 18 |                                            |  7 |
| [2016-09-14 Wed] | 19 |                                            |  8 |
| [2016-09-15 Thu] | 20 |                                            |  9 |
| [2016-09-16 Fri] | 21 |                                            | 10 |
| [2016-09-17 Sat] | 22 |                                            | 11 |
| [2016-09-18 Sun] | 23 |                                            | 12 |
| [2016-09-19 Mon] | 24 |                                            | 13 |
| [2016-09-20 Tue] | 25 |                                            | 14 |
| [2016-09-21 Wed] | 26 |                                            | 15 |
| [2016-09-22 Thu] | 27 |                                            | 16 |
| [2016-09-23 Fri] | 28 | Measure                                    | 17 |
| [2016-09-24 Sat] | 29 | Measure                                    | 18 |
| [2016-09-25 Sun] | 30 | Measure                                    | 19 |
| [2016-09-26 Mon] | 31 |                                            | 20 |
| [2016-09-27 Tue] | 32 |                                            | 21 |
| [2016-09-28 Wed] | 33 | Dissect and scan                           | 22 |
| [2016-09-29 Thu] | 34 |                                            | 23 |
| [2016-09-30 Fri] | 35 | Dissect and scan                           | 24 |
| [2016-10-01 Sat] | 36 |                                            | 25 |
| [2016-10-02 Sun] | 37 |                                            | 26 |
| [2016-10-03 Mon] | 38 | Collect tissue                             | 27 |
| [2016-10-04 Tue] | 39 |                                            | 28 |
| [2016-10-05 Wed] | 40 |                                            | 29 |
| [2016-10-06 Thu] | 41 |                                            | 30 |
| [2016-10-07 Fri] | 42 |                                            | 31 |
| [2016-10-08 Sat] | 43 |                                            | 32 |
| [2016-10-09 Sun] | 44 | Pull empty pots                            | 33 |
| [2016-10-10 Mon] | 45 | Stop photo capture, photograph color cards | 34 |


* Initial (human-readable) start on MIAPPE 1.0 information
Per Ćwiek-Kupczyńska et al. 2016
(http://doi.org/10.1186/s13007-016-0144-4).
See http://cropnet.pl/phenotypes/wp-content/uploads/2016/04/MIAPPE.pdf

1. General
   a. Unique identifier: 2016-08-26_zip
   b. Title: Transgenic complementation experiment 2016-08-26_zip
   c. Description: Tested for leaf shape defects in 8 genotypes.
2. Timing and location
   a. Start of experiment: 2016-08-26 (removal of seed from stratification)
   b. Duration: 45 days, until 2016-10-10
   c. Geographic location: Danforth Center Plant Growth Facility
3. Biosource
   a. Organism: Arabidopsis thaliana
   c. Seed source: Carrington lab seed stocks,
      derived from lines from the ABRC (Col-0 CS60000)
      and the R.S. Poethig lab (zip-1 mutant; Hunter et al. 2003)

4. Environment
   a. Growth facility: DDPSC PGF Chamber 129
      (Conviron MTR25 reach-in chamber
      with PolyLux fluorescent bulbs.)
   c. Air humidity: 50%
   d. Daily photon flux: 200 µmol/(s * m^2) * (60 * 60 * 8) = 5.76 mol/(s * m^2)
   e. Temperature: 22 °C
   f. Rooting medium: ProMix FPX growth mix
   g. Container type: "1801 pot", fifteen pots per flat.
   j. Number of plants per container: 1 plant per pot
   k. pH: approximately neutral---not measured directly.

5. Treatments: none.
6. Design: see 10-dps_randomize-positions-for-transplanting.R
