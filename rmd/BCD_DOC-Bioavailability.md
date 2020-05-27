BCD and DOC Bioavailability
================
Nicholas Baetge
5/26/2020

# Intro

This document shows plots and tables of the merged field-experiment
NAAMES DOC data.

``` r
library(tidyverse) 
library(rmarkdown)
library(knitr)
library(readxl)
library(data.table) 
library(scales)
library(zoo)
library(oce)
library(patchwork)
#rmarkdown tables
library(stargazer)
library(pander)
library(growthcurver)
#stat tests
library(lmtest)
library(lmodel2)
library(rstatix)
library(ggpubr)
```

# Import Data

``` r
doc <- read_rds("~/naames_bioav_ms/Output/processed_bioavailability.rds") %>% 
   select(-c(Bottle, doc, interp_doc)) %>% 
  arrange(Cruise, Station, Hours) %>% 
  distinct() %>% 
  mutate(bge = ifelse(Cruise == "AT39" & Station == 4, T, F),
         bge = ifelse(Cruise == "AT34" & Station %in% c(1, 2, 3), T, bge),
         bge = ifelse(Cruise == "AT38" & Station %in% c(3, 6), T, bge),
         degree_bin = as.character(degree_bin))  %>% 
  distinct()  %>% 
  #for 3 stations, 
  filter(int_delta_DOC_ez > 0)

bcd <- read_rds("~/naames_bioav_ms/Output/processed_integrated_BCD.rds")
```

Units for imported data frames are currently:

  - BP, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - BCD, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - BA, cells m<sup>-3</sup>
  - BC, µmol C m<sup>-3</sup>
  - mew, d<sup>-1</sup>
  - NPP, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - ∆DOC (from export MS, integrated to Ez), mol C
m<sup>-2</sup>

# Bar plots: NPP, BP, BA, µ, ∆DOC

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

Error bars for µ represent standard deviation from mean of values
calculated using different CCFs to convert BA to BC (Global Initial CCF,
Global Stationary CCF, Season-Specific Initial CCF, Season-Specific
Stationary CCF). All other error bars represent standard deviation from
mean of station values.

# Table: NPP and Bacterial Parameters

``` r
bcd_table <- bcd %>% 
  select(Cruise:degree_bin, int.bcd.global_bge, int.bcd.season_bge) %>% 
  distinct() %>% 
  mutate_at(vars(int.bcd.global_bge, int.bcd.season_bge), funs(./1000)) %>% 
  pivot_longer(cols = c(int.bcd.global_bge, int.bcd.season_bge), names_to = "var", values_to = "bcd") %>% 
  select(-var) %>% 
  distinct() %>% 
  group_by(Cruise, Season, Subregion, degree_bin, Station) %>% 
  summarise(ave_BCD = round(mean(bcd), 2),
            sd_BCD = round(sd(bcd), 2)) %>% 
  ungroup() %>% 
  left_join(., bcd %>% 
              select(Cruise:Station, degree_bin,ave_Ez, sd_Ez, int.NPP, sd_int.NPP) %>% 
              distinct() %>% 
              mutate_at(vars(contains("NPP")), funs(round(./10^3, 2)))
            ) %>% 
  left_join(., bcd %>% 
              select(Cruise:Station, degree_bin, bcd.npp.global_bge, bcd.npp.season_bge) %>% 
              distinct() %>% 
              pivot_longer(cols = c( bcd.npp.global_bge, bcd.npp.season_bge), names_to = "var", values_to = "bcd_npp") %>% 
              select(-var) %>% 
              distinct() %>% 
              group_by(Cruise, Season, Subregion, degree_bin, Station) %>% 
              summarise(ave_BCD_NPP = mean(bcd_npp),
                        sd_BCD_NPP = sd(bcd_npp)) %>% 
              ungroup()
            ) %>% 
  arrange(factor(Season, levels = levels)) %>% 
  group_by(Cruise) %>% 
  mutate(Cruise_ave_Ez = mean(ave_Ez),
         Cruise_sd_Ez = sd(ave_Ez),
         Cruise_ave_NPP = mean(int.NPP),
         Cruise_sd_NPP = sd(int.NPP),
         Cruise_ave_BCD = mean(ave_BCD),
         Cruise_sd_BCD = sd(ave_BCD),
         Cruise_ave_BCD_NPP = mean(ave_BCD_NPP),
         Cruise_sd_BCD_NPP = sd(ave_BCD_NPP)
         ) %>% 
  ungroup() %>% 
  group_by(Cruise, Subregion) %>% 
  mutate(Cruise_SR_ave_NPP = mean(int.NPP),
         Cruise_SR_sd_NPP = mean(int.NPP),
         Cruise_SR_ave_BCD = mean(ave_BCD),
         Cruise_SR_sd_BCD = sd(ave_BCD),
         Cruise_SR_ave_BCD_NPP = mean(ave_BCD_NPP),
         Cruise_SR_sd_BCD_NPP = sd(ave_BCD_NPP)
         ) %>% 
  ungroup() %>% 
  mutate_at(vars(contains("BCD_NPP")), round) %>% 
  mutate_at(vars(Cruise_ave_NPP:Cruise_sd_BCD, Cruise_SR_ave_NPP:Cruise_SR_sd_BCD ), round, 2) %>% 
  arrange(factor(Season, levels = levels), factor(Subregion, levels = levels))
```

| Cruise | Season       | Subregion   | degree\_bin | Station | ave\_BCD | sd\_BCD | ave\_Ez | sd\_Ez | int.NPP | sd\_int.NPP | ave\_BCD\_NPP | sd\_BCD\_NPP | Cruise\_ave\_Ez | Cruise\_sd\_Ez | Cruise\_ave\_NPP | Cruise\_sd\_NPP | Cruise\_ave\_BCD | Cruise\_sd\_BCD | Cruise\_ave\_BCD\_NPP | Cruise\_sd\_BCD\_NPP | Cruise\_SR\_ave\_NPP | Cruise\_SR\_sd\_NPP | Cruise\_SR\_ave\_BCD | Cruise\_SR\_sd\_BCD | Cruise\_SR\_ave\_BCD\_NPP | Cruise\_SR\_sd\_BCD\_NPP |
| :----- | :----------- | :---------- | ----------: | ------: | -------: | ------: | ------: | -----: | ------: | ----------: | ------------: | -----------: | --------------: | -------------: | ---------------: | --------------: | ---------------: | --------------: | --------------------: | -------------------: | -------------------: | ------------------: | -------------------: | ------------------: | ------------------------: | -----------------------: |
| AT39   | Early Spring | GS/Sargasso |          39 |       1 |     0.12 |    0.01 |     106 |     NA |    0.60 |          NA |            19 |            1 |       112.50000 |       12.79323 |             0.72 |            0.40 |             0.14 |            0.11 |                    19 |                    6 |                 0.96 |                0.96 |                 0.22 |                0.13 |                        22 |                        3 |
| AT39   | Early Spring | GS/Sargasso |          39 |       2 |     0.31 |    0.02 |      98 |     NA |    1.31 |          NA |            24 |            1 |       112.50000 |       12.79323 |             0.72 |            0.40 |             0.14 |            0.11 |                    19 |                    6 |                 0.96 |                0.96 |                 0.22 |                0.13 |                        22 |                        3 |
| AT39   | Early Spring | Subtropical |          44 |       3 |     0.09 |    0.01 |     120 |     NA |    0.41 |          NA |            22 |            1 |       112.50000 |       12.79323 |             0.72 |            0.40 |             0.14 |            0.11 |                    19 |                    6 |                 0.48 |                0.48 |                 0.07 |                0.02 |                        17 |                        8 |
| AT39   | Early Spring | Subtropical |          44 |       4 |     0.06 |    0.00 |     126 |     NA |    0.56 |          NA |            11 |            1 |       112.50000 |       12.79323 |             0.72 |            0.40 |             0.14 |            0.11 |                    19 |                    6 |                 0.48 |                0.48 |                 0.07 |                0.02 |                        17 |                        8 |
| AT34   | Late Spring  | Subtropical |          44 |       5 |     0.53 |    0.02 |      91 |     NA |    0.99 |          NA |            53 |            2 |        79.33333 |       23.64459 |             1.89 |            1.12 |             0.40 |            0.17 |                    25 |                   14 |                 0.90 |                0.90 |                 0.35 |                0.25 |                        37 |                       23 |
| AT34   | Late Spring  | Subtropical |          48 |       4 |     0.17 |    0.00 |     116 |     33 |    0.80 |        0.08 |            21 |            1 |        79.33333 |       23.64459 |             1.89 |            1.12 |             0.40 |            0.17 |                    25 |                   14 |                 0.90 |                0.90 |                 0.35 |                0.25 |                        37 |                       23 |
| AT34   | Late Spring  | Temperate   |          50 |       3 |     0.48 |    0.01 |      52 |      1 |    3.56 |        0.32 |            13 |            0 |        79.33333 |       23.64459 |             1.89 |            1.12 |             0.40 |            0.17 |                    25 |                   14 |                 3.56 |                3.56 |                 0.48 |                  NA |                        13 |                       NA |
| AT34   | Late Spring  | Subpolar    |          54 |       0 |     0.24 |    0.01 |      87 |     NA |    1.47 |          NA |            17 |            0 |        79.33333 |       23.64459 |             1.89 |            1.12 |             0.40 |            0.17 |                    25 |                   14 |                 1.99 |                1.99 |                 0.42 |                0.18 |                        21 |                        5 |
| AT34   | Late Spring  | Subpolar    |          54 |       2 |     0.60 |    0.02 |      58 |      5 |    2.98 |        0.45 |            20 |            1 |        79.33333 |       23.64459 |             1.89 |            1.12 |             0.40 |            0.17 |                    25 |                   14 |                 1.99 |                1.99 |                 0.42 |                0.18 |                        21 |                        5 |
| AT34   | Late Spring  | Subpolar    |          56 |       1 |     0.41 |    0.01 |      72 |      1 |    1.53 |        0.24 |            27 |            1 |        79.33333 |       23.64459 |             1.89 |            1.12 |             0.40 |            0.17 |                    25 |                   14 |                 1.99 |                1.99 |                 0.42 |                0.18 |                        21 |                        5 |
| AT38   | Early Autumn | GS/Sargasso |          42 |       1 |     0.12 |    0.01 |     236 |     11 |    0.14 |        0.03 |            84 |            5 |       179.16667 |       47.04643 |             0.30 |            0.20 |             0.11 |            0.03 |                    48 |                   24 |                 0.14 |                0.14 |                 0.12 |                  NA |                        84 |                       NA |
| AT38   | Early Autumn | Subtropical |          44 |       2 |     0.14 |    0.01 |     207 |     NA |    0.23 |          NA |            62 |            4 |       179.16667 |       47.04643 |             0.30 |            0.20 |             0.11 |            0.03 |                    48 |                   24 |                 0.19 |                0.19 |                 0.10 |                0.04 |                        53 |                       14 |
| AT38   | Early Autumn | Subtropical |          47 |       3 |     0.07 |    0.00 |     200 |     NA |    0.18 |          NA |            37 |            2 |       179.16667 |       47.04643 |             0.30 |            0.20 |             0.11 |            0.03 |                    48 |                   24 |                 0.19 |                0.19 |                 0.10 |                0.04 |                        53 |                       14 |
| AT38   | Early Autumn | Subtropical |          49 |       4 |     0.09 |    0.01 |     174 |     21 |    0.16 |        0.07 |            60 |            4 |       179.16667 |       47.04643 |             0.30 |            0.20 |             0.11 |            0.03 |                    48 |                   24 |                 0.19 |                0.19 |                 0.10 |                0.04 |                        53 |                       14 |
| AT38   | Early Autumn | Temperate   |          52 |       5 |     0.11 |    0.01 |     157 |     NA |    0.42 |          NA |            27 |            2 |       179.16667 |       47.04643 |             0.30 |            0.20 |             0.11 |            0.03 |                    48 |                   24 |                 0.42 |                0.42 |                 0.11 |                  NA |                        27 |                       NA |
| AT38   | Early Autumn | Subpolar    |          53 |       6 |     0.13 |    0.01 |     101 |      7 |    0.65 |        0.11 |            20 |            1 |       179.16667 |       47.04643 |             0.30 |            0.20 |             0.11 |            0.03 |                    48 |                   24 |                 0.65 |                0.65 |                 0.13 |                  NA |                        20 |                       NA |
| AT32   | Late Autumn  | Subtropical |          40 |       7 |     0.15 |    0.01 |     112 |     19 |    0.32 |        0.13 |            48 |            3 |       101.00000 |       22.46775 |             0.19 |            0.07 |             0.13 |            0.01 |                    76 |                   21 |                 0.21 |                0.21 |                 0.14 |                0.01 |                        70 |                       16 |
| AT32   | Late Autumn  | Subtropical |          43 |       6 |     0.14 |    0.01 |     104 |      1 |    0.18 |        0.01 |            78 |            5 |       101.00000 |       22.46775 |             0.19 |            0.07 |             0.13 |            0.01 |                    76 |                   21 |                 0.21 |                0.21 |                 0.14 |                0.01 |                        70 |                       16 |
| AT32   | Late Autumn  | Subtropical |          44 |       5 |     0.14 |    0.01 |     103 |     NA |    0.17 |          NA |            85 |            5 |       101.00000 |       22.46775 |             0.19 |            0.07 |             0.13 |            0.01 |                    76 |                   21 |                 0.21 |                0.21 |                 0.14 |                0.01 |                        70 |                       16 |
| AT32   | Late Autumn  | Subtropical |          46 |       4 |     0.12 |    0.01 |     126 |     NA |    0.17 |          NA |            70 |            4 |       101.00000 |       22.46775 |             0.19 |            0.07 |             0.13 |            0.01 |                    76 |                   21 |                 0.21 |                0.21 |                 0.14 |                0.01 |                        70 |                       16 |
| AT32   | Late Autumn  | Temperate   |          51 |       3 |     0.12 |    0.01 |      59 |     NA |    0.11 |          NA |           109 |            7 |       101.00000 |       22.46775 |             0.19 |            0.07 |             0.13 |            0.01 |                    76 |                   21 |                 0.11 |                0.11 |                 0.12 |                  NA |                       109 |                       NA |
| AT32   | Late Autumn  | Subpolar    |          54 |       2 |     0.12 |    0.01 |     102 |      4 |    0.19 |        0.02 |            64 |            4 |       101.00000 |       22.46775 |             0.19 |            0.07 |             0.13 |            0.01 |                    76 |                   21 |                 0.19 |                0.19 |                 0.12 |                  NA |                        64 |                       NA |

BCD &
NPP

# Line plots: DOC Decay Curves

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

Black vertical dashed and dotted lines indicate the 7 and 30-day marks,
respectively. Dashed decay lines indicate experiments in which BGEs
could be calculated.

# Table: ∆DOC, Bioavailability, Persistence

``` r
bioav.table <- doc %>% 
  select(Season:Subregion, degree_bin, int_delta_DOC_ez, contains("bioav_accm_doc"), contains("bioav_doc"), contains("per_bioav"), contains("persis"), contains("ddoc")) %>% 
  distinct() %>% 
  group_by(Season, Subregion, degree_bin) %>% 
  summarise(ave_accm = round(mean(int_delta_DOC_ez, na.rm = T), 2),
            sd_accm = round(sd(int_delta_DOC_ez, na.rm = T), 2),
            
            ave_total.deldoc = round(mean(total.bioav_accm_doc, na.rm = T), 2),
            sd_total.deldoc = round(sd(total.bioav_accm_doc, na.rm = T), 2),
            ave_st.deldoc = round(mean(st.bioav_accm_doc, na.rm = T), 2),
            sd_st.deldoc = round(sd(st.bioav_accm_doc, na.rm = T), 2),
            ave_total.per_bioav = round(mean(total.per_bioav, na.rm = T), 2),
            sd_total.per_bioav = round(sd(total.per_bioav, na.rm = T), 2),
            ave_st.per_bioav = round(mean(st.per_bioav, na.rm = T), 2),
            sd_st.per_bioav = round(sd(st.per_bioav, na.rm = T), 2),
            
            ave_per_persis = round(mean(per_persis, na.rm = T), 2),
            sd_per_persis = round(sd(per_persis, na.rm = T), 2),
            
            ave_total.ddoc = round(mean(total.ddoc, na.rm = T), 2),
            sd_total.ddoc = round(sd(total.ddoc, na.rm = T), 2),
            ave_st.ddoc = round(mean(st.ddoc, na.rm = T), 2),
            sd_st.ddoc = round(sd(st.ddoc, na.rm = T), 2),
            
            
            ) %>%
   arrange(factor(Season, levels = levels), degree_bin) %>% 
  ungroup()
```

| Season       | Subregion   | degree\_bin | ave\_accm | sd\_accm | ave\_total.deldoc | sd\_total.deldoc | ave\_st.deldoc | sd\_st.deldoc | ave\_total.per\_bioav | sd\_total.per\_bioav | ave\_st.per\_bioav | sd\_st.per\_bioav | ave\_per\_persis | sd\_per\_persis | ave\_total.ddoc | sd\_total.ddoc | ave\_st.ddoc | sd\_st.ddoc |
| :----------- | :---------- | :---------- | --------: | -------: | ----------------: | ---------------: | -------------: | ------------: | --------------------: | -------------------: | -----------------: | ----------------: | ---------------: | --------------: | --------------: | -------------: | -----------: | ----------: |
| Early Spring | GS/Sargasso | 39          |      0.18 |     0.04 |              3.90 |             1.13 |           1.15 |          1.20 |                 100.5 |                16.26 |               28.5 |             27.58 |            \-0.5 |           16.26 |            46.5 |          13.44 |        164.5 |      171.83 |
| Early Spring | Subtropical | 44          |      0.16 |     0.22 |              2.45 |             0.92 |           0.35 |          0.49 |                  70.5 |                 2.12 |               13.5 |             19.09 |             29.5 |            2.12 |            28.0 |          14.14 |         50.0 |       70.71 |
| Late Spring  | Subtropical | 44          |      0.32 |       NA |              4.00 |               NA |           2.70 |            NA |                  34.0 |                   NA |               23.0 |                NA |             66.0 |              NA |            38.0 |             NA |        386.0 |          NA |
| Late Spring  | Subtropical | 48          |      0.24 |       NA |               NaN |               NA |           0.70 |            NA |                   NaN |                   NA |               25.0 |                NA |              NaN |              NA |             NaN |             NA |        100.0 |          NA |
| Late Spring  | Temperate   | 50          |      0.20 |       NA |               NaN |               NA |           2.30 |            NA |                   NaN |                   NA |               29.0 |                NA |              NaN |              NA |             NaN |             NA |        329.0 |          NA |
| Late Spring  | Subpolar    | 54          |      0.17 |     0.01 |               NaN |               NA |           1.00 |            NA |                   NaN |                   NA |               14.0 |                NA |              NaN |              NA |             NaN |             NA |        143.0 |          NA |
| Late Spring  | Subpolar    | 56          |      0.03 |       NA |               NaN |               NA |           2.50 |            NA |                   NaN |                   NA |               76.0 |                NA |              NaN |              NA |             NaN |             NA |        357.0 |          NA |
| Early Autumn | GS/Sargasso | 42          |      0.43 |       NA |              4.70 |               NA |           2.30 |            NA |                  39.0 |                   NA |               19.0 |                NA |             61.0 |              NA |            68.0 |             NA |        329.0 |          NA |
| Early Autumn | Subtropical | 44          |      0.62 |     0.39 |              5.90 |               NA |           1.60 |            NA |                  31.0 |                   NA |                8.0 |                NA |             69.0 |              NA |            87.0 |             NA |        229.0 |          NA |
| Early Autumn | Subtropical | 47          |      0.22 |       NA |              4.50 |               NA |           1.60 |            NA |                  32.0 |                   NA |               11.0 |                NA |             68.0 |              NA |            68.0 |             NA |        229.0 |          NA |
| Early Autumn | Subtropical | 49          |      0.97 |       NA |              3.80 |               NA |           0.90 |            NA |                  25.0 |                   NA |                6.0 |                NA |             75.0 |              NA |            59.0 |             NA |        129.0 |          NA |
| Early Autumn | Temperate   | 50          |      0.34 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Early Autumn | Temperate   | 52          |      0.76 |       NA |              5.30 |               NA |           0.90 |            NA |                  37.0 |                   NA |                6.0 |                NA |             63.0 |              NA |            85.0 |             NA |        129.0 |          NA |
| Early Autumn | Subpolar    | 53          |      0.39 |     0.00 |              5.00 |               NA |           1.80 |            NA |                  56.0 |                   NA |               20.0 |                NA |             44.0 |              NA |            85.0 |             NA |        257.0 |          NA |
| Late Autumn  | Subtropical | 40          |      0.72 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Late Autumn  | Subtropical | 41          |      0.65 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Late Autumn  | Subtropical | 43          |      0.67 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Late Autumn  | Subtropical | 44          |      0.57 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Late Autumn  | Subtropical | 46          |      0.39 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Late Autumn  | Temperate   | 51          |      0.36 |     0.08 |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |
| Late Autumn  | Subpolar    | 54          |      0.16 |       NA |               NaN |               NA |            NaN |            NA |                   NaN |                   NA |                NaN |                NA |              NaN |              NA |             NaN |             NA |          NaN |          NA |

Seasonal Accumulated DOC Bioavailability and
Persistance

# Bar plots: Experiment ∆DOC and %Bioavailability

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

# Bar plots: BCD and BCD:NPP

We’ll convert BCD and NPP to mmol C m<sup>-3</sup> d<sup>-1</sup> before
plotting.

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />
<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

# Merge DOC and BCD data

Here we merge the integrated and depth-normalized data. We’ll also
convert:

  - BCD to mmol C m<sup>-3</sup> d<sup>-1</sup>
  - BC to mmol C m<sup>-3</sup>
  - NPP to mmol C m<sup>-3</sup> d<sup>-1</sup>

<!-- end list -->

``` r
bcd_bioav <- bcd %>% 
  select(Cruise:degree_bin, int.bcd.season_bge, int.bcd.global_bge) %>% distinct() %>% 
  pivot_longer(cols = c(int.bcd.season_bge, int.bcd.global_bge), names_to = "var", values_to = "bcd") %>% 
  select(-var) %>% distinct() %>% 
  group_by(Cruise,Subregion, degree_bin, Station) %>% 
  summarise(ave_BCD = mean(bcd),
            sd_BCD = sd(bcd)) %>% ungroup() %>% 
  left_join(., bcd %>% 
              select(Cruise:degree_bin, int.bc.global_i:int.bc.season_s) %>% distinct() %>% 
              pivot_longer(cols = c(int.bc.global_i:int.bc.season_s), names_to = "var", values_to = "bc") %>% 
              select(-var) %>% distinct() %>% 
              group_by(Cruise,Subregion, degree_bin, Station) %>% 
              summarise(ave_BC = mean(bc),
                        sd_BC = sd(bc)) %>% ungroup()
  ) %>% 
  left_join(., bcd %>% 
              select(Cruise:degree_bin, mew.global_i:mew.season_s) %>% distinct() %>%
              pivot_longer(cols = c(mew.global_i:mew.season_s), names_to = "var", values_to = "mew") %>% 
              select(-var) %>% distinct() %>% 
              group_by(Cruise,Subregion, degree_bin, Station) %>%
              summarise(ave_mew = mean(mew),
                        sd_mew = sd(mew)) %>% 
              ungroup()
  ) %>% 
  left_join(., bcd %>% 
              select(Cruise:degree_bin, bcd.npp.global_bge, bcd.npp.season_bge) %>%
              distinct() %>%
              pivot_longer(cols = c(bcd.npp.global_bge, bcd.npp.season_bge), names_to = "var", values_to = "bcd.npp") %>% 
              select(-var) %>% distinct() %>% 
              group_by(Cruise,Subregion, degree_bin, Station) %>%
              summarise(ave_bcd.npp = mean(bcd.npp),
                        sd_bcd.npp = sd(bcd.npp)) %>% 
              ungroup()
  ) %>% 
  left_join(bcd %>% 
              select(Cruise:sd_Ez, int.bp:sd_int.ba, int.NPP:sd_int.NPP) %>%
              distinct(), 
            .) %>% 
  mutate(Station = as.character(Station)) %>% 
  left_join(., doc %>% 
              select(Cruise:Subregion, degree_bin, int_delta_DOC_ez, accm_doc:bge) %>%
              distinct() %>% 
              mutate(degree_bin = as.numeric(degree_bin))
            ) %>% 
  select(Cruise:ave_Ez, everything()) %>% 
  mutate_at(vars(int.NPP, int.NPP:sd_BC), funs(./1000))
```

# Regressions: Property-Property

## \*\* BCD v NPP

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_BCD ~ int.NPP, data = bcd_bioav, nperm =
    ## 99)
    ## 
    ## n = 22   r = 0.8432318   r-square = 0.7110398 
    ## Parametric P-values:   2-tailed = 8.331438e-07    1-tailed = 4.165719e-07 
    ## Angle between the two OLS regression lines = 3.257019 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method  Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.09139499 0.1441207        8.201040              0.01
    ## 2     MA 0.09045965 0.1453222        8.268469              0.01
    ## 3    SMA 0.07053676 0.1709147        9.698976                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS     0.04034571      0.14244427  0.1012668   0.1869746
    ## 2     MA     0.05658255      0.12392227  0.1023369   0.1888400
    ## 3    SMA     0.03305794      0.09977855  0.1333514   0.2190592
    ## 
    ## Eigenvalues: 0.8695995 0.007042267 
    ## 
    ## H statistic used for computing C.I. of MA: 0.001790768

## BCD v Seasonally Accumulated DOC

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_BCD ~ int_delta_DOC_ez, data = bcd_bioav,
    ## nperm = 99)
    ## 
    ## n = 22   r = -0.3622026   r-square = 0.1311907 
    ## Parametric P-values:   2-tailed = 0.09761773    1-tailed = 0.04880887 
    ## Angle between the two OLS regression lines = 46.11499 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept      Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.2838786 -0.2095363       -11.83433              0.03
    ## 2     MA 0.3142083 -0.2886885       -16.10282              0.03
    ## 3    SMA 0.4252611 -0.5785059       -30.04964                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS      0.1665008       0.4012565 -0.4610501  0.04197743
    ## 2     MA      0.1865100       0.4694859 -0.6939206  0.04456922
    ## 3    SMA      0.3489297       0.5416807 -0.8823293 -0.37930181
    ## 
    ## Eigenvalues: 0.07884381 0.02038421 
    ## 
    ## H statistic used for computing C.I. of MA: 0.1023136

## BCD v % DOC Bioavailability

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_BCD ~ total.per_bioav, data = bcd_bioav,
    ## nperm = 99)
    ## 
    ## n = 11   r = -0.03618315   r-square = 0.00130922 
    ## Parametric P-values:   2-tailed = 0.9158866    1-tailed = 0.4579433 
    ## Angle between the two OLS regression lines = 7.711631 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept         Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.1703331 -0.0001775209     -0.01017120              0.45
    ## 2     MA 0.1703333 -0.0001775251     -0.01017144              0.45
    ## 3    SMA 0.4265402 -0.0049061748     -0.28110085                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept   2.5%-Slope  97.5%-Slope
    ## 1    OLS    -0.05317205       0.3938383 -0.003874611  0.003519569
    ## 2     MA    -0.02998780       0.3706547 -0.003874732  0.003519677
    ## 3    SMA     0.29324991       0.6938801 -0.009840300 -0.002446120
    ## 
    ## Eigenvalues: 790.9637 0.019014 
    ## 
    ## H statistic used for computing C.I. of MA: 1.366913e-05

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_BCD ~ st.per_bioav, data = bcd_bioav, nperm
    ## = 99)
    ## 
    ## n = 15   r = 0.2848781   r-square = 0.08115552 
    ## Parametric P-values:   2-tailed = 0.3034158    1-tailed = 0.1517079 
    ## Angle between the two OLS regression lines = 1.772885 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method  Intercept       Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.17012379 0.002734089       0.1566514              0.16
    ## 2     MA 0.17011884 0.002734320       0.1566646              0.16
    ## 3    SMA 0.02324893 0.009597401       0.5498737                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept   2.5%-Slope 97.5%-Slope
    ## 1    OLS     0.01357132       0.3266763 -0.002778185 0.008246363
    ## 2     MA     0.05214262       0.2880915 -0.002778421 0.008247228
    ## 3    SMA    -0.12617939       0.1097460  0.005555484 0.016580032
    ## 
    ## Eigenvalues: 373.5456 0.03161451 
    ## 
    ## H statistic used for computing C.I. of MA: 3.038986e-05

## BCD v DOC removal rate

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_BCD ~ total.ddoc, data = bcd_bioav, nperm =
    ## 99)
    ## 
    ## n = 11   r = -0.2906552   r-square = 0.08448044 
    ## Parametric P-values:   2-tailed = 0.3858909    1-tailed = 0.1929454 
    ## Angle between the two OLS regression lines = 1.076953 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept        Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.2614861 -0.001734719     -0.09939199              0.24
    ## 2     MA 0.2614894 -0.001734776     -0.09939523              0.24
    ## 3    SMA 0.5074190 -0.005968306     -0.34195471                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept   2.5%-Slope  97.5%-Slope
    ## 1    OLS   -0.006064326       0.5290366 -0.006040842  0.002571404
    ## 2     MA    0.011334183       0.5116484 -0.006041112  0.002571496
    ## 3    SMA    0.338092626       0.8383858 -0.011665701 -0.003053454
    ## 
    ## Eigenvalues: 534.4925 0.01743046 
    ## 
    ## H statistic used for computing C.I. of MA: 1.854379e-05

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_BCD ~ st.ddoc, data = bcd_bioav, nperm =
    ## 99)
    ## 
    ## n = 15   r = 0.35897   r-square = 0.1288595 
    ## Parametric P-values:   2-tailed = 0.1888503    1-tailed = 0.09442513 
    ## Angle between the two OLS regression lines = 0.2143474 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method  Intercept        Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  0.1162594 0.0005533841      0.03170657              0.14
    ## 2     MA  0.1162592 0.0005533853      0.03170664              0.14
    ## 3    SMA -0.0844120 0.0015415888      0.08832646                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept    2.5%-Slope 97.5%-Slope
    ## 1    OLS    -0.08546416      0.31798301 -0.0003087381 0.001415506
    ## 2     MA    -0.05880961      0.29132783 -0.0003087387 0.001415510
    ## 3    SMA    -0.30510792      0.04502866  0.0009041594 0.002628404
    ## 
    ## Eigenvalues: 14478.07 0.02997338 
    ## 
    ## H statistic used for computing C.I. of MA: 7.432574e-07

## \*\* BCD:NPP v Seasonally Accumulated DOC

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_bcd.npp ~ int_delta_DOC_ez, data =
    ## bcd_bioav, nperm = 99)
    ## 
    ## n = 22   r = 0.4907742   r-square = 0.2408593 
    ## Parametric P-values:   2-tailed = 0.02038641    1-tailed = 0.01019321 
    ## Angle between the two OLS regression lines = 0.8488135 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  24.57823  51.23432        88.88183              0.01
    ## 2     MA -37.29245 212.69992        89.73063              0.01
    ## 3    SMA   4.20806 104.39490        89.45118                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS       4.778622      44.3778319    8.80832    93.66033
    ## 2     MA    -429.832680      -0.3730376  116.35032  1237.12281
    ## 3    SMA     -15.226027      17.2877194   70.26056   155.11256
    ## 
    ## Eigenvalues: 810.2684 0.05643823 
    ## 
    ## H statistic used for computing C.I. of MA: 1.515613e-05

## \*\* BCD:NPP v % DOC Bioavailability

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_bcd.npp ~ total.per_bioav, data =
    ## bcd_bioav, nperm = 99)
    ## 
    ## n = 11   r = -0.6558056   r-square = 0.430081 
    ## Parametric P-values:   2-tailed = 0.02845251    1-tailed = 0.01422625 
    ## Angle between the two OLS regression lines = 23.09672 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept      Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  67.41761 -0.5401644       -28.37634              0.01
    ## 2     MA  78.55377 -0.7456976       -36.71180              0.01
    ## 3    SMA  82.77821 -0.8236654       -39.47710                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS       39.07196        95.76326  -1.009041  -0.0712874
    ## 2     MA       49.79679       137.08523  -1.825976  -0.2149479
    ## 3    SMA       64.09788       114.90710  -1.416648  -0.4788942
    ## 
    ## Eigenvalues: 1109.563 218.0097 
    ## 
    ## H statistic used for computing C.I. of MA: 0.1730358

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_bcd.npp ~ st.per_bioav, data = bcd_bioav,
    ## nperm = 99)
    ## 
    ## n = 15   r = -0.2446094   r-square = 0.05983374 
    ## Parametric P-values:   2-tailed = 0.3795976    1-tailed = 0.1897988 
    ## Angle between the two OLS regression lines = 62.39348 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Confidence interval = NA when the limits of the confidence interval
    ## cannot be computed. This happens when the correlation is 0
    ## or the C.I. includes all 360 deg. of the plane (H >= 1)
    ## 
    ## Regression results
    ##   Method Intercept      Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  39.18973 -0.2702043       -15.12049              0.19
    ## 2     MA  65.23668 -1.4873517       -56.08564              0.19
    ## 3    SMA  57.04656 -1.1046360       -47.84623                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS       20.96308        57.41637 -0.9119718   0.3715632
    ## 2     MA             NA              NA         NA          NA
    ## 3    SMA       47.01270        74.48035 -1.9192988  -0.6357638
    ## 
    ## Eigenvalues: 523.6656 305.682 
    ## 
    ## H statistic used for computing C.I. of MA: 1.209453

## BCD:NPP v % DOC removal rate

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_bcd.npp ~ total.ddoc, data = bcd_bioav,
    ## nperm = 99)
    ## 
    ## n = 11   r = 0.3352604   r-square = 0.1123995 
    ## Parametric P-values:   2-tailed = 0.3135141    1-tailed = 0.156757 
    ## Angle between the two OLS regression lines = 52.93142 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Confidence interval = NA when the limits of the confidence interval
    ## cannot be computed. This happens when the correlation is 0
    ## or the C.I. includes all 360 deg. of the plane (H >= 1)
    ## 
    ## Regression results
    ##   Method Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  18.63639 0.3359241        18.56844              0.18
    ## 2     MA -20.28408 1.0059164        45.16899              0.18
    ## 3    SMA -20.05539 1.0019797        45.05666                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS      -25.59071       62.863482 -0.3758943    1.047742
    ## 2     MA             NA              NA         NA          NA
    ## 3    SMA      -74.59829        8.102065  0.5172661    1.940903
    ## 
    ## Eigenvalues: 715.1016 355.9986 
    ## 
    ## H statistic used for computing C.I. of MA: 1.122486

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_bcd.npp ~ st.ddoc, data = bcd_bioav, nperm
    ## = 99)
    ## 
    ## n = 15   r = 0.331706   r-square = 0.1100289 
    ## Parametric P-values:   2-tailed = 0.2271205    1-tailed = 0.1135602 
    ## Angle between the two OLS regression lines = 24.77448 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept      Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 21.455750 0.05885556        3.368290              0.13
    ## 2     MA 21.112490 0.06054594        3.464797              0.13
    ## 3    SMA -2.623353 0.17743289       10.061442                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept  2.5%-Slope 97.5%-Slope
    ## 1    OLS    -2.01171128        44.92321 -0.04143921   0.1591503
    ## 2     MA    -0.08319039        42.04457 -0.04253390   0.1649239
    ## 3    SMA   -28.34766166        12.38539  0.10352248   0.3041120
    ## 
    ## Eigenvalues: 14529.66 404.2127 
    ## 
    ## H statistic used for computing C.I. of MA: 0.01056753

## \*\* µ v NPP

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_mew ~ int.NPP, data = bcd_bioav, nperm =
    ## 99)
    ## 
    ## n = 22   r = 0.4559422   r-square = 0.2078833 
    ## Parametric P-values:   2-tailed = 0.03295421    1-tailed = 0.0164771 
    ## Angle between the two OLS regression lines = 1.710908 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method  Intercept       Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.02399992 0.007841369       0.4492681              0.03
    ## 2     MA 0.02399849 0.007843206       0.4493734              0.03
    ## 3    SMA 0.01671597 0.017198164       0.9852851                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept   2.5%-Slope 97.5%-Slope
    ## 1    OLS     0.01549505      0.03250479 0.0007018590  0.01498088
    ## 2     MA     0.01843886      0.02955749 0.0007022413  0.01498497
    ## 3    SMA     0.01005031      0.02116604 0.0114817002  0.02576072
    ## 
    ## Eigenvalues: 0.8518126 0.0001995465 
    ## 
    ## H statistic used for computing C.I. of MA: 5.099022e-05

## µ v Seasonally Accumulated DOC

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_mew ~ int_delta_DOC_ez, data = bcd_bioav,
    ## nperm = 99)
    ## 
    ## n = 22   r = -0.2687469   r-square = 0.07222492 
    ## Parametric P-values:   2-tailed = 0.2265249    1-tailed = 0.1132624 
    ## Angle between the two OLS regression lines = 11.32543 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method  Intercept       Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.03609874 -0.01564422      -0.8962747              0.13
    ## 2     MA 0.03611764 -0.01569355      -0.8991002              0.13
    ## 3    SMA 0.05240983 -0.05821172      -3.3315260                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept  2.5%-Slope 97.5%-Slope
    ## 1    OLS     0.02389346      0.04830402 -0.04179734  0.01050890
    ## 2     MA     0.02606533      0.04617824 -0.04194897  0.01054026
    ## 3    SMA     0.04453622      0.06457901 -0.08996996 -0.03766373
    ## 
    ## Eigenvalues: 0.07436479 0.0002336779 
    ## 
    ## H statistic used for computing C.I. of MA: 0.0006879666

## µ v % DOC Bioavailability

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_mew ~ total.per_bioav, data = bcd_bioav,
    ## nperm = 99)
    ## 
    ## n = 11   r = 0.2926405   r-square = 0.08563844 
    ## Parametric P-values:   2-tailed = 0.3825046    1-tailed = 0.1912523 
    ## Angle between the two OLS regression lines = 0.1313827 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method   Intercept        Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  0.01401907 0.0002147669      0.01230524              0.23
    ## 2     MA  0.01401907 0.0002147670      0.01230524              0.23
    ## 3    SMA -0.01410815 0.0007338935      0.04204899                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept    2.5%-Slope  97.5%-Slope
    ## 1    OLS    -0.01797142     0.046009574 -0.0003144010 0.0007439349
    ## 2     MA    -0.01465223     0.042690366 -0.0003144012 0.0007439354
    ## 3    SMA    -0.05203812     0.005304448  0.0003756074 0.0014339433
    ## 
    ## Eigenvalues: 790.9637 0.0003895297 
    ## 
    ## H statistic used for computing C.I. of MA: 2.80019e-07

    ## RMA was not requested: it will not be computed.

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ave_mew ~ st.per_bioav, data = bcd_bioav, nperm
    ## = 99)
    ## 
    ## n = 15   r = 0.2282662   r-square = 0.05210544 
    ## Parametric P-values:   2-tailed = 0.4132078    1-tailed = 0.2066039 
    ## Angle between the two OLS regression lines = 0.2320426 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method   Intercept        Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS 0.023974921 0.0002226234      0.01275538              0.23
    ## 2     MA 0.023974917 0.0002226236      0.01275540              0.23
    ## 3    SMA 0.007868068 0.0009752801      0.05587942                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept    2.5%-Slope  97.5%-Slope
    ## 1    OLS    0.007816658      0.04013319 -0.0003463154 0.0007915623
    ## 2     MA    0.011799611      0.03615022 -0.0003463157 0.0007915632
    ## 3    SMA   -0.007598929      0.01675165  0.0005601594 0.0016980370
    ## 
    ## Eigenvalues: 373.5429 0.00033679 
    ## 
    ## H statistic used for computing C.I. of MA: 3.23692e-07

# Plots: Property-Property

<img src="BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-46-1.png" style="display: block; margin: auto;" />
