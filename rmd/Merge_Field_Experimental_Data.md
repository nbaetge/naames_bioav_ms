Merge Field and Experimental Data
================
Nicholas Baetge
5/21/2020

# Intro

This document shows how processed **individual bottle** data from NAAMES
remineralization bioassays were combined with field data to calculate
derived variables such as DOC bioavailability and BCD.

``` r
library(tidyverse) 
library(lubridate)
library(readxl)
library(zoo)
library(oce)
```

# Merge Field NPP, Bacterial Productivity and Abundance Data

We’ll convert NPP from mg C m<sup>-2</sup> d to µmol C m<sup>-2</sup>
d<sup>-1</sup>

``` r
#NPP data
#we'll convert NPP from mg C/ m2 / d to umol C/ m2 / d 
npp <- read_rds("~/GITHUB/naames_bioav_ms/Input/master/processed_bf.2.2020.rds") %>% 
  #omitting CampCN 98 because only measurements at 5 m and below 1500 m were taken. the interpolated values for this cast were as a result, not reliable
  filter(!CampCN == 98) %>% 
  #the multiday N2S4 drifted from the 48˚Lat bin to the 47˚Lat bin, but was a lagrangian study (followed float), so we'll denote all casts within station 4 as being in the 48˚N bin
  mutate(degree_bin = ifelse(Cruise == "AT34" & Station == "4", 48, degree_bin),
         Max_MLD = ifelse(Cruise == "AT34" & Station == "4", 508, Max_MLD)) %>% 
  select(Cruise, Latitude, Longitude, Date, degree_bin, Station, Season, Subregion, CampCN) %>% 
  left_join(., read_xlsx("~/GITHUB/naames_bioav_ms/Input/master/NPP.xlsx", sheet = 2) %>% 
              filter(!station %in% c("NA","N4S2RF")) %>% 
              separate(station, into = c("Cruise", "Station"), sep = "S") %>% 
              select(-Cruise) %>% 
              mutate(Cruise = ifelse(Bloom_phase == "Accumulation phase", "AT39-6", NA),
         Cruise = ifelse(Bloom_phase == "Climax transition", "AT34", Cruise),
         Cruise = ifelse(Bloom_phase == "Equilibrium phase", "AT38", Cruise),
         Cruise = ifelse(Bloom_phase == "Winter transition", "AT32", Cruise)) %>%
           mutate(date = ymd(date),
                  Ez_NPP = (Ez_NPP * 10^3)/12) %>% 
           select(Cruise, date, Station, Ez, Ez_NPP) %>% 
           mutate(Station = ifelse(Station %in% c("4a", "4b", "4c", "4d"), 4, Station),
                  Station = ifelse(Station %in% c("6a", "6b", "6c", "6d", "6e"), 6, Station),
                  Station = as.numeric(Station)) %>% 
           rename(Date = date)) %>% 
  mutate(Cruise = gsub("AT39-6", "AT39", Cruise)) %>% 
  distinct() 


#Bacterial Production and Abundance Data, join with NPP
bact_npp <- read_rds("~/GITHUB/naames_bioav_ms/Input/master/processed_bf.2.2020.rds") %>%
  filter(Target_Z %in% c(0, 5, 10, 25, 50, 75, 100, 150, 200, 300)) %>% 
  select(Cruise:Station, degree_bin, CampCN, Target_Z, BactProd_C, BactProd_C_sd, BactAbund, BactAbund_sd) %>%  
  mutate(Cruise = gsub("AT39-6", "AT39", Cruise)) %>% 
  mutate(degree_bin = ifelse(Cruise == "AT34" & Station == 4, 48, degree_bin)) %>% 
  mutate_at(vars(BactProd_C:BactProd_C_sd), round, 2) %>% 
  filter(!is.na(BactProd_C) | Target_Z == 0) %>% 
  group_by(Cruise, Station, CampCN) %>% 
  fill(BactProd_C:BactAbund_sd, .direction = "up") %>% 
  drop_na(BactProd_C) %>% 
  rename(Depth = Target_Z,
         bp = BactProd_C, 
         sd_bp = BactProd_C_sd,
         ba = BactAbund,
         sd_ba = BactAbund_sd) %>% 
  ungroup() %>% 
  select(Cruise:CampCN, Depth,  everything()) %>% 
  left_join(., npp) %>% 
  select(Cruise:CampCN, Date, Latitude, Longitude,Ez:Ez_NPP, Depth:sd_ba) %>% 
  arrange(CampCN) %>% 
  mutate(ba = ifelse(ba == "NaN", NA, ba),
         sd_ba = ifelse(sd_ba == "NaN", NA, sd_ba))
```

# Interpolate Data for Euphotic Zone Depth

Euphotic zone depths were determined from model output in Fox et al
2020.

``` r
#split the df by CampCN  
add_ez.list <- split(bact_npp, bact_npp$CampCN)

#create a function to add an empty row to each cast, then add the Ez to the Target Z column 
add.func <- function(morty){
  morty[nrow(morty) + 1,] <- NA
  morty$Depth[is.na(morty$Depth)] <- morty$Ez
  rick <- morty %>% 
    fill(., Cruise:Ez_NPP, .direction = c("updown")) %>% 
    arrange(CampCN, Depth)
 }

#apply function to list 
added_ez.list <- lapply(add_ez.list, add.func)

#save the list as a data frame 
added_ez.df <- plyr::ldply(added_ez.list, data.frame) %>% 
  drop_na(Depth) %>% 
  select(-.id)  %>% 
  group_by(Cruise, Station, CampCN, Depth) %>% 
  fill(bp:sd_ba, .direction = "downup") %>% 
  ungroup() %>% 
  distinct()

#split the data frame into lists based on the campaign cast number
to_interpolate.list <- split(added_ez.df, added_ez.df$CampCN)

#create a function that will linearly interpolate each VOI according to the depth intervals of the casts 
interpolate.func <- function(copper) {
to_interpolate.df <- copper %>% 
  select(Depth:ncol(.)) %>% 
  zoo(., order.by = .$Depth) 

interp_bp <- as.numeric(na.approx(to_interpolate.df$bp, na.rm = F))
interp_ba <- as.numeric(na.approx(to_interpolate.df$ba, na.rm = F))
Depth <- to_interpolate.df$Depth
interpolations.df <- data.frame(Depth, interp_bp, interp_ba)
}

#apply function to list 
interpolations.list <- lapply(to_interpolate.list, interpolate.func)

#save the list as a data frame 
interpolations.df <- plyr::ldply(interpolations.list, data.frame) %>% 
  rename(., CampCN = .id) 

#combine the interpolated and non-interpolated data frames
interpolations.df$CampCN <- as.numeric(interpolations.df$CampCN)
interpolated.df <- right_join(bact_npp, interpolations.df) %>% 
  group_by(CampCN) %>% 
  fill(Cruise:degree_bin, Date:Ez_NPP, .direction = "downup") %>% 
  ungroup() %>% 
  mutate_at(vars(interp_bp), round, 2)
```

# Estimate BCD and Specific Growth Rates

We’ll estimate BCD for each cast of each station of all the cruises
using bacterial production, the seasonal averages of BGE, and the
station specific average BGEs (use season ave if not available). We’ll
also estimate specific growth rates (mew) by first converting bacterial
abundance to bacterial carbon using seasonal averages of CCFs and then
dividing bacterial production by bacterial carbon.

BCD and specific growth rates will be integrated through the euphotic
zone.

``` r
bge <- read_rds("~/GITHUB/naames_bioav_ms/Output/bge_summary.rds") %>% 
  select(Season, Cruise, Station, contains("bge")) %>% 
  distinct() %>% 
  left_join(., read_rds("~/GITHUB/naames_bioav_ms/Output/bge_cruise_means.rds") %>% select(Cruise, bge_mean) %>% rename(cruise_bge = bge_mean)) %>% 
  left_join(., read_rds("~/GITHUB/naames_bioav_ms/Output/bge_station_means.rds") %>% select(Cruise, Station, bge_mean) %>% rename(station_bge = bge_mean))
  

#calculate bcd and specific growth rates
bcd <- interpolated.df %>% 
  mutate_at(vars(Station), as.character) %>% 
  left_join(., bge) %>% 
  mutate(Spring_Autumn = ifelse(Season %in% c("Early Spring", "Late Spring"), "Spring", "Autumn")) %>% 
  group_by(Spring_Autumn) %>% 
  fill(c(contains("bge")), .direction = "downup") %>% 
  ungroup() %>% 
  mutate(bcd = round(interp_bp/cruise_bge),
         bcd_global = round(interp_bp/0.26))  #units are  µmol C / m3 / d, bp in nmol C / L / d is equivalent to  µmol C / m3 / d 
```

``` r
int_bcd <- bcd %>% 
  group_by(CampCN) %>% 
  filter(Depth <= Ez) %>% 
  mutate(int.bcd = integrateTrapezoid(Depth, bcd, type = "A"),
         int.bcd_global = integrateTrapezoid(Depth, bcd_global, type = "A"),
         int.bp = integrateTrapezoid(Depth, interp_bp, type = "A"),
         int.ba = integrateTrapezoid(Depth, interp_ba, type = "A")) %>% 
  mutate_at(vars(contains("int.bcd")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains("int.")), funs(./Ez)) %>% 
  mutate(int.NPP = Ez_NPP/Ez,
         bp.npp = int.bp/int.NPP * 100,
         bcd.npp = int.bcd/int.NPP * 100,
         bcd.npp_global = int.bcd_global/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()
```

    ## Warning: funs() is soft deprecated as of dplyr 0.8.0
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## This warning is displayed once per session.

# Merge Experimental DOC with Export MS Data

Here, we put our experimental DOC data in the context of the results
from the export MS to explore remineralization of the seasonally
accumulated DOC pool:

  - We first calculate treatment averages for each experiment
  - For the surface values at each station, we then subtract the DOC
    concentration of the winter/early springtime deeply mixed condition
    (from the export MS) from the DOC concentrations in the experiments
    conducted at the respective station (**norm.doc**).
  - We estimate the concentration of the seasonally accumulated DOC pool
    as the difference between the initial DOC concentration in each
    experiment and the DOC concentration of the mixed condition
    (**accm.doc**)
  - We then estimate the total DOC drawdown for the experiments
    (**bioav.doc**), the total remaining DOC after
    drawdown(**persis.doc**), the percent DOC bioavalability and
    persistence (**per.bioav** and **per.persis**, respectively)
  - Lastly we calculate the rate of DOC drawdown (nmol C L<sup>-1</sup>
    d<sup>-1</sup>,
**ddoc**)

<!-- end list -->

``` r
export <- readRDS("~/GITHUB/naames_bioav_ms/Input/master/processed_export_for_bioavMS.6.7.20.rds") %>% 
  mutate(Cruise = gsub("AT39-6", "AT39", Cruise)) %>% 
  select(Cruise, Season, degree_bin,  Subregion, Station, ave_Ez, sd_Ez, Target_Z,  redis_DOC_vol, int_delta_DOC_ez, int_delta_DOC_100, NCP_mol_ez, NCP_mol_100) %>% 
  distinct() %>% 
  mutate_at(vars(redis_DOC_vol), round, 1) %>% 
  mutate(redis_DOC_vol = ifelse(degree_bin == 44, 55.9, redis_DOC_vol),
         redis_DOC_vol = ifelse(degree_bin == 42, 55.4, redis_DOC_vol),
         redis_DOC_vol = ifelse(degree_bin == 39, 55.5, redis_DOC_vol),
         redis_DOC_vol = ifelse(degree_bin == 51, 52, redis_DOC_vol),
         ) %>%  #following Baetge et al., 2020 use the redistribution of early autumn
  mutate_at(vars(int_delta_DOC_ez, NCP_mol_ez, int_delta_DOC_100, NCP_mol_100), round, 2) %>%
  mutate_at(vars(Station), as.character) %>% 
  drop_na(int_delta_DOC_100) %>% 
  group_by(Cruise, Station) %>% 
  mutate(ave_int_delta_DOC_ez = mean(int_delta_DOC_ez, na.rm = T),
         sd_int_delta_DOC_ez = sd(int_delta_DOC_ez, na.rm = T),
         ave_NCP_mol_ez = mean(NCP_mol_ez, na.rm = T),
         sd_NCP_mol_ez = sd(NCP_mol_ez, na.rm = T)) %>%
  mutate_at(vars(ave_Ez, sd_Ez), round) %>%
  mutate_at(vars(ave_int_delta_DOC_ez, sd_int_delta_DOC_ez, ave_NCP_mol_ez, sd_NCP_mol_ez), round, 2) %>% 
  select(-c(int_delta_DOC_ez, NCP_mol_ez)) %>% 
  distinct() %>% 
  ungroup() %>% 
  select(-Target_Z) %>% 
  distinct()
  

exp_doc <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_bge.rds") %>% 
  left_join(., export %>% select(Cruise, Station, degree_bin) %>% distinct()) %>% 
  mutate(degree_bin = ifelse(is.na(degree_bin), 39, degree_bin)) %>% 
  filter(Treatment == "Control") %>% 
  select(Cruise, Station, degree_bin, Bottle, Hours, Days, stationary.harvest, interp_ptoc, sd_ptoc, interp_doc, sd_doc, i.poc, s.poc, contains("del.doc")) %>% 
  mutate(doc_star = ifelse(Cruise != "AT34" & Hours < stationary.harvest, round(interp_ptoc - i.poc,1), round(interp_ptoc - s.poc,1)),
         doc_star = ifelse(Cruise != "AT34" & Hours > stationary.harvest, interp_ptoc, doc_star),
         combined_doc = ifelse(Cruise == "AT34", interp_doc, doc_star),
          sd_combined_doc = ifelse(Cruise == "AT34", sd_doc, sd_ptoc))  %>%  
  group_by(Cruise, Station, Bottle) %>% 
  mutate(end.remin = last(Hours),
         initial.doc = first(combined_doc),
         longterm.del.doc = ifelse(Cruise != "AT34", first(interp_ptoc) - last(na.omit(interp_ptoc)), first(interp_doc) - last(na.omit(interp_doc))),
         longterm.del.doc = ifelse(longterm.del.doc < 2, NA, longterm.del.doc),
         shortterm.del.doc = ifelse(Cruise == "AT34", del.doc, del.doc.star)) 

export_bioav <- export %>% 
  full_join(., exp_doc %>% 
  select(Cruise, Station, Bottle, stationary.harvest, end.remin, initial.doc, longterm.del.doc,
         shortterm.del.doc) %>% 
  distinct()) %>% 
  mutate(redis_DOC_vol = ifelse(Station %in% c("S2RD", "S2RF"), 55.5, redis_DOC_vol),
         degree_bin = ifelse(is.na(degree_bin), 39, degree_bin),
         Season = ifelse(is.na(Season), "Early Spring", Season)) %>% 
  filter(!Cruise == "AT32") %>% 
  group_by(Cruise, Station) %>% 
  mutate(accm.doc = mean(initial.doc, na.rm = T) - redis_DOC_vol) %>% 
  ungroup() %>% 
  mutate(shortterm.bioav.doc =  ifelse(!is.na(shortterm.del.doc), shortterm.del.doc/accm.doc * 100, NA),
         longterm.bioav.doc =  ifelse(!is.na(longterm.del.doc), longterm.del.doc/accm.doc * 100, NA), 
         persis.doc = accm.doc - longterm.del.doc,
         percent_persis.doc = persis.doc/accm.doc * 100,
         shortterm.ddoc = (shortterm.del.doc/(stationary.harvest/24)) * 1000,
         longterm.ddoc = (longterm.del.doc/(end.remin/24)) * 1000,
         stationary.harvest = stationary.harvest/24,
         end.remin = end.remin/24)

example_doc_remin <- exp_doc %>% 
  filter(degree_bin == 44) %>% 
  left_join(., export_bioav %>% select(Season, Cruise, Station, degree_bin, accm.doc, redis_DOC_vol)) %>% 
  select(Season, Cruise:Days, combined_doc, sd_combined_doc, accm.doc, redis_DOC_vol) %>% 
  distinct() %>% 
  mutate(norm.doc = combined_doc - redis_DOC_vol) %>% 
  group_by(Cruise, Station, Bottle) %>% 
  mutate(del.bioav.doc = norm.doc - last(norm.doc),
         per.pers.doc = norm.doc/accm.doc * 100)
```

# Save Data

``` r
# Integrated BCD
saveRDS(int_bcd, "~/GITHUB/naames_bioav_ms/Output/processed_BCD.rds")
# Integrated DOC and bioavailability
saveRDS(export_bioav, "~/GITHUB/naames_bioav_ms/Output/processed_DOC_bioavailability.rds")
# Lat 44
saveRDS(example_doc_remin, "~/GITHUB/naames_bioav_ms/Output/processed_lat44_remins.rds")
```
