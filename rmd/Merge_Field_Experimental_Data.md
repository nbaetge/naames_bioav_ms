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

# Merge Field NPP, Bacterial Productivity and Abundance Data

``` r
#NPP data
#we'll convert NPP from mg C/ m2 / d to umol C/ m2 / d 
npp <- read_xlsx("~/naames_bioav_ms/Input/NPP for Nick.xlsx", sheet = 2) %>% 
  filter(!station %in% c("NA","N4S2RF")) %>% 
  separate(station, into = c("Cruise", "Station"), sep = "S") %>% 
  select(-Cruise) %>% 
  mutate(Cruise = ifelse(Bloom_phase == "Accumulation phase", "AT39", NA),
         Cruise = ifelse(Bloom_phase == "Climax transition", "AT34", Cruise),
         Cruise = ifelse(Bloom_phase == "Equilibrium phase", "AT38", Cruise),
         Cruise = ifelse(Bloom_phase == "Winter transition", "AT32", Cruise)) %>% 
  select(Cruise, Station, Ez, Ez_NPP) %>% 
  mutate(Station = ifelse(Station %in% c("4a", "4b", "4c", "4d"), 4, Station),
         Station = ifelse(Station %in% c("6a", "6b", "6c", "6d", "6e"), 6, Station),
         Station = as.numeric(Station)) %>% 
  group_by(Cruise, Station) %>% 
  mutate(ave_Ez = mean(Ez),
         sd_Ez = sd(Ez),
         Ez_NPP = Ez_NPP * (10^3/12),
         ave_NPP = mean(Ez_NPP),
         sd_NPP = sd(Ez_NPP)) %>%
  ungroup() %>% 
  select(-c(Ez, Ez_NPP)) %>% 
  distinct() %>% 
  mutate_at(vars(ave_Ez:sd_NPP), round)


#Bacterial Production and Abundance Data, join with NPP
bact_npp <- read_rds("~/naames_bioav_ms/Input/processed_bf.2.2020.rds") %>% 
  filter(Target_Z %in% c(0, 5, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500, 750, 1000, 1250, 1500 )) %>% 
  select(Cruise:Station, degree_bin, CampCN, Target_Z, BactProd_C, BactProd_C_sd, BactAbund, BactAbund_sd) %>%  
  group_by(Cruise, Station, CampCN, Target_Z) %>% 
  ungroup() %>% 
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
  mutate(max_depth = max(Depth, na.rm = T)) %>% 
  ungroup() %>% 
  select(Cruise:CampCN, max_depth, Depth,  everything()) %>% 
  left_join(., npp) %>% 
  select(Cruise:max_depth, ave_Ez:sd_NPP, everything()) 
```

# Interpolate Data for Euphotic Zone Depth

Euphotic zone depths were determined from model output in Fox er al
2020.

``` r
#split the df by CampCN  
add_ez.list <- split(bact_npp, bact_npp$CampCN)

#create a function to add an empty row to each cast, then add the Ez to the Target Z column 
add.func <- function(morty){
  morty[nrow(morty) + 1,] <- NA
  morty$Depth[is.na(morty$Depth)] <- morty$ave_Ez
  rick <- morty %>% 
    fill(., Cruise:sd_NPP, .direction = c("updown")) %>% 
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
  fill(Cruise:max_depth, .direction = "downup") %>% 
  group_by(Cruise, Station, CampCN) %>% 
  fill(ave_Ez:sd_NPP, .direction = "updown") %>% 
  ungroup() %>% 
  mutate_at(vars(interp_bp), round, 2)
```

# Estimate BCD and Specific Growth Rates

We’ll estimate BCD for each station of all the cruises using bacterial
production and the seasonal and global averages of BGE. We’ll also
estimate specific growth rates (mew) by first converting bacterial
abundance to bacterial carbon using seasonal and global averages of CCFs
and then dividing bacterial production by bacterial carbon.

BCD and specific growth rates will be integrated through the euphotic
zone.

``` r
#calculate bcd and specific growth rates
bcd <- interpolated.df %>% 
  left_join(.,   readRDS("~/naames_bioav_ms/Output/processed_bge.rds") %>%
              select(Season, season_bge, global_bge, global_initial_ccf:ave_stationary_ccf) %>%
              distinct()) %>% 
  select(-c(bp:sd_ba)) %>% 
  fill(contains("global"), .direction = "downup") %>% 
  mutate_at(vars(global_initial_ccf:ave_stationary_ccf), round) %>% 
  mutate(season_bge = ifelse(is.na(season_bge), 0.22, season_bge ),
         ave_initial_ccf = ifelse(is.na(ave_initial_ccf), 55, ave_initial_ccf),
         ave_stationary_ccf = ifelse(is.na(ave_stationary_ccf), 31, ave_stationary_ccf),
         
         ###########BCD##################
         #units are  µmol C / m3 / d 
         #bp in nmol C / L / d is equivalent to  µmol C / m3 / d 
         bcd.season_bge = interp_bp/season_bge,
         bcd.global_bge = interp_bp/global_bge,
         
         ##########Bacterial Carbon###########
         # bact carbon in µmol C/m3
         bc.global_i = interp_ba * (10^3) * global_initial_ccf * (1/10^15) * (1/12) * 10^6,
         bc.global_s = interp_ba * (10^3) * global_stationary_ccf * (1/10^15) * (1/12) * 10^6,
         bc.season_i = interp_ba * (10^3) * ave_initial_ccf * (1/10^15) * (1/12) * 10^6,
         bc.season_s = interp_ba * (10^3) * ave_stationary_ccf * (1/10^15) * (1/12) * 10^6)  %>%
  mutate_at(vars(bcd.season_bge, bcd.global_bge), round) %>% 
  ##########Average BA, BP, BC, BCD for each depth of each cruise and station#######
  group_by(Cruise, Station, Depth) %>% 
  mutate(ave_bcd.season_bge = mean(bcd.season_bge, na.rm = T),
         ave_bcd.global_bge = mean(bcd.global_bge, na.rm = T),
         #convert to BA cells / m3
         ave_ba = mean(interp_ba, na.rm = T) * 1000,
         sd_ba = sd(interp_ba, na.rm = T) * 1000,
         ave_bp = mean(interp_bp, na.rm = T),
         sd_bp = sd(interp_bp, na.rm = T),
         ave_bc.global_i = mean(bc.global_i, na.rm = T),
         sd_bc.global_i = sd(bc.global_i, na.rm = T),
         ave_bc.global_s = mean(bc.global_s, na.rm = T),
         sd_bc.global_i = sd(bc.global_i, na.rm = T),
         ave_bc.season_i = mean(bc.season_i, na.rm = T),
         sd_bc.season_i = sd(bc.season_i, na.rm = T),
         ave_bc.season_s = mean(bc.season_s, na.rm = T),
         sd_bc.season_i = sd(bc.season_i, na.rm = T) ) %>% 
  mutate_at(vars(ave_bcd.season_bge, ave_bcd.global_bge), round) %>% 
  ungroup() 
  #########Integrate to Ez##############

int_bcd <- bcd %>% 
  group_by(Cruise, Station) %>% 
  filter(Depth <= ave_Ez) %>% 
  mutate(int.bcd.season_bge = integrateTrapezoid(Depth, ave_bcd.season_bge, type = "A"),
         int.bcd.global_bge = integrateTrapezoid(Depth, ave_bcd.global_bge, type = "A"),
         int.bp = integrateTrapezoid(Depth, ave_bp, type = "A"),
         sd_int.bp = integrateTrapezoid(Depth, sd_bp, type = "A"),
         int.ba = integrateTrapezoid(Depth, ave_ba, type = "A"),
         sd_int.ba = integrateTrapezoid(Depth, sd_ba, type = "A"),
         int.bc.global_i = integrateTrapezoid(Depth, ave_bc.global_i, type = "A"),
         int.bc.global_s = integrateTrapezoid(Depth, ave_bc.global_s, type = "A"),
         int.bc.season_i = integrateTrapezoid(Depth, ave_bc.season_i, type = "A"),
         int.bc.season_s = integrateTrapezoid(Depth, ave_bc.season_s, type = "A")) %>% 
  mutate_at(vars(contains("int.bcd")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains("int.")), funs(./ave_Ez)) %>% 
  mutate(int.NPP = ave_NPP/ave_Ez,
         sd_int.NPP = sd_NPP/ave_Ez, 
  ######## Estimate Specific Growth Rates #####       
         mew.global_i = int.bp/int.bc.global_i,
         mew.global_s = int.bp/int.bc.global_s,
         mew.season_i = int.bp/int.bc.season_i,
         mew.season_s = int.bp/int.bc.season_s,
  ######### Estimate BCD:NPP #########
         bcd.npp.global_bge = int.bcd.global_bge/int.NPP * 100,
         bcd.npp.season_bge = int.bcd.season_bge/int.NPP * 100) %>% 
  select(-c(CampCN, interp_bp, interp_ba, bcd.season_bge:bc.season_s)) %>% 
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
    conducted at the respective station (**norm\_doc**).
  - We estimate the concentration of the seasonally accumulated DOC pool
    as the difference between the initial DOC concentration in each
    experiment and the DOC concentration of the mixed condition
    (**accm\_doc**)
  - We then estimate the total DOC drawdown for the experiments
    (**bioav\_doc**), the total remaining DOC after
    drawdown(**persis\_doc**), the percent DOC bioavalability and
    persistence (**per\_bioav** and **per\_persis**, respectively)
  - Lastly we calculate the rate of DOC drawdown (nmol C L<sup>-1</sup>
    d<sup>-1</sup>,
**ddoc**)

<!-- end list -->

``` r
export <- readRDS("~/naames_bioav_ms/Input/processed_export_for_bioavMS.5.14.20.rds") %>% 
  mutate(Cruise = gsub("AT39-6", "AT39", Cruise)) %>% 
  select(Cruise, Season, degree_bin, Subregion, Station, Target_Z, Ez, ave_N, sd_N, ave_PO4, sd_PO4, ave_Pro:sd_Nano, redis_DOC_vol, NCP_mol_100, doc_ncp_100, int_delta_DOC_100, doc_don_100, int_delta_N_100, int_N_100_vol, int_N_ez, int_PO4_100_vol, int_PO4_ez) %>% 
  distinct() %>% 
  mutate_at(vars(redis_DOC_vol, int_N_100_vol, int_PO4_100_vol), round, 1) %>% 
  mutate_at(vars(NCP_mol_100:doc_ncp_100, int_delta_DOC_100, int_delta_N_100, int_N_ez, int_PO4_ez), round, 2) %>%
  mutate_at(vars(doc_don_100), round) %>% 
  mutate_at(vars(Station), as.character) %>% 
  #error associated with the approach in calculating seasonally accumulated DON led to 17%  of total samples (5 of 29) being negative.  As a result these data were removed from further analyses.
  mutate(doc_don_100 = ifelse(doc_don_100 < 0, NA, doc_don_100), 
         total_phyto = ave_Pro + ave_Syn + ave_Pico + ave_Nano)


export_bioav <- left_join(readRDS("~/naames_bioav_ms/Output/processed_bge.rds") %>% 
                            filter(!Station == "U", Depth == 10, Treatment == "Control") %>% 
                            select(Cruise, Season, Station, Subregion, Bottle, Hours, Days, doc, interp_doc) %>% 
                            distinct() %>% 
                            drop_na(interp_doc), 
                          export %>% 
                            select(Cruise:Station, redis_DOC_vol) %>%  
                            distinct()) %>% 
  mutate(redis_DOC_vol = ifelse(Station %in% c("S2RD", "S2RF"), 55.0, redis_DOC_vol),
         degree_bin = ifelse(is.na(degree_bin), 39, degree_bin)) %>% 
  group_by(Cruise, Station, Hours) %>% 
  mutate(ave_doc = ifelse(!is.na(interp_doc), round(mean(interp_doc, na.rm = T),1), NA),
         sd_doc = ifelse(!is.na(doc), round(sd(doc, na.rm = T),1), NA)) %>% 
  ungroup() %>% 
  group_by(Cruise, Station) %>% 
  mutate(norm_doc = ave_doc - redis_DOC_vol,
         norm_doc_from_t0 = ifelse(!is.na(norm_doc), first(norm_doc) - norm_doc, NA),
         accm_doc = first(ave_doc) - redis_DOC_vol,
         total.bioav_accm_doc = ifelse(last(Days) > 30, last(ave_doc) - redis_DOC_vol, NA),
         st.bioav_accm_doc = ifelse(Days == 7, ave_doc - redis_DOC_vol, NA),
         lt.bioav_accm_doc = total.bioav_accm_doc - st.bioav_accm_doc,
         time.total.bioav_accm_doc = ifelse(last(Days) > 30, last(Days), NA),
         time.lt.bioav_accm_doc = time.total.bioav_accm_doc - 7,
         total.bioav_doc = ifelse(last(Days) > 30, first(norm_doc) - last(norm_doc), NA),
         st.bioav_doc = ifelse(Days == 7, first(norm_doc) - norm_doc, NA),
         lt.bioav_doc = total.bioav_doc - st.bioav_doc,
         total.per_bioav = round((total.bioav_doc/accm_doc * 100)),
         st.per_bioav = round((st.bioav_doc/accm_doc * 100)),
         lt.per_bioav = round((lt.bioav_doc/accm_doc * 100)),
         persis_doc = accm_doc - total.bioav_doc,
         per_persis = round((persis_doc/accm_doc * 100)),
         st.ddoc = round((st.bioav_doc/7) * 1000),
         lt.ddoc = round((lt.bioav_doc/time.lt.bioav_accm_doc) * 1000),
         total.ddoc = round((total.bioav_doc/time.total.bioav_accm_doc) * 1000)) %>%
  fill(accm_doc:lt.ddoc, .direction = "downup") %>% 
  ungroup() 
```

# Save Data

``` r
# Integrated Bacteria
saveRDS(int_bcd, "~/naames_bioav_ms/Output/processed_integrated_BCD.rds")
# Field Data
saveRDS(export, "~/naames_bioav_ms/Output/processed_field.rds")
# Bioavailability
saveRDS(export_bioav, "~/naames_bioav_ms/Output/processed_bioavailability.rds")
```