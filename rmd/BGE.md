BGE
================
Nicholas Baetge
5/20/2020

# Intro

This document shows how processed **individual bottle** data from NAAMES
remineralization bioassays were combined, QC’d, and analyzed. Here, we
calculate derived variables, including BGE.

``` r
library(tidyverse) 
library(readxl)
library(zoo)
library(patchwork)
#stat tests
library(lmtest)
library(lmodel2)
library(rstatix)
library(ggpubr)
```

# Merge and Tidy Processed Data

``` r
ba <- readRDS("~/GITHUB/naames_bioav_ms/Input/Bact_Abund_Input.rds") %>% 
  ungroup() %>% 
  select(-Depth)

oc <- readRDS("~/GITHUB/naames_bioav_ms/Output/filt_processed_doc.rds") %>% 
  select(-c(Depth, Days, key:facet_bottle, doc_from_last:pdoc_from_t0)) 

bc <- readRDS("~/GITHUB/naames_bioav_ms/Output/filt_processed_bacterial_carbon.rds") %>% 
  select(-c(Depth, i.poc.cf1, s.poc.cf1)) %>% 
  rename(i.poc = i.poc.c1.um,
         s.poc = s.poc.c1.um,
         i.ccf = i.ccf.c1,
         s.ccf = s.ccf.c1) # gf75 POC corrected by subtracting universal blank (average of 0.2 and TFF filtrate)


merge <- full_join(ba, oc) %>% 
  left_join(., bc) %>% 
  arrange(Cruise, Station, Treatment, Bottle, Hours) %>%
  mutate(Days = round(Hours/24, 2)) %>% 
  drop_na(Hours) %>% 
  select(Cruise:Treatment, Hours, Days, everything()) %>% 
  group_by(Cruise, Station, Treatment, Bottle) %>% 
  mutate(s.timepoint = ifelse(s.timepoint == Timepoint, Hours, NA)) %>% 
  fill(s.timepoint, .direction = "downup") %>% 
  rename(stationary.harvest = s.timepoint) %>% 
  select(Season, Cruise:Treatment, stationary.harvest, Hours, Days, everything()) %>% 
   select(-c(Timepoint, key)) 
```

# Interpolations

## Split

``` r
to_interpolate <- merge %>% 
  select(Cruise, Station, Treatment, Bottle, Hours, cells, pdoc, doc, toc, ptoc) %>% 
  group_by(Cruise, Station, Treatment, Bottle) # first we'll define the grouping of our dataframe

list <- to_interpolate %>% 
  group_split()  #then we can convert the dataframe into a list, broken up by the groups (list elements)

keys <- to_interpolate %>% 
  group_keys() %>%
  mutate(key = paste(Cruise, Station, Treatment, Bottle, sep = "."))

names(list) <- keys$key
```

## Write the function

``` r
interp.func <- function(x) {
  y <- zoo(x, order.by = x$Hours)
  interp_cells <- round(as.numeric(na.approx(y$cells, na.rm = F)))
  interp_toc <- round(as.numeric(na.approx(y$toc, na.rm = F)), 1)
  interp_ptoc <- round(as.numeric(na.approx(y$ptoc, na.rm = F)), 1)
  interp_doc <- round(as.numeric(na.approx(y$doc, na.rm = F)), 1)
  interp_pdoc <- round(as.numeric(na.approx(y$pdoc, na.rm = F)), 1)
  z <- cbind(y, interp_cells,  interp_toc, interp_ptoc, interp_doc, interp_pdoc)
  as_tibble(z)
}
```

## Apply and Combine

``` r
interpolated <- lapply(list, interp.func) %>% 
  plyr::ldply(., as.data.frame) %>% 
  select(-.id) %>% 
  mutate_at(vars(Hours:interp_pdoc), as.numeric) %>% 
  left_join(merge, .)
```

    ## Joining, by = c("Cruise", "Station", "Bottle", "Treatment", "Hours", "cells", "doc", "toc", "ptoc", "pdoc")

# Calculate Derived Variables

We need to readjust our initial bacterial carbon numbers before
calculating other parameters. Since the initial condition was 1.2-µm
filtrate, the carbon content of that filtrate would be higher than the
30:70 incubation mix at the beginning of the experiment. As a result
we’ll correct the initial carbon number to account for this: Bacterial
POC<sub>initial</sub> = 0.3 \* Bacterial POC~1.2 µm filtrate~

We’ll then calculate and define:

  - BOC (bacterial cell carbon) in µmol C L<sup>-1</sup> by applying
    timepoint specific CCFs
  - DOC\_star = TOC or PTOC - BOC
  - Bacterial growth efficiencies (BGE), only where ∆DOC is resolvable:
      - ∆BC and ∆DOC from T0 to stationary

<!-- end list -->

``` r
calcs <- interpolated %>%
  mutate_at(vars(i.poc), funs(.*0.3)) %>% 
  group_by(Cruise, Station, Treatment, Bottle) %>% 
######cell carbon based on CCFs####
mutate(boc_i.ccf = ifelse(Hours == 0, round((interp_cells * i.ccf) / (12*10^9), 1), NA),
       boc_s.ccf = ifelse(Hours == stationary.harvest, round((interp_cells * s.ccf) / (12*10^9), 1), NA)) %>% 
######delta cell carbon####
  group_by(Cruise, Station, Treatment, Bottle) %>% 
  fill(contains("boc"), .direction = "updown") %>% 
  ungroup() %>% 
  group_by(Cruise, Station, Treatment) %>% 
  fill(contains("boc_i"), .direction = "updown") %>% 
  ungroup() %>% 
  group_by(Cruise, Station, Treatment, Bottle) %>% 
  mutate(del.poc = s.poc - i.poc,
         # del.poc = ifelse(Cruise != "AT34", NA, del.poc),
         del.boc = boc_s.ccf - boc_i.ccf) %>% 
######DOC star######
mutate(doc.star_i = ifelse(Hours == 0, round(ptoc - boc_i.ccf,1), NA),
       doc.star_s = ifelse(Hours == stationary.harvest, round(interp_ptoc - boc_s.ccf,1), NA)) %>% 
  fill(contains("star"), .direction = "updown") %>% 
######delta DOC####  
mutate(del.doc = first(doc) - ifelse(Hours == stationary.harvest, interp_doc, NA),
       del.doc = ifelse(Cruise != "AT34", first(doc) - ifelse(Hours == stationary.harvest, interp_pdoc, NA), del.doc),
       del.doc = ifelse(del.doc < 1.5, NA, del.doc),
       del.doc.star = doc.star_i - doc.star_s,
       del.doc.star = ifelse(del.doc.star < 1.5, NA, del.doc.star),
       del.poc = ifelse(!is.na(del.doc), del.poc, NA)) %>%
  fill(contains("del.doc"), .direction = "downup") %>% 
  ######BGE, points####  
mutate(bge = ifelse(Cruise == "AT34", del.boc/(del.doc), NA),
       bge = ifelse(Cruise != "AT34", del.boc/(del.doc.star), NA),
       bge.poc = del.poc/del.doc) %>% 
  mutate_at(vars(contains("bge")), round, 2) %>% 
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

# Summary Data, Table, and Plots

``` r
bge_summary <- calcs %>% 
  select(Season:Treatment, i.poc, s.poc, i.ccf, boc_i.ccf:bge.poc) %>% 
  distinct() 

bge_summary.table <- bge_summary %>% 
  filter(Treatment == "Control") %>% 
  select(Season, Station, Bottle, del.poc, del.boc, del.doc, del.doc.star, bge, bge.poc, i.ccf) %>% 
  pivot_longer(c(bge, bge.poc), names_to = "calc", values_to = "bge") %>% 
  drop_na(bge) %>% 
  mutate(ave_bge = mean(bge),
         sd_bge = sd(bge)) %>% # global ave bge
  group_by(Season) %>% # cruise average bges
  mutate(ave_ccf_cruise = mean(i.ccf, na.rm = T),
         sd_ccf_cruise = sd(i.ccf, na.rm = T),
         ave_bge_cruise = mean(bge),
         sd_bge_cruise = sd(bge),
         min_bge_cruise = min(bge),
         max_bge_cruise = max(bge)) %>% 
  add_count() %>% 
  ungroup() %>% 
  group_by(Season, Station) %>% # station average bges
  mutate(ave_bge_station = mean(bge),
         sd_bge_station = sd(bge),
         ave_del.poc_station = mean(del.poc, na.rm = T),
         sd_del.poc_station = sd(del.poc, na.rm = T),
         ave_del.boc_station = mean(del.boc, na.rm = T),
         sd_del.boc_station = sd(del.boc, na.rm = T),
         ave_del.doc_station = mean(del.doc, na.rm = T),
         sd_del.doc_station = sd(del.doc, na.rm = T),
         ave_del.doc.star_station = mean(del.doc.star, na.rm = T),
         sd_del.doc.star_station = sd(del.doc.star, na.rm = T),
         ave_ccf_station = mean(i.ccf, na.rm = T),
         sd_ccf_station = sd(i.ccf, na.rm = T)) %>% 
  ungroup() %>% 
  arrange(Season) %>% 
  mutate_at(vars(ave_bge:sd_del.doc.star_station), round, 2) 

bge_summary.table 
```

    ## # A tibble: 24 x 31
    ##    Season Station Bottle del.poc del.boc del.doc del.doc.star i.ccf calc    bge
    ##    <chr>  <chr>   <chr>    <dbl>   <dbl>   <dbl>        <dbl> <dbl> <chr> <dbl>
    ##  1 Early… 2       A        NA      2.5     NA            4.2     77 bge    0.6 
    ##  2 Early… 2       B        NA      2.70    NA            4.20    77 bge    0.64
    ##  3 Early… 3       A        NA      0.80    NA            2.90    60 bge    0.28
    ##  4 Early… 3       B        NA      1.60    NA            2.10    60 bge    0.76
    ##  5 Early… 4       A        NA      0.900   NA            1.9     59 bge    0.47
    ##  6 Early… 5       A        NA      0.4     NA            1.7     37 bge    0.24
    ##  7 Early… 6       A        NA      0.7     NA            2.70    43 bge    0.26
    ##  8 Early… 6       B        NA      0.3     NA            2.40    43 bge    0.13
    ##  9 Early… 1       A         0.98  NA        3.10        NA       NA bge.…  0.32
    ## 10 Early… 3       B         0.74   0.9      2.90        NA        7 bge.…  0.26
    ## # … with 14 more rows, and 21 more variables: ave_bge <dbl>, sd_bge <dbl>,
    ## #   ave_ccf_cruise <dbl>, sd_ccf_cruise <dbl>, ave_bge_cruise <dbl>,
    ## #   sd_bge_cruise <dbl>, min_bge_cruise <dbl>, max_bge_cruise <dbl>, n <dbl>,
    ## #   ave_bge_station <dbl>, sd_bge_station <dbl>, ave_del.poc_station <dbl>,
    ## #   sd_del.poc_station <dbl>, ave_del.boc_station <dbl>,
    ## #   sd_del.boc_station <dbl>, ave_del.doc_station <dbl>,
    ## #   sd_del.doc_station <dbl>, ave_del.doc.star_station <dbl>,
    ## #   sd_del.doc.star_station <dbl>, ave_ccf_station <dbl>, sd_ccf_station <dbl>

# Save Data

``` r
saveRDS(calcs, "~/GITHUB/naames_bioav_ms/Output/processed_bge.rds")

saveRDS(bge_summary.table,  "~/GITHUB/naames_bioav_ms/Output/processed_bge_summary.rds")
```
