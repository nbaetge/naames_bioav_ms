NAAMES DOC Remineralization Bioassays - Bottles
================
Nicholas Baetge
2/29/2020

# Intro

This document shows how **individual bottle** data from NAAMES DOC
remineralization bioassays were processed, QC’d, and analyzed.

``` r
library(tidyverse) 
library(rmarkdown)
library(knitr)
library(data.table) 
library(scales)
library(zoo)
library(oce)
#rmarkdown tables
library(stargazer)
library(pander)
library(growthcurver)
#stat tests
library(lmtest)
library(lmodel2)
```

# Bacterial Abundance

Here, growth curves are fitted to the logistic equation using the
GrowthCurver package. The timing of the onset of the stationary growth
phase is determined from the model
fit.

## Import and Tidy Data

``` r
ba.df <- read_csv("~/naames_bioav_ms/Input/N2-4_BactA_Remin_Master.csv") %>% 
  mutate(Season = Cruise,
         Season = gsub("AT34", "Late Spring", Season),
         Season = gsub("AT38", "Early Autumn", Season),
         Season = gsub("AT39", "Early Spring", Season),
         Treatment = ifelse(Bottle == "Niskin", "Niskin", Treatment),
         Treatment = ifelse(Bottle == "GF75", "GF75", Treatment),
         Treatment = gsub("TW15", "TW20", Treatment),
         Treatment = ifelse(is.na(Treatment), Bottle, Treatment),
         Treatment = ifelse(Treatment %in% c("100mls", "250mls", "500mls", "1L"), "Volume", Treatment),
         Treatment = gsub("Surface", "Parallel", Treatment),
         Treatment_Btl = ifelse(Bottle == "Niskin", "Niskin", Treatment_Btl),
         Treatment_Btl = ifelse(Bottle == "GF75" & is.na(Treatment_Btl), "GF75", Treatment_Btl),
         Treatment_Btl = gsub("TW15", "TW20", Treatment_Btl),
         Treatment_Btl = gsub("TW12-A", "TW12-C", Treatment_Btl),
         Treatment_Btl = gsub("TW12-B", "TW12-D", Treatment_Btl),
         Treatment_Btl = gsub("TW13-A", "TW12-E", Treatment_Btl),
         Treatment_Btl = gsub("TW13-B", "TW12-F", Treatment_Btl),
         Treatment_Btl = ifelse(is.na(Treatment_Btl), Bottle, Treatment_Btl),
         Bottle = gsub("SA", "A", Bottle),
         Bottle = gsub("SB", "B", Bottle),
         Bottle = gsub("DA", "C", Bottle),
         Bottle = gsub("DB", "D", Bottle),
         Bottle = ifelse(Bottle == "GF75", Treatment_Btl, Bottle),
         Bottle = gsub("Control-", "", Bottle),
         Bottle = gsub("TW12-", "", Bottle),
         Bottle = gsub("TW13-", "", Bottle),
         Bottle = gsub("MixDS-", "", Bottle),
         Bottle = gsub("MixSD-", "", Bottle),
         Bottle = gsub("SynExd-", "", Bottle),
         Bottle = gsub("TWExd-", "", Bottle),
         Bottle = gsub("SynLys-", "", Bottle),
         Bottle = gsub("TW5-", "", Bottle),
         Bottle = gsub("TW10-", "", Bottle),
         Bottle = gsub("TW20-", "", Bottle),
         Bottle = gsub("12", "", Bottle),
         Bottle = gsub("13", "", Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW12" & Bottle == "A", gsub("A", "C", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW12" & Bottle == "B", gsub("B", "D", Bottle), Bottle),
          Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW13" & Bottle == "A", gsub("A", "E", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW13" & Bottle == "B", gsub("B", "F", Bottle), Bottle),
          Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "Control" & Bottle == "C", gsub("C", "G", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "Control" & Bottle == "D", gsub("D", "H", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW12" & Bottle == "C", gsub("C", "I", Bottle),Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW12" & Bottle == "D", gsub("D", "J", Bottle), Bottle),
          Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW13" & Bottle == "C", gsub("C", "K", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW13" & Bottle == "D", gsub("D", "L", Bottle), Bottle),
          Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "GF75" & Bottle == "A", gsub("A", "G", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "GF75" & Bottle == "B", gsub("B", "H", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "GF75" & Bottle == "C", gsub("C", "I", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "GF75" & Bottle == "D", gsub("D", "J", Bottle), Bottle),
          Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "GF75" & Bottle == "E", gsub("E", "K", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "GF75" & Bottle == "F", gsub("F", "L", Bottle), Bottle),
          Bottle = ifelse(Cruise == "AT34" & Station == 1 & Depth == 200 & Treatment == "GF75" & Bottle == "A", gsub("A", "C", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 1 & Depth == 200 & Treatment == "GF75" & Bottle == "B", gsub("B", "D", Bottle), Bottle),
           Bottle = ifelse(Cruise == "AT34" & Station == 4 & Depth == 200 & Treatment == "GF75" & Bottle == "A", gsub("A", "C", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 4 & Depth == 200 & Treatment == "GF75" & Bottle == "B", gsub("B", "D", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 5 & Depth == 200 & Treatment == "GF75" & Bottle == "A", gsub("A", "C", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 5 & Depth == 200 & Treatment == "GF75" & Bottle == "B", gsub("B", "D", Bottle), Bottle),
         Depth = gsub(5, 10, Depth),
         Depth = as.numeric(Depth)) %>% 
  rename(cells = cellsperL,
         sd_cells = sd_cellsperL,
         p_cells = p_cellsperL,
         sd_p_cells = p_sdcellsperL) %>% 
  drop_na(cells) %>% 
  select(-Treatment_Btl) %>% 
  select(Season, everything()) %>% 
  arrange(Cruise, Station, Depth, Treatment, Bottle, Hours)

ba.df$Season <- factor(ba.df$Season, levels = levels)
```

## Inspect Growth Curves

To be able to calculate derived variables, such as carbon per cell and
BGE, we’ll need to inspect the growth curves and define the stationary
phase for each experiment.

We are using the growthcurver package to fit growth curve data to the
standard form of the logistic equation and then return a data table with
population-level information. Specifically, the model outputs we are
really only interested in are:

  - t\_mid, the time at which 1/2 carrying capacity is reached. We will
    double this estimate to attain the time at which carrying capacity
    is reached (stationary)
  - sigma, a measure of the goodnesss of fit of the parameters of the
    logistic equation for the data; it is the residual sum of squares
    from the nonlinear regression model. Smaller sigma values indicate a
    better fit of the logistic curve to the data than larger values.
  - df, degrees of freedom

<!-- end list -->

``` r
gc_input <- ba.df %>% 
  filter(!Treatment %in% c("Niskin","GF75", "Volume", "Parallel", "TFF-Ret")) %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) %>% 
  #calculate the running delta of nat log cells for each exp
  mutate(lncells = round(log(cells), 2),
         delta_lncells = lncells - lncells[which.min(Hours)],
         delta_lncells = ifelse(delta_lncells < 0, 0, delta_lncells)) 

gc_input_keys <- gc_input %>% 
  group_keys() %>% 
  mutate(key = paste(Cruise, ", S", Station, ", Z =", Depth, ",", Treatment, ",", Bottle))
gc_input_header <- gc_input_keys$key

gc_input_list <- gc_input %>% 
  group_split()
names(gc_input_list) <- gc_input_header
```

``` r
gcplot.func <- function(x){
  gc_fit <- SummarizeGrowth(x$Hours, x$delta_lncells) 
}
gcplot.list <- lapply(gc_input_list, gcplot.func)
```

For the sake of keeping the length of this document as short as
possible, the growth curve model plots are not included here. However,
they are included below in the section “Re-evaluate Growth Curves”.

## Re-evaluate Growth Curves

To improve the model fits and thus, the timing estimates of the
stationary transition for each experiment, we need to evaluate each
growth curve, manually filtering out death/secondary growth phases and
poorly replicated samples that skew the model results.

``` r
eval_gc_input <- gc_input %>% 
  mutate(key = paste(Cruise, Station, Depth, Treatment, Bottle, Timepoint, sep = ".")) %>% 
  select(Season:Treatment, key, everything()) %>% 
  filter(!Hours >= 300) %>% 
  filter(!key %in% c("AT34.1.10.Control.A.1","AT34.1.10.Control.B.1", "AT34.1.200.Control.C.3", "AT34.1.200.Control.D.3","AT34.2.10.Control.A.4", "AT34.2.10.Control.A.9", "AT34.2.10.Control.A.10", "AT34.2.10.Control.A.11","AT34.2.10.Control.A.4", "AT34.2.10.Control.B.7", "AT34.2.10.Control.B.9", "AT34.2.10.Control.A.11", "AT34.2.10.TW12.C.5",
  "AT34.2.10.TW13.E.2", "AT34.2.10.TW13.F.3", "AT34.2.200.Control.G.2", "AT34.2.200.Control.G.5", "AT34.2.200.Control.G.6", "AT34.2.200.Control.G.7", "AT34.2.200.Control.H.5", "AT34.2.200.Control.H.6", "AT34.2.200.Control.H.7", "AT34.2.200.TW12.I.2", "AT34.2.200.TW12.J.2", "AT34.2.200.TW13.K.2", "AT34.3.10.Control.A.8", "AT34.3.10.Control.A.9", "AT34.3.10.Control.A.10", "AT34.3.10.Control.B.8", "AT34.3.10.Control.B.9", "AT34.3.10.Control.B.10", "AT34.3.200.Control.C.4", "AT34.3.200.Control.C.6", "AT34.3.200.Control.C.9", "AT34.3.200.Control.D.7", "AT34.4.200.Control.C.5", "AT34.4.200.Control.C.6", "AT34.4.200.Control.C.7", "AT34.4.200.Control.D.5", "AT34.4.200.Control.D.6", "AT34.4.200.Control.D.7", "AT34.5.10.Control.A.5", "AT34.5.10.Control.B.5", "AT34.U.10.Control.A.5", "AT34.U.10.Control.B.8", "AT34.U.10.Control.B.1", "AT38.1.10.Control.A.1", "AT38.1.10.Control.B.7", "AT38.1.10.MixSD.G.6", "AT38.1.10.SynExd.I.8", "AT38.1.10.SynExd.I.9", "AT38.1.10.SynExd.I.10", "AT38.1.10.SynExd.J.8", "AT38.1.10.SynExd.J.9", "AT38.1.10.SynExd.J.10", "AT38.1.10.SynLys.M.9", "AT38.1.10.SynLys.M.10", "AT38.1.10.SynLys.N.9", "AT38.1.10.SynLys.N.10", "AT38.1.200.Control.F.5", "AT38.1.200.Control.F.6", "AT38.1.200.Control.E.5", "AT38.1.200.Control.E.6", "AT38.1.200.MixDS.C.5", "AT38.1.200.MixDS.C.8", "AT38.1.200.MixDS.D.5", "AT38.1.200.MixDS.D.6", "AT38.4.10.Control.A.7", "AT38.4.10.Control.B.7", "AT38.4.10.Control.B.8", "AT38.5.10.Control.A.9", "AT38.5.10.Control.A.10", "AT38.5.10.Control.A.11", "AT38.5.10.Control.B.7"))

eval_gc_input_keys <- eval_gc_input %>% 
  group_keys() %>% 
  mutate(key = paste(Cruise, ", S", Station, ", Z =", Depth, ",", Treatment, ",", Bottle)) 
eval_gc_input_header <- eval_gc_input_keys$key

eval_gc_input_list <- eval_gc_input %>% 
  group_split()
names(eval_gc_input_list) <- eval_gc_input_header

eval_gcplot.list <- lapply(eval_gc_input_list, gcplot.func)
```

We’ve omitted 88 of the 768 (~11%) of the observations to improve the
model fit. 47% of those omitted observations were taken after 300 hours,
which for many experiments constituted the death phase. For the sake of
space, we won’t include the growthcurver plots in this document, but
they can be found in the .rmd

## Wrangle GrowthCurver Data

We’ll add the model output data to the original dataframe, but first
we’ll tidy the model output and look at the distribution of the model
fits.

We’ll determine thresholds for sigma (smaller sigma values indicate a
better fit) and df for omission of model fit data.

``` r
gcdata.func <- function(y){
  gc.fit <- SummarizeGrowth(y$Hours, y$delta_lncells)
  gcdata.df <- as.data.frame(unlist(gc.fit$vals), stringsAsFactors = FALSE)
}

gcdata <- eval_gc_input_list %>% 
  lapply(., gcdata.func) %>% 
  #save the list as a data frame %>% 
  data.frame(.) %>% 
  t(.) %>% 
  as_data_frame(.) %>% 
  mutate_at(., vars(k:auc_e), as.numeric) %>% 
  mutate_at(vars(sigma, t_mid), round, 2) %>% 
  #calculate time to stationary for each curve
  mutate(stationary = round(t_mid*2)) %>% 
  mutate(key = eval_gc_input_header) %>% 
  ungroup() 
```

    ## Warning: `as_data_frame()` is deprecated, use `as_tibble()` (but mind the new semantics).
    ## This warning is displayed once per session.

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

Based on these two histograms, we’ll omit model fit data where sigma \>
0.21 and df \< 4. In these experiments, we cannot reliably distinguish
the onset of stationary.

``` r
gcmerge.df <- eval_gc_input %>% 
  mutate(key = paste(Cruise, ", S", Station, ", Z =", Depth, ",", Treatment, ",", Bottle)) %>% 
  full_join(., gcdata %>% filter(!sigma > 0.21 | !df < 4)) %>% 
  ungroup() 
```

## Add Stationary Timepoint

We’ll add a row to each experiment that includes the stationary
timepoint.

``` r
gcmerge.list <- gcmerge.df %>% 
  #Split the dataframe into lists
  split(.,  .$key)

#add row to each list element
newrow.func <- function(morty){
  morty[nrow(morty) + 1,] <- NA
  rick <- morty %>% 
    fill(Season:Bottle, Treatment, key, stationary)
  rick
}

#apply function to list
gcstat.list <- lapply(gcmerge.list, newrow.func)

gcstat.df <- plyr::ldply(gcstat.list, data.frame) %>% 
  select(Season:Hours, stationary, cells:auc_e) %>%
  mutate(Hours = ifelse(is.na(Hours), stationary, Hours))
```

# Organic Carbon

## Import and Tidy Data

We will need to tidy the data in the same way the bacterial abundance
dataframe was tidied before merging the datasets. We’ll then omit data
low quality data i.e. poor replication or possible contamination. We
observed contamination in many of the NAAMES 2 experiments after
incubations were transferred to
WHOI.

``` r
oc.df <- read_csv("~/naames_bioav_ms/Input/N2-4_DOC_Remin_Master.csv") %>%
  mutate(Season = Cruise,
         Season = gsub("AT34", "Late Spring", Season),
         Season = gsub("AT38", "Early Autumn", Season),
         Season = gsub("AT39", "Early Spring", Season),
         Treatment = ifelse(Bottle == "Niskin", "Niskin", Treatment),
         Treatment = ifelse(Bottle == "GF75", "GF75", Treatment),
         Treatment = gsub("TW15", "TW20", Treatment),
         Treatment = ifelse(is.na(Treatment), Bottle, Treatment),
         Treatment = ifelse(Treatment %in% c("Surface", "Deep"), "Parallel", Treatment),
         Treatment = ifelse(Treatment %in% c("100mls", "250mls", "500mls", "1L"), "Volume", Treatment),
         Treatment_Btl = ifelse(Bottle == "Niskin", "Niskin", Treatment_Btl),
         Treatment_Btl = ifelse(Bottle == "GF75" & is.na(Treatment_Btl), "GF75", Treatment_Btl),
         Treatment_Btl = gsub("TW15", "TW20", Treatment_Btl),
         Treatment_Btl = gsub("TW12-A", "TW12-C", Treatment_Btl),
         Treatment_Btl = gsub("TW12-B", "TW12-D", Treatment_Btl),
         Treatment_Btl = gsub("TW13-A", "TW12-E", Treatment_Btl),
         Treatment_Btl = gsub("TW13-B", "TW12-F", Treatment_Btl),
         Treatment_Btl = ifelse(is.na(Treatment_Btl), Bottle, Treatment_Btl),
         Bottle = gsub("SA", "A", Bottle),
         Bottle = gsub("SB", "B", Bottle),
         Bottle = gsub("DA", "C", Bottle),
         Bottle = gsub("DB", "D", Bottle),
         Bottle = ifelse(Bottle == "GF75", Treatment_Btl, Bottle),
         Bottle = gsub("Control-", "", Bottle),
         Bottle = gsub("TW12-", "", Bottle),
         Bottle = gsub("TW13-", "", Bottle),
         Bottle = gsub("MixDS-", "", Bottle),
         Bottle = gsub("MixSD-", "", Bottle),
         Bottle = gsub("SynExd-", "", Bottle),
         Bottle = gsub("TWExd-", "", Bottle),
         Bottle = gsub("SynLys-", "", Bottle),
         Bottle = gsub("TW5-", "", Bottle),
         Bottle = gsub("TW10-", "", Bottle),
         Bottle = gsub("TW20-", "", Bottle),
         Bottle = gsub("12", "", Bottle),
         Bottle = gsub("13", "", Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW12" & Bottle == "A", gsub("A", "C", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW12" & Bottle == "B", gsub("B", "D", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW13" & Bottle == "A", gsub("A", "E", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW13" & Bottle == "B", gsub("B", "F", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "Control" & Bottle == "C", gsub("C", "G", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "Control" & Bottle == "D", gsub("D", "H", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW12" & Bottle == "C", gsub("C", "I", Bottle),Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW12" & Bottle == "D", gsub("D", "J", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW13" & Bottle == "C", gsub("C", "K", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW13" & Bottle == "D", gsub("D", "L", Bottle), Bottle),
         Depth = gsub(5, 10, Depth),
         Depth = as.numeric(Depth),
         key = paste(Cruise, Station, Depth, Treatment, Bottle, Timepoint, sep = "."),
         PTOC_ave = ifelse(Treatment == "Parallel", TOC_ave, PTOC_ave),
         PTOC_sd = ifelse(Treatment == "Parallel", TOC_sd, PTOC_sd)) %>% 
  rename(doc = DOC_ave,
         sd_doc = DOC_sd,
         toc = TOC_ave, 
         sd_toc = TOC_sd, 
         ptoc = PTOC_ave, 
         sd_ptoc = PTOC_sd, 
         pdoc = PDOC_ave, 
         sd_pdoc = PDOC_sd) %>%
  filter(!Treatment %in% c("Niskin", "GF75", "TFF-Ret", "Volume")) %>% 
  select(-Treatment_Btl) %>% 
  select(Season, everything()) %>% 
  arrange(Cruise, Station, Depth, Treatment, Bottle, Hours)  %>% 
  filter(!is.na(Hours))
```

Because a different suite of organic carbon samples were taken on each
cruise, we’ll separate the cruise datasets, drop NA values, quality
filter them, and then recombine the datasets. To quality filter, we’ll
calculate rolling ∆DOC values i.e. the diffence between the one
observation and the one before it. We’ll omit the data where ∆DOC
indicates contamination due to transport etc. after the 21 d shipboard
occupation (increase greater than our detection limit, 1.5 µmol C
L<sup>-1</sup>).

``` r
n2_oc <- oc.df %>% 
  filter(Cruise == "AT34") %>% 
  drop_na(doc) %>%
  group_by(Station, Depth, Treatment, Bottle) %>%
  mutate(doc_from_last = lag(doc) - doc) %>% 
  filter(!doc_from_last <= -1.5 | is.na(doc_from_last) | !Hours > 504)  %>% 
  mutate(doc_from_last = lag(doc) - doc) %>% 
  filter(!doc_from_last <= -1.5 | is.na(doc_from_last) | !Hours > 504) %>% 
  mutate(doc_from_last = lag(doc) - doc) %>% 
  filter(!doc_from_last <= -1.5 | is.na(doc_from_last) | !Hours > 504) %>% 
  mutate(doc_from_last = lag(doc) - doc,
         doc_from_t0 = first(doc) - doc) %>% 
  filter(!Station == 4 | !Hours == 192) %>% 
  ungroup()

n3_4_oc <- oc.df %>% 
  filter(!Cruise == "AT34") %>% 
  drop_na(ptoc) %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) %>% 
  mutate(ptoc_from_last = lag(ptoc) - ptoc) %>% 
  filter(!ptoc_from_last <= -1.5 | is.na(ptoc_from_last) | !Hours > 504) %>% 
  mutate(ptoc_from_last = lag(ptoc) - ptoc,
         ptoc_from_t0 = first(ptoc) - ptoc,
         doc_from_t0 = first(doc) - doc,
         toc_from_t0 = first(toc) - toc, 
         pdoc_from_t0 = first(pdoc) - pdoc) %>%
  ungroup()

oc_p <- bind_rows(n2_oc, n3_4_oc) %>% 
  mutate(Days = round(Hours/24, 1),
         facet_depth = paste(Depth, "m"), 
         facet_treatment = Treatment, 
         facet_treatment = gsub("MixSD", "Surface Comm., Deep DOM", facet_treatment),
         facet_treatment = gsub("MixDS", "Deep Comm., Surface DOM", facet_treatment),
          facet_treatment = gsub("SynExd", "Synechococcus Exudate", facet_treatment),
         facet_treatment = gsub("SynLys", "Synechococcus Lysate", facet_treatment),
          facet_treatment = gsub("TWExd", "T. Weissflogii Exudate", facet_treatment),
         facet_treatment = gsub("TW5", "+5 µmol C/L Exudate", facet_treatment),
          facet_treatment = gsub("TW10", "+10 µmol C/L Exudate", facet_treatment),
          facet_treatment = gsub("TW20", "+20 µmol C/L Exudate", facet_treatment),
           facet_treatment = gsub("TW12", "T. Weissflogii 12C Exudate", facet_treatment),
          facet_treatment = gsub("TW13", "T. Weissflogii 13C Exudate", facet_treatment),
         facet_bottle = Bottle,
         facet_bottle = gsub("W", "Whole Seawater", facet_bottle),
          facet_bottle = gsub("1.2", "1.2 µm Filtrate", facet_bottle), 
          facet_bottle = gsub("NV", "1.2 µm Filtrate: TFF Filtrate (3:7)", facet_bottle)) 

oc_p$Treatment <- factor(oc_p$Treatment, levels = levels)
oc_p$facet_treatment <- factor(oc_p$facet_treatment, levels = levels)
oc_p$facet_bottle <- factor(oc_p$facet_bottle, levels = levels)
```

## Plot Curves

### NAAMES 2

#### No Addition: Long-term

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

#### No Addition: Short-term

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

#### Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

There looks to be some poor replication between the experiments
conducted at station 4 and U:

  - S4 bottle A measurements look suspect. additionally, the increase
    between day 6 and 7 for both bottles are suspect
  - SU bottle B measurements at 24 h is suspect relative to the downward
    trend in both bottles

These points will affect BGE calculations and my need to be omitted. The
BGE from S4 Bottle A may be unreliable.

BGEs cannot be calculated for the 200 m control experiments given that
not DOC drawdown was
observed.

### NAAMES 3

#### No Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

#### Mixing

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

#### Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

### NAAMES 4

#### No Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-31-1.png" style="display: block; margin: auto;" />

The replication between experiments at Station 4 is
poor.

#### Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-33-1.png" style="display: block; margin: auto;" />

#### Variable Set-Up

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-34-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-35-1.png" style="display: block; margin: auto;" />

# Bacterial Carbon

POC from DOC remineralization bioassays was estimated by filtering cells
onto double-stacked GF75 filters (0.3 µm pore size) and then running
each filter on an elemental analyzer at Bigelow (Maine). This document
shows how bacterial carbon estimates were calculated from those
measurements.

Upon importing the dataset, we’ll add two categorical variables that
indicate what type of filtrate was used as the “blank” water. On NAAMES
2, 0.2 µm-filtered water was used and on NAAMES 4, TFF-filtered water
was used. Unfortunately, blanks were not taken on NAAMES
3.

``` r
bigelow <- read_csv("~/naames_bioav_ms/Input/N2-4_BactC_Processing.csv") %>% 
  select(-c(Treatment_Btl:Days)) %>% 
  rename(Filter = Rep) %>% 
  mutate(BlnkFiltrate = ifelse(Cruise == "AT34", "0.2 µm", "TFF"),
         BlnkFiltrate = ifelse(Cruise == "AT38", NA, BlnkFiltrate),
         facetlabel = paste(Depth, " m"),
         Season = Cruise,
         Season = gsub("AT34", "Late Spring", Season),
         Season = gsub("AT38", "Early Autumn", Season),
         Season = gsub("AT39", "Early Spring", Season),
         Treatment = ifelse(is.na(Treatment), "Niskin", Treatment),
         Treatment = ifelse(Bottle %in% c("SA13", "SB13", "DA13", "DB13"), "TW13", Treatment),
         Bottle = gsub("SA", "A", Bottle),
         Bottle = gsub("SB", "B", Bottle),
         Bottle = gsub("DA", "C", Bottle),
         Bottle = gsub("DB", "D", Bottle), 
         Timepoint = ifelse(Treatment == "Niskin", 0, Timepoint)
         ) %>% 
  mutate_at(vars(BactN_µg, BactC_µg), round, 1)

bigelow$Season <- factor(bigelow$Season, levels = levels)
```

Here, we plot all of the bacterial POC data from all of the experiments
across all of the cruises, including blanks. We’ve cut off the max value
to 25 µg C, but note there are samples with much greater
values.

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-37-1.png" style="display: block; margin: auto;" />
There is a large spread in the data points for the second filter of the
non-blank samples on the early autumn cruise, with higher values
suggesting more material (colloidal? bacteria? detritus?) slipped
through the first filter onto the second filter. It is unfortunate we
don’t have blanks for this cruise to help diagnose this. It may be
inappropriate to apply a universal blank correction determined from the
spring cruises to this cruise.

### Determine Blank Values

``` r
blank.data <- bigelow %>% 
  filter(Bottle == "Blank") %>% 
  #filter out outliers
  # filter(!Season == "Late Spring" | !Station == 1 | !Depth == 10) %>% 
  # filter(!Season == "Early Spring" |  !Depth == 200) %>% 
  mutate(BactC_µM = round(BactC_µg/12,2),
         BactN_µM = round(BactN_µg/14,2),
         BactC_N = ifelse(BactN_µM > 0, round(BactC_µM/BactN_µM,1), NA))
```

The 0.2 µm filter 2 value of 6.46 and the TFF filter 1 value of 9.8
appears to be outliers relative to all the other filters (across cruises
and filtrate types), indicating possible contamination or the
compromised integrity of the first filter. These samples (filter 1 & 2)
are going to be perserved in the analysis for now.

#### Filter Summary Statistics

##### t-tests

``` r
t.test.data <- bigelow %>% 
  filter(!Bottle %in% c("100mls", "250mls", "500mls")) %>%
  # filter(!Season == "Late Spring" | !Station == 1 | !Depth == 10 | !Bottle == "Blank") %>% 
  # filter(!Season == "Early Spring" |  !Depth == 200 | !Bottle == "Blank") %>% 
  mutate(Treatment = ifelse(Bottle == "Blank", "Blank", Treatment),
         Treatment = ifelse(Treatment == "Blank" & BlnkFiltrate == "TFF", "TFF", Treatment),
         Treatment = ifelse(Treatment == "Blank" & BlnkFiltrate == "0.2 µm", "0.2 µm", Treatment),
         Treatment = gsub("Niskin", "1.2 µm", Treatment),
         Type = ifelse(Bottle == "Blank", "Blank", "Non-blank")) %>% 
  select(Season, Treatment, Type, facetlabel, Filter, BactC_µg)  %>% 
  rename(Filtrate = Treatment,
         Depth = facetlabel)
```

###### 0.2 µm filtrate, filter 1 v filter 2

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Filter
    ## t = 1.1028, df = 17.827, p-value = 0.2848
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5619966  1.8019966
    ## sample estimates:
    ## mean in group 1 mean in group 2 
    ##            3.56            2.94

**Not significantly different**

###### TFF filtrate, filter 1 v filter 2

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Filter
    ## t = 2.7311, df = 6.2211, p-value = 0.03291
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.266534 4.504895
    ## sample estimates:
    ## mean in group 1 mean in group 2 
    ##        5.385714        3.000000

**Significantly different**

###### 0.2 µm and TFF filtrate, filter 1

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Filtrate
    ## t = -1.9334, df = 8.299, p-value = 0.08794
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -3.9897395  0.3383109
    ## sample estimates:
    ## mean in group 0.2 µm    mean in group TFF 
    ##             3.560000             5.385714

**Not signficantly different**

###### 0.2 µm and TFF filtrate, filter 2

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Filtrate
    ## t = -0.1386, df = 10.39, p-value = 0.8924
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.019692  0.899692
    ## sample estimates:
    ## mean in group 0.2 µm    mean in group TFF 
    ##                 2.94                 3.00

**Not significantly different**

###### All blanks, filter 1 v filter 2

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Filter
    ## t = 2.5785, df = 24.297, p-value = 0.0164
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2695235 2.4245942
    ## sample estimates:
    ## mean in group 1 mean in group 2 
    ##        4.311765        2.964706

**Significantly different**

###### Blank and Non-blank (initial condition), filter 1

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Type
    ## t = -5.4128, df = 27.22, p-value = 9.812e-06
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -10.843868  -4.884132
    ## sample estimates:
    ##     mean in group Blank mean in group Non-blank 
    ##                   3.560                  11.424

**Significantly different**

###### Blank and Non-blank (initial condition), filter 2

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Type
    ## t = -2.2595, df = 25.323, p-value = 0.0327
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -6.4116635 -0.2989247
    ## sample estimates:
    ##     mean in group Blank mean in group Non-blank 
    ##                2.964706                6.320000

**Significantly different**

###### Blank and Non-blank, filter 2, without N3

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  BactC_µg by Type
    ## t = -1.8834, df = 31.003, p-value = 0.06906
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.60097600  0.06372109
    ## sample estimates:
    ##     mean in group Blank mean in group Non-blank 
    ##                2.964706                3.733333

**Not signficantly different**

Based on these results, it seems like filter 2 of the blanks and the
initial conditions could be a good blank.

For the filters that had 0.2 µm filtrate run through them:

\-Any particles remaining in the 0.2 µm filtrate are equally caught on
the top and bottom GF75 filters -Handling of the filters is the same
-DOC adsorption is the same on both filters

We would expect to see the same for the TFF filters, but the top and
bottom filters are significantly different. As a result these filters
drive a significant difference between filter 1 and filter 2 when all
blanks are included.

However, given that filter 2 of the blanks and non-NAAMES 3\* 1.2 µm
filtrate are not signficantly different, it seems like filter 2 could
represent a good blank correction.

\*The lower filter cell retention and the higher filter 2 POC values for
NAAMES3 would suggest slippage, making omission of this data reasonable
if blank correction factors are based on blank and non-blank samples.

##### Summary Stats

``` r
blank.stats <- blank.data %>% 
  group_by(Season, Depth, Filter) %>% 
  select(Season, Depth, Filter, BactC_µg, BactN_µg, BactC_µM:BactC_N) %>% 
  summarize(ave_c_µg = mean(BactC_µg, na.rm = T),
            sd_c_µg = sd(BactC_µg, na.rm = T),
            max_c_µg = max(BactC_µg, na.rm = T),
            min_c_µg = min(BactC_µg, na.rm = T),
            ave_n_µg = mean(BactN_µg, na.rm = T),
            sd_n_µg = sd(BactN_µg, na.rm = T),
            max_n_µg = max(BactN_µg, na.rm = T),
            min_n_µg = min(BactN_µg, na.rm = T),
            
            ave_c_µM = mean(BactC_µM, na.rm = T),
            sd_c_µM = sd(BactC_µM, na.rm = T),
            max_c_µM = max(BactC_µM, na.rm = T),
            min_c_µM = min(BactC_µM, na.rm = T),
            ave_n_µM = mean(BactN_µM, na.rm = T),
            sd_n_µM = sd(BactN_µM, na.rm = T),
            max_n_µM = max(BactN_µM, na.rm = T),
            min_n_µM = min(BactN_µM, na.rm = T),
            
            ave_c_n_µM = mean(BactC_N, na.rm = T),
            sd_c_n_µM = sd(BactC_N, na.rm = T),
            max_c_n_µM = max(BactC_N, na.rm = T),
            min_c_n_µM = min(BactC_N, na.rm = T),
            
            n()) %>% 
  mutate_at(vars(ave_c_µg:min_c_n_µM), round, 2) %>% 
  arrange(Filter, Depth,Season)
```

| Season       | Depth | Filter | ave\_c\_µg | sd\_c\_µg | max\_c\_µg | min\_c\_µg | ave\_n\_µg | sd\_n\_µg | max\_n\_µg | min\_n\_µg | ave\_c\_µM | sd\_c\_µM | max\_c\_µM | min\_c\_µM | ave\_n\_µM | sd\_n\_µM | max\_n\_µM | min\_n\_µM | ave\_c\_n\_µM | sd\_c\_n\_µM | max\_c\_n\_µM | min\_c\_n\_µM | n() |
| :----------- | ----: | -----: | ---------: | --------: | ---------: | ---------: | ---------: | --------: | ---------: | ---------: | ---------: | --------: | ---------: | ---------: | ---------: | --------: | ---------: | ---------: | ------------: | -----------: | ------------: | ------------: | --: |
| Early Spring |    10 |      1 |       4.65 |      1.32 |        6.5 |        2.6 |       0.27 |      0.32 |        0.9 |        0.1 |       0.39 |      0.11 |       0.54 |       0.22 |       0.02 |      0.02 |       0.06 |       0.01 |         28.00 |        13.75 |          44.0 |           9.0 |   6 |
| Late Spring  |    10 |      1 |       4.20 |      1.03 |        5.4 |        3.2 |       0.28 |      0.40 |        1.0 |        0.1 |       0.35 |      0.08 |       0.45 |       0.27 |       0.02 |      0.03 |       0.07 |       0.01 |         27.80 |        14.24 |          45.0 |           6.0 |   5 |
| Early Spring |   200 |      1 |       9.80 |        NA |        9.8 |        9.8 |       0.40 |        NA |        0.4 |        0.4 |       0.82 |        NA |       0.82 |       0.82 |       0.03 |        NA |       0.03 |       0.03 |         27.30 |           NA |          27.3 |          27.3 |   1 |
| Late Spring  |   200 |      1 |       2.92 |      1.06 |        4.5 |        1.8 |       0.24 |      0.26 |        0.7 |        0.1 |       0.24 |      0.09 |       0.38 |       0.15 |       0.02 |      0.02 |       0.05 |       0.01 |         21.36 |        12.94 |          38.0 |           3.8 |   5 |
| Early Spring |    10 |      2 |       2.97 |      0.33 |        3.4 |        2.5 |       0.20 |      0.20 |        0.6 |        0.1 |       0.25 |      0.03 |       0.28 |       0.21 |       0.02 |      0.01 |       0.04 |       0.01 |         21.53 |         7.99 |          28.0 |           6.2 |   6 |
| Late Spring  |    10 |      2 |       3.28 |      1.84 |        6.5 |        2.1 |       0.26 |      0.36 |        0.9 |        0.1 |       0.27 |      0.15 |       0.54 |       0.18 |       0.02 |      0.02 |       0.06 |       0.01 |         23.06 |        18.44 |          54.0 |           4.3 |   5 |
| Early Spring |   200 |      2 |       3.20 |        NA |        3.2 |        3.2 |       0.10 |        NA |        0.1 |        0.1 |       0.27 |        NA |       0.27 |       0.27 |       0.01 |        NA |       0.01 |       0.01 |         27.00 |           NA |          27.0 |          27.0 |   1 |
| Late Spring  |   200 |      2 |       2.60 |      0.48 |        3.1 |        2.0 |       0.22 |      0.27 |        0.7 |        0.1 |       0.22 |      0.04 |       0.26 |       0.17 |       0.02 |      0.02 |       0.05 |       0.01 |         18.08 |         8.64 |          26.0 |           4.4 |   5 |

GF75 Blank Filters

#### Non-Blank Filter Statistics

``` r
nonblank.data <- bigelow %>% 
  filter(!Bottle == "Blank") %>% 
  mutate(BactC_µM = round(BactC_µg/12,2),
         BactN_µM = round(BactN_µg/14,2),
         BactC_N = ifelse(BactN_µM > 0, round(BactC_µM/BactN_µM,1), NA))

nonblank.stats <- nonblank.data %>% 
  mutate(Type = ifelse(Treatment %in% c("Control", "Niskin"), "Control", "Non-Control"),
         Time = ifelse(Timepoint == 0, "Initial", "Stationary")) %>% 
  group_by(Season, Depth, Type, Time, Filter) %>% 
  select(Season, Depth, Type, Time, Filter, BactC_µg, BactN_µg, BactC_µM:BactC_N) %>% 
  summarize(ave_c_µg = mean(BactC_µg, na.rm = T),
            sd_c_µg = sd(BactC_µg, na.rm = T),
            max_c_µg = max(BactC_µg, na.rm = T),
            min_c_µg = min(BactC_µg, na.rm = T),
            ave_n_µg = mean(BactN_µg, na.rm = T),
            sd_n_µg = sd(BactN_µg, na.rm = T),
            max_n_µg = max(BactN_µg, na.rm = T),
            min_n_µg = min(BactN_µg, na.rm = T),
            
            ave_c_µM = mean(BactC_µM, na.rm = T),
            sd_c_µM = sd(BactC_µM, na.rm = T),
            max_c_µM = max(BactC_µM, na.rm = T),
            min_c_µM = min(BactC_µM, na.rm = T),
            ave_n_µM = mean(BactN_µM, na.rm = T),
            sd_n_µM = sd(BactN_µM, na.rm = T),
            max_n_µM = max(BactN_µM, na.rm = T),
            min_n_µM = min(BactN_µM, na.rm = T),
            
            ave_c_n_µM = mean(BactC_N, na.rm = T),
            sd_c_n_µM = sd(BactC_N, na.rm = T),
            max_c_n_µM = max(BactC_N, na.rm = T),
            min_c_n_µM = min(BactC_N, na.rm = T),
            n()) %>% 
  mutate_at(vars(ave_c_µg:min_c_n_µM), round, 2) %>% 
  arrange(Filter, Depth, Time, Type, Season) 
```

| Season       | Depth | Type        | Time       | Filter | ave\_c\_µg | sd\_c\_µg | max\_c\_µg | min\_c\_µg | ave\_n\_µg | sd\_n\_µg | max\_n\_µg | min\_n\_µg | ave\_c\_µM | sd\_c\_µM | max\_c\_µM | min\_c\_µM | ave\_n\_µM | sd\_n\_µM | max\_n\_µM | min\_n\_µM | ave\_c\_n\_µM | sd\_c\_n\_µM | max\_c\_n\_µM | min\_c\_n\_µM | n() |
| :----------- | ----: | :---------- | :--------- | -----: | ---------: | --------: | ---------: | ---------: | ---------: | --------: | ---------: | ---------: | ---------: | --------: | ---------: | ---------: | ---------: | --------: | ---------: | ---------: | ------------: | -----------: | ------------: | ------------: | --: |
| Early Spring |    10 | Control     | Initial    |      1 |       6.76 |      1.81 |       10.0 |        4.8 |       0.79 |      0.48 |        1.4 |        0.2 |       0.56 |      0.15 |       0.83 |       0.40 |       0.06 |      0.04 |       0.10 |       0.01 |         15.91 |        13.00 |          40.0 |           6.7 |   7 |
| Late Spring  |    10 | Control     | Initial    |      1 |      15.38 |      4.40 |       19.7 |        8.7 |       2.50 |      0.83 |        3.9 |        1.8 |       1.28 |      0.37 |       1.64 |       0.72 |       0.18 |      0.06 |       0.28 |       0.13 |          7.24 |         1.60 |           8.8 |           5.5 |   5 |
| Early Autumn |    10 | Control     | Initial    |      1 |      15.41 |      8.10 |       29.0 |        5.3 |       1.79 |      1.07 |        3.2 |        0.4 |       1.28 |      0.68 |       2.42 |       0.44 |       0.13 |      0.08 |       0.23 |       0.03 |         12.94 |         9.87 |          37.3 |           5.0 |  10 |
| Early Spring |    10 | Control     | Stationary |      1 |      14.61 |      5.24 |       23.6 |        8.6 |       2.64 |      1.49 |        5.5 |        1.0 |       1.22 |      0.44 |       1.97 |       0.72 |       0.19 |      0.11 |       0.39 |       0.07 |          7.43 |         1.93 |          11.0 |           3.5 |  14 |
| Late Spring  |    10 | Control     | Stationary |      1 |      17.35 |      8.64 |       35.8 |        9.5 |       3.05 |      1.47 |        5.8 |        0.8 |       1.45 |      0.72 |       2.98 |       0.79 |       0.22 |      0.10 |       0.41 |       0.06 |          9.56 |        12.74 |          49.7 |           3.7 |  12 |
| Early Autumn |    10 | Control     | Stationary |      1 |      17.05 |      6.22 |       26.2 |        6.5 |       2.08 |      0.92 |        3.1 |        0.4 |       1.42 |      0.52 |       2.18 |       0.54 |       0.15 |      0.06 |       0.22 |       0.03 |         11.49 |         7.24 |          31.7 |           6.5 |  12 |
| Early Spring |    10 | Non-Control | Stationary |      1 |      76.03 |     35.88 |      122.3 |       36.3 |      18.07 |      8.50 |       28.6 |        8.5 |       6.34 |      2.99 |      10.19 |       3.02 |       1.29 |      0.61 |       2.04 |       0.61 |          4.92 |         0.10 |           5.0 |           4.8 |   6 |
| Late Spring  |    10 | Non-Control | Stationary |      1 |      67.60 |     22.88 |       93.0 |       45.3 |      17.68 |      6.28 |       23.3 |       11.7 |       5.63 |      1.91 |       7.75 |       3.77 |       1.26 |      0.45 |       1.66 |       0.84 |          4.50 |         0.28 |           4.7 |           4.1 |   4 |
| Early Autumn |    10 | Non-Control | Stationary |      1 |      43.48 |     22.26 |       89.2 |       22.3 |       9.31 |      5.45 |       19.2 |        4.5 |       3.62 |      1.85 |       7.43 |       1.86 |       0.66 |      0.39 |       1.37 |       0.32 |          5.72 |         1.40 |           9.1 |           4.1 |   9 |
| Early Spring |   200 | Control     | Initial    |      1 |       9.50 |        NA |        9.5 |        9.5 |       0.30 |        NA |        0.3 |        0.3 |       0.79 |        NA |       0.79 |       0.79 |       0.02 |        NA |       0.02 |       0.02 |         39.50 |           NA |          39.5 |          39.5 |   1 |
| Late Spring  |   200 | Control     | Initial    |      1 |       4.98 |      1.10 |        6.5 |        3.6 |       0.42 |      0.45 |        1.2 |        0.1 |       0.42 |      0.09 |       0.54 |       0.30 |       0.03 |      0.03 |       0.09 |       0.01 |         24.06 |        18.13 |          46.0 |           6.0 |   5 |
| Early Autumn |   200 | Control     | Initial    |      1 |      12.10 |        NA |       12.1 |       12.1 |       1.10 |        NA |        1.1 |        1.1 |       1.01 |        NA |       1.01 |       1.01 |       0.08 |        NA |       0.08 |       0.08 |         12.60 |           NA |          12.6 |          12.6 |   1 |
| Late Spring  |   200 | Control     | Stationary |      1 |       6.57 |      1.21 |        8.6 |        5.2 |       0.90 |      0.54 |        1.6 |        0.1 |       0.55 |      0.10 |       0.72 |       0.43 |       0.06 |      0.04 |       0.11 |       0.01 |         17.99 |        22.24 |          72.0 |           4.8 |  10 |
| Early Autumn |   200 | Control     | Stationary |      1 |      18.35 |      2.19 |       19.9 |       16.8 |       3.65 |      0.21 |        3.8 |        3.5 |       1.53 |      0.18 |       1.66 |       1.40 |       0.26 |      0.01 |       0.27 |       0.25 |          5.85 |         0.35 |           6.1 |           5.6 |   2 |
| Late Spring  |   200 | Non-Control | Stationary |      1 |      47.30 |     22.60 |       68.0 |       27.4 |      11.30 |      6.70 |       17.3 |        5.5 |       3.94 |      1.89 |       5.67 |       2.28 |       0.81 |      0.48 |       1.24 |       0.39 |          5.22 |         0.79 |           6.0 |           4.4 |   4 |
| Early Autumn |   200 | Non-Control | Stationary |      1 |      36.63 |     21.32 |       60.7 |       20.1 |       7.37 |      6.52 |       14.9 |        3.5 |       3.06 |      1.78 |       5.06 |       1.68 |       0.52 |      0.46 |       1.06 |       0.25 |          6.93 |         2.26 |           9.3 |           4.8 |   3 |
| Early Spring |    10 | Control     | Initial    |      2 |       3.47 |      0.55 |        4.2 |        2.7 |       0.27 |      0.19 |        0.6 |        0.1 |       0.29 |      0.05 |       0.35 |       0.22 |       0.02 |      0.01 |       0.04 |       0.01 |         18.71 |        10.13 |          35.0 |           8.2 |   7 |
| Late Spring  |    10 | Control     | Initial    |      2 |       4.60 |      1.44 |        6.0 |        2.3 |       0.50 |      0.38 |        0.9 |        0.1 |       0.38 |      0.12 |       0.50 |       0.19 |       0.03 |      0.03 |       0.06 |       0.01 |         21.84 |        19.42 |          47.0 |           3.2 |   5 |
| Early Autumn |    10 | Control     | Initial    |      2 |      11.06 |     10.14 |       37.0 |        4.0 |       0.70 |      0.44 |        1.3 |        0.1 |       0.92 |      0.84 |       3.08 |       0.33 |       0.05 |      0.03 |       0.09 |       0.01 |         27.18 |        25.72 |          83.0 |           4.0 |  10 |
| Early Spring |    10 | Control     | Stationary |      2 |       4.37 |      1.37 |        7.4 |        2.9 |       0.28 |      0.17 |        0.5 |        0.1 |       0.36 |      0.11 |       0.62 |       0.24 |       0.02 |      0.01 |       0.04 |       0.01 |         21.74 |         9.33 |          37.0 |           9.0 |  14 |
| Late Spring  |    10 | Control     | Stationary |      2 |       3.84 |      0.76 |        5.4 |        2.7 |       1.43 |      1.61 |        6.2 |        0.0 |       0.32 |      0.06 |       0.45 |       0.22 |       0.10 |      0.11 |       0.44 |       0.00 |          4.02 |         2.29 |           8.5 |           0.6 |  12 |
| Early Autumn |    10 | Control     | Stationary |      2 |       8.60 |      5.41 |       20.8 |        3.8 |       0.74 |      0.72 |        2.6 |        0.1 |       0.72 |      0.45 |       1.73 |       0.32 |       0.05 |      0.05 |       0.19 |       0.01 |         29.71 |        27.29 |          78.0 |           3.6 |  12 |
| Early Spring |    10 | Non-Control | Stationary |      2 |       4.93 |      0.83 |        6.1 |        3.6 |       0.90 |      1.67 |        4.3 |        0.1 |       0.41 |      0.07 |       0.51 |       0.30 |       0.06 |      0.12 |       0.31 |       0.01 |         25.93 |        16.59 |          43.0 |           1.4 |   6 |
| Late Spring  |    10 | Non-Control | Stationary |      2 |       5.70 |      1.96 |        7.4 |        3.5 |       0.85 |      0.54 |        1.5 |        0.2 |       0.48 |      0.17 |       0.62 |       0.29 |       0.06 |      0.04 |       0.11 |       0.01 |         15.43 |        15.57 |          38.0 |           2.6 |   4 |
| Early Autumn |    10 | Non-Control | Stationary |      2 |       9.06 |      4.81 |       16.3 |        3.9 |       0.81 |      0.90 |        2.4 |        0.1 |       0.76 |      0.40 |       1.36 |       0.32 |       0.06 |      0.06 |       0.17 |       0.01 |         36.41 |        44.64 |         132.0 |           6.1 |   9 |
| Early Spring |   200 | Control     | Initial    |      2 |       7.00 |        NA |        7.0 |        7.0 |       0.10 |        NA |        0.1 |        0.1 |       0.58 |        NA |       0.58 |       0.58 |       0.01 |        NA |       0.01 |       0.01 |         58.00 |           NA |          58.0 |          58.0 |   1 |
| Late Spring  |   200 | Control     | Initial    |      2 |       2.58 |      0.52 |        3.4 |        2.0 |       0.22 |      0.27 |        0.7 |        0.1 |       0.21 |      0.04 |       0.28 |       0.17 |       0.02 |      0.02 |       0.05 |       0.01 |         18.68 |         9.15 |          28.0 |           3.4 |   5 |
| Early Autumn |   200 | Control     | Initial    |      2 |       4.20 |        NA |        4.2 |        4.2 |       0.10 |        NA |        0.1 |        0.1 |       0.35 |        NA |       0.35 |       0.35 |       0.01 |        NA |       0.01 |       0.01 |         35.00 |           NA |          35.0 |          35.0 |   1 |
| Late Spring  |   200 | Control     | Stationary |      2 |       2.85 |      0.65 |        3.7 |        1.5 |       0.44 |      0.32 |        0.9 |        0.1 |       0.24 |      0.06 |       0.31 |       0.12 |       0.03 |      0.02 |       0.06 |       0.01 |         11.85 |         8.66 |          27.0 |           3.8 |  10 |
| Early Autumn |   200 | Control     | Stationary |      2 |       6.55 |      2.47 |        8.3 |        4.8 |       0.40 |      0.28 |        0.6 |        0.2 |       0.54 |      0.21 |       0.69 |       0.40 |       0.02 |      0.02 |       0.04 |       0.01 |         39.50 |        41.72 |          69.0 |          10.0 |   2 |
| Late Spring  |   200 | Non-Control | Stationary |      2 |       4.12 |      0.35 |        4.3 |        3.6 |       0.57 |      0.62 |        1.4 |        0.1 |       0.34 |      0.03 |       0.36 |       0.30 |       0.04 |      0.04 |       0.10 |       0.01 |         19.20 |        16.19 |          36.0 |           3.6 |   4 |
| Early Autumn |   200 | Non-Control | Stationary |      2 |       6.93 |      4.63 |       12.2 |        3.5 |       0.53 |      0.75 |        1.4 |        0.1 |       0.58 |      0.39 |       1.02 |       0.29 |       0.04 |      0.05 |       0.10 |       0.01 |         27.07 |        15.99 |          42.0 |          10.2 |   3 |

GF75 Non-blank Filters

``` r
all_filter_stats <- blank.stats %>%
  mutate(Sample = "Blank",
         Type = "Blank",
         Time = "Initial") %>%
  bind_rows(., nonblank.stats %>% mutate(Sample = "Non-blank")) %>%
  select(Sample, Type, everything()) 
  
#write_csv(all_filter_stats, "~/Desktop/NAAMES_GF75_Stats.csv")
```

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-53-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-54-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-55-1.png" style="display: block; margin: auto;" />

This plot shows greater variability in the POC values from the first
blank filter (red points). Whether using surface or deep water, 0.2 µm
filtrate or TFF filtrate, the small range in the POC values of the
second filter indicates that the second filter could provide a
reasonable estimate of DOC adsorption onto GF75s (at least for the
spring cruises) and represent a decent blank correction.

### Blank Corrections and Bacterial Carbon Estimates

We’ll report non-blank corrected POC estimates as well as
blank-corrected POC estimates. Because we will be calculating delta
bacterial carbon to then estimate BGE, a blank correction may not be
necessary. However, if we were to blank correct, we’ll use universal
correction factors based on the data/plots above. The detection limit
for Bigelow’s Costech ECS 4010 elemental analysis is 0.1 µg C.

We’ll make 3 universal correction factors:

  - poc.cf1 = the mean of filter 2 for all blanks
  - poc.cf2 = the mean of filter 2 of only TFF blanks
  - poc.cf3 = the mean of filter 1 and 2 of only TFF blanks

Each of these will be subtracted from each sample filter to estimate
bacterial carbon in µg. Values below the detection limits for N (\<0.1)
are converted to 0. The addition of the corrected value from filter 1
and from filter 2 of each sample will provide the bacterial carbon
estimate for that sample.

Bacterial carbon estimates will then be converted to µmol
C<sup>-1</sup>. The same will be done for PON.

``` r
corrfactor <- blank.data %>%
  group_by(Filter) %>%
  mutate(poc.cf1 = ifelse(Filter == 2, round(mean(BactC_µg),1), NA),
         pon.cf1 = ifelse(Filter == 2, round(mean(BactN_µg),1), NA)) %>% 
  ungroup() %>% 
  group_by(BlnkFiltrate, Filter) %>%
  mutate(poc.cf2 = ifelse(BlnkFiltrate == "TFF" & Filter == 2, round(mean(BactC_µg),1), NA),
         pon.cf2 = ifelse(BlnkFiltrate == "TFF" & Filter == 2, round(mean(BactN_µg),1), NA)) %>% 
  ungroup() %>% 
  group_by(BlnkFiltrate) %>% 
  mutate(poc.cf3 = ifelse(BlnkFiltrate == "TFF", round(mean(BactC_µg),1), NA),
         pon.cf3 = ifelse(BlnkFiltrate == "TFF", round(mean(BactN_µg),1), NA)) %>% 
  ungroup() %>% 
  fill(poc.cf1:pon.cf3, .direction = "updown") %>% 
  select(poc.cf1:pon.cf3) %>% 
  distinct()

corrfactor
```

    ## # A tibble: 1 x 6
    ##   poc.cf1 pon.cf1 poc.cf2 pon.cf2 poc.cf3 pon.cf3
    ##     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1       3     0.2       3     0.2     4.2     0.2

The carbon correction factors for the mean of filter 2 for all blanks
and the mean of filter 2 of only TFF blanks are the same. So, we’ll
calculate 2 corrected values:

  - poc.c1.ug = POC (µg) - poc.cf1 or 2
  - poc.c3.ug = POC (µg) - poc.cf3

and convert the calculated vales to µmol C L<sup>-1</sup>:

  - poc.c1.um
  - poc.c3.um

Since the nitrogen carbon conversion factors are all the same, we’ll
only calculate 1 corrected value

  - pon.c.ug = PON (µg) - pon.cf1, 2, or 3

and convert the calculated values to µmol N L <sup>-1</sup>:

  - pon.c.um

We’ll also calculate C:N ratios:

  - cn = poc.um/pon.um
  - cn.c1 = poc.c1.um/pon.um
  - cn.c3 = poc.c3.um/pon.m

<!-- end list -->

``` r
bactcarbon <- bigelow %>%
  filter(!Bottle == "Blank") %>%
  mutate(Treatment = ifelse(Bottle %in% c("1L", "100mls", "250mls", "500mls"), "Volume", Treatment)) %>% 
  rename(pon.ug = BactN_µg,
         poc.ug = BactC_µg) %>% 
  cbind(., corrfactor) %>% 
  mutate(poc.c1.ug = poc.ug - poc.cf1,
         poc.c3.ug = poc.ug - poc.cf3,
         pon.c.ug = pon.ug - pon.cf1,
         poc.c1.ug = ifelse(poc.c1.ug < 0, 0, poc.c1.ug),
         poc.c3.ug = ifelse(poc.c3.ug < 0, 0, poc.c3.ug),
         pon.c.ug = ifelse(pon.c.ug < 0, 0, pon.c.ug)) %>%
  group_by(Cruise, Station, Depth, Bottle, Timepoint) %>%
  mutate(poc.ug = round((sum(poc.ug)),1),
         pon.ug = round((sum(pon.ug)),1),
         poc.c1.ug = round((sum(poc.c1.ug)),1),
         poc.c3.ug = round((sum(poc.c3.ug)),1),
         pon.c.ug = round((sum(pon.c.ug)),1)) %>% 
  ungroup() %>% 
  mutate(poc.um = round(poc.ug/12, 1),
         poc.c1.um = round(poc.c1.ug/12, 1),
         poc.c3.um = round(poc.c3.ug/12, 1),
         pon.um = round(pon.ug/14, 1),
         pon.c.um = round(pon.c.ug/14, 1),
         cn = round(poc.um/pon.um),
         cn.c1 = round(poc.c1.um/pon.um),
         cn.c3 = round(poc.c3.um/pon.um)) %>% 
  select(Season, Cruise, Station, Depth, facetlabel, Treatment, Bottle, Timepoint, poc.ug, pon.ug, poc.cf1, poc.cf2, poc.cf3, pon.cf1, pon.cf2, pon.cf3, poc.c1.ug, poc.c3.ug, pon.c.ug, poc.um:cn.c3 ) %>% 
  distinct() %>% 
  arrange(Season, Cruise, Station, Depth, Treatment, Bottle, Timepoint) 
```

### Sampling Volume Effect on Bacterial Carbon

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-58-1.png" style="display: block; margin: auto;" />

This test needs to be done again \#flattenthecurve

### Calculate Carbon Per Cell

Here we calculate fg C or fg N cell<sup>-1</sup> by:

  - Estimating the cell abundance captured on the 0.3 µm GF75 filters
    (whole water - filtrate, cells L<sup>-1</sup>)
  - Dividing the GF75 POC or PON value (µmol C or N L<sup>-1</sup>) by
    the cell abundance on the GF75 filter, then applying the subsequent
    conversions

The GF75 POC data are tidied before calculations are made (e.g renaming
bottle identifiers so they are congruent with datasets to be merged
with, duplicating values for initial timepoints where one POC value was
taken for replicate bottles). For ease, initial POC and cell abundance
data are handled separately from data collected at stationary
timepoints.

#### Calculate Filter Cell Abundance and Retention

``` r
#initial cell abundance
i.gf75 <- ba.df %>% 
  select(Season:cells) %>% 
  filter(Timepoint == 0, Bottle %in% c("GF75", "Niskin")) %>% 
  group_by(Cruise, Station, Depth) %>% 
  mutate(filtrate.cells = ifelse(Bottle == "GF75", cells, NA)) %>% 
  fill(filtrate.cells, .direction = "downup") %>% 
  filter(Bottle == "Niskin") %>% 
  mutate(i.gf75.cells = cells - filtrate.cells,
         i.gf75.ret = round((cells - filtrate.cells)/cells, 2)) %>% 
  ungroup() %>% 
  select(Season:Depth, i.gf75.cells, i.gf75.ret)

#stationary cell abundance
s.gf75 <- ba.df %>% 
  select(Season:cells) %>% 
  filter(!Timepoint == 0) %>% 
  arrange(Cruise, Station, Depth, Bottle, Timepoint) %>% 
  group_by(Cruise, Station,  Depth, Bottle) %>%
  mutate(gf75.timepoint = ifelse(Treatment == "GF75", Timepoint, NA)) %>% 
  fill(gf75.timepoint, .direction = "updown") %>% 
  filter(Timepoint == gf75.timepoint) %>% 
  mutate(filtrate.cells = ifelse(Treatment == "GF75", cells, NA)) %>% 
  fill(filtrate.cells, .direction = "updown") %>%  
  filter(!Treatment == "GF75") %>% 
  mutate(s.gf75.cells = cells - filtrate.cells,
         s.gf75.ret = round((cells - filtrate.cells)/cells, 2)) %>% 
  ungroup() %>% 
  select(Season:Bottle, Treatment, Timepoint, s.gf75.cells, s.gf75.ret) 
```

### Combine Cell Abundance and POC Data, Calculate CCF

We’ll calculate 4 CCFs (fg C cell<sup>-1</sup>) from the POC data:

  - i.ccf.c1 = poc.c1.um/cells on filter \* 12 \* 10<sup>9</sup>
    (initial)
  - s.ccf.c1 = poc.c1.um/cells on filter \* 12 \* 10<sup>9</sup>
    (stationary)
  - i.ccf.c3 = poc.c3.um/cells on filter \* 12 \* 10<sup>9</sup>
    (initial)
  - s.ccf.c3 = poc.c3.um/cells on filter \* 12 \* 10<sup>9</sup>
    (stationary)

We’ll also add 3 other CCFs from the literature:

  - ccf.white = 6.5 (SAR11), White et al 2019
      - also close to 6.3 of Carlson et al 1996 and 7 of Zubkov 2000
  - ccf.fukuda = 12.3 (average oceanic Pacific), Fukuda et al 1998
      - also close to 10 generally used by many papers Wear et al 2020,
        Christian and Karl 1994, Caron et al 1995, James et al 2017
  - ccf.lee = 20 (Long Island Sound), Lee and Fuhrman 1987
      - generally used as reported by Gunderson et al 2002

<!-- end list -->

``` r
#initial timepoints
i.fg_cell <-  bactcarbon %>% 
  filter(!Treatment == "Volume", Timepoint == 0) %>% 
  select_at(vars(poc.ug:cn.c3), .funs = funs(paste0("i.", .))) %>% 
  bind_cols(bactcarbon %>% filter(!Treatment == "Volume", Timepoint == 0) %>% select(Season:facetlabel), .) %>% 
  mutate(Station = gsub("2RD", "S2RD", Station),
         Station = gsub("2RF", "S2RF", Station)) %>% 
  left_join(., i.gf75) %>% 
  mutate(i.ccf.c1 = round(i.poc.c1.um/i.gf75.cells * 12 * 10^9),
         i.ccf.c3 = round(i.poc.c3.um/i.gf75.cells * 12 * 10^9)) 
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

``` r
#stationary timepoints
s.fg_cell <- bactcarbon %>% 
  filter(!Treatment == "Volume", !Timepoint == 0) %>% 
  select_at(vars(poc.ug:cn.c3), .funs = funs(paste0("s.", .))) %>% 
  bind_cols(bactcarbon %>% filter(!Treatment == "Volume", !Timepoint == 0) %>% select(Season:facetlabel, Treatment, Bottle, Timepoint), .) %>% 
  mutate(Station = gsub("2RD", "S2RD", Station),
         Station = gsub("2RF", "S2RF", Station),
         Treatment = gsub("TW15", "TW20", Treatment),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW12" & Bottle == "A12", gsub("A12", "C", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW12" & Bottle == "B12", gsub("B12", "D", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW13" & Bottle == "A13", gsub("A13", "E", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 10 & Treatment == "TW13" & Bottle == "B13", gsub("B13", "F", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "Control" & Bottle == "C", gsub("C", "G", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "Control" & Bottle == "D", gsub("D", "H", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW12" & Bottle == "C12", gsub("C12", "I", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW12" & Bottle == "D12", gsub("D12", "J", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW13" & Bottle == "C13", gsub("C13", "K", Bottle), Bottle),
         Bottle = ifelse(Cruise == "AT34" & Station == 2 & Depth == 200 & Treatment == "TW13" & Bottle == "D13", gsub("D13", "L", Bottle), Bottle)) %>% 
  add_row(Season = "Early Autumn", Cruise = "AT38", Station = 1, Depth = 10, Treatment = "SynExd", Bottle = "I") %>% 
    add_row(Season = "Early Autumn", Cruise = "AT38", Station = 1, Depth = 10, Treatment = "SynLys", Bottle = "N") %>% 
  add_row(Season = "Early Autumn", Cruise = "AT38", Station = 1, Depth = 10, Treatment = "TWExd", Bottle = "K") %>% 
  add_row(Season = "Early Autumn", Cruise = "AT38", Station = 1, Depth = 200, Treatment = "SynLys", Bottle = "P") %>% 
  arrange(Cruise, Station, Depth, Treatment, Bottle) %>% 
  group_by(Cruise, Station, Depth, Treatment) %>% 
  fill(s.poc.ug:s.cn.c3, .direction = "updown") %>% 
  mutate(Bottle = gsub("JI", "J", Bottle),
         Bottle = gsub("MN", "M", Bottle),
         Bottle = gsub("LK", "L", Bottle),
         Bottle = gsub("OP", "O", Bottle)) %>% 
  left_join(., s.gf75) %>% 
  mutate(s.ccf.c1 = round(s.poc.c1.um/s.gf75.cells * 12 * 10^9),
         s.ccf.c3 = round(s.poc.c3.um/s.gf75.cells * 12 * 10^9)) %>% 
  ungroup() %>% 
  rename(s.timepoint = Timepoint)

fg_cell <- left_join(s.fg_cell, i.fg_cell)  %>% 
   group_by(Cruise, Station, Depth, Treatment) %>% 
  fill(c(s.timepoint,i.poc.ug:i.ccf.c3), .direction = "updown") %>% 
  mutate(ccf.white = 6.5,
         ccf.fukuda = 12.3,
         ccf.lee = 20) %>% 
  ungroup()

fg_cell$Season <- factor(fg_cell$Season, levels = levels)
```

#### Retention of GF75 Filters

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-61-1.png" style="display: block; margin: auto;" />

``` r
fg_cell %>% 
  select(Season, Cruise, Station, Bottle, Depth, Treatment, i.gf75.ret, s.gf75.ret) %>% 
  mutate(facet_depth = paste(Depth, "m")) %>% 
  gather(key, value, i.gf75.ret:s.gf75.ret) %>% 
  mutate(key = gsub("i.gf75.ret", "Initial", key),
         key = gsub("s.gf75.ret", "Stationary", key),
         value = ifelse(key == "Initial" & Bottle %in% c("B", "D", "F"), NA, value)) %>% 
  filter(key == "Initial", Treatment == "Control") %>% 
  drop_na() %>% 
  ungroup() %>% 
  rename(Initial_Summary = value) %>% 
  select(Initial_Summary) %>% 
  summary()
```

    ##  Initial_Summary 
    ##  Min.   :0.1000  
    ##  1st Qu.:0.5575  
    ##  Median :0.7350  
    ##  Mean   :0.6567  
    ##  3rd Qu.:0.8000  
    ##  Max.   :0.8800

This figure shows that while GF75 filter retention is fairly consistent
across the samples taken during the stationary phase of cell growth
(median 85%, mean 83%), they are pretty variable for the samples taken
at the initial condition (median 74%, mean 66%, range 10 - 88%, 1st Qu =
~56%). Those samples with filter retentions \<= 0.58 (1st Qu.) likely
reflect faulty procedure (i.e. ripped filters) and may be associated
with questionable measurements. We’ll flag these data.

``` r
fg_cell.qc <- fg_cell %>% 
  mutate(i.gf75.flag = ifelse(i.gf75.ret <= 0.56, "Questionable", "Acceptable"),
         i.gf75.flag = ifelse(is.na(i.gf75.flag), "ND", i.gf75.flag),
         s.gf75.flag = ifelse(s.gf75.ret <= 0.56, "Questionable", "Acceptable"),
         s.gf75.flag = ifelse(is.na(s.gf75.flag), "ND", s.gf75.flag),
         gf75.flag = ifelse(i.gf75.flag == "Questionable", "Initial Questionable", NA),
         gf75.flag = ifelse(s.gf75.flag == "Questionable" & is.na(gf75.flag), "Stationary Questionable", gf75.flag),
         gf75.flag = ifelse(is.na(gf75.flag) & is.na(i.gf75.flag) & is.na(s.gf75.flag), "ND", gf75.flag),
         gf75.flag = ifelse(is.na(gf75.flag) & i.gf75.flag == "ND" & s.gf75.flag == "ND", "ND", gf75.flag),
         gf75.flag = ifelse(is.na(gf75.flag), "Acceptable", gf75.flag),
         gf75.flag = ifelse(i.gf75.flag == "ND", paste(gf75.flag, "ND Initial", sep = "/" ), gf75.flag),
          gf75.flag = ifelse(s.gf75.flag == "ND", paste(gf75.flag, "ND Stationary", sep = "/" ), gf75.flag)) %>% 
  select(-c(i.gf75.flag, s.gf75.flag))
```

#### CCFs and C:N ratios

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-64-1.png" style="display: block; margin: auto;" />

There’s a lot going on this plot. “Questionable” refers to samples with
retentions \>= 56%, where cell slippage may have occured. “ND” refers to
points where there is missing filter retention data. We’ve faceted the
plot by (1) control or non-control, (2) the type of blank correction,
and (3) depth.

Overall, there are mostly questionable measurements at 200 m. It is
peculiar that the acceptable measurements at that depth horizon have
CCFs \>20. There are a lot of acceptable measurments at 10 m, with
outliers at the initial condition being questionalble measurements. The
CCFs for the Early Autumn still seem high even excluding the outliers.

The POC correction using both filters of the TFF blanks appears to
dampen the
CCFs.

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-65-1.png" style="display: block; margin: auto;" />

C:N ratios of POC-corrected samples appear to be more realistic than
uncorrected values. However, even exluding questionable measurements, N3
still has an outlier for the initial condition. The points on the line
intersecting the max of the y-axis represent “Inf” values, where PON
values are 0.

Summary Statistics including N3: i.cn i.cn.c1 i.cn.c3  
Min. : 4.0 Min. :2.0 Min. : 1  
1st Qu.: 8.0 1st Qu.:5.0 1st Qu.: 3  
Median :11.5 Median :8.5 Median : 6  
Mean : Inf Mean :Inf Mean :Inf  
3rd Qu.: Inf 3rd Qu.:Inf 3rd Qu.: 21  
Max. : Inf Max. :Inf Max. :Inf  
NA’s :2 NA’s :2 NA’s :8

``` 
i.ccf.c1        i.ccf.c3    
```

Min. : 5.00 Min. : 0.00  
1st Qu.: 8.25 1st Qu.: 4.00  
Median :16.00 Median : 9.00  
Mean :25.65 Mean :19.09  
3rd Qu.:41.50 3rd Qu.:35.00  
Max. :77.00 Max. :67.00  
NA’s :4 NA’s :4

Summary Statistics excluding N3: i.cn i.cn.c1 i.cn.c3  
Min. : 4 Min. : 2 Min. :1.0  
1st Qu.: 8 1st Qu.: 3 1st Qu.:2.0  
Median : 10 Median : 6 Median :4.5  
Mean :Inf Mean :Inf Mean :Inf  
3rd Qu.:Inf 3rd Qu.:Inf 3rd Qu.:6.0  
Max. :Inf Max. :Inf Max. :Inf  
NA’s :2 NA’s :2 NA’s :8

``` 
i.ccf.c1        i.ccf.c3     
```

Min. : 5.00 Min. : 0.000  
1st Qu.: 7.00 1st Qu.: 3.000  
Median :11.00 Median : 7.000  
Mean :14.94 Mean : 8.188  
3rd Qu.:16.75 3rd Qu.: 9.250  
Max. :52.00 Max. :35.000  
NA’s :4 NA’s :4

Overall (including AT38), the median corrected bacterial C:N ratios for
the in situ community (initial condition) estimated here are are
comparable to those reported for oceanic bacteria by Fukuda et al.
(1998): \[6.8 \pm 1.2\]

While the medians are within the range, the mean aveage contents of
carbon are higher than those reported by Fukuda et al. (1998):
\[12.4 \pm 6.3\] fg C cell<sup>-1</sup>

Deviations from Fukuda et al. may be in part due to differences in
sampling locations (NA v Pacific and Southern Ocean) and time of
collection (range of seasons for NAAMES, winter in both Pacific and
Southern Ocean). They could also be due to differences in growth
conditions (i.e nutrient limitation) and differences in bacterial
communities (Vrede et al 2002). The estimates become more in agreement
with Fukuda et al when N3 data are
excluded.

#### Cell Carbon v GF75 Retention

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-70-1.png" style="display: block; margin: auto;" />

#### Cell Carbon v GF75 Filter Cell Abundance

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-71-1.png" style="display: block; margin: auto;" />

#### Cell Carbon v POC

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-72-1.png" style="display: block; margin: auto;" />

#### Cell Carbon v C:N

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-73-1.png" style="display: block; margin: auto;" />

## Biovolume

Cell biovolume was determined from 10 images captured with a digital
camera (Retiga Exi-QImaging, Surrey, BC, Canada). Images were processed
using ImageJ software according to established image analyses protocols
(Baldwin and Bankston 1988; Sieracki et al. 1989, Ducklow et al., 1995).
In addition to calculating cell carbon using GF75 POC-derived CCFs,
we’ll also caclulate cell carbon based on the biovolume measurements
we have. Biovolume-based cell carbon could be be calculated as:
**Bacterial Abundance \* Cell Biovolume \* CCF of 148 fg C
µm<sup>-3</sup>** However, this CCF, which comes from Gunderson et al,
is likely only really applicable to the Sargasso Sea. We will need to
determine the best way to estimate biovolume-based cell carbon for this
dataset.

``` r
biovol <- read_csv("~/naames_bioav_ms/Input/Biovolume_Summary.csv") %>% 
  select(Cruise:Timepoint, Biovol_Mean_All) %>% 
  rename(biovol = Biovol_Mean_All) %>% 
  mutate_at(vars(biovol), round, 3) %>% 
  mutate(phase = ifelse(Timepoint == 0, "initial", "later"),
         Station = gsub("2RD", "S2RD", Station),
         Station = gsub("2RF", "S2RF", Station)) %>% 
  group_by(Cruise, phase) %>% 
  mutate(mean_biovol = round(mean(biovol, na.rm = T),3),
         sd_biovol = round(sd(biovol, na.rm = T),3)) %>% 
  add_tally() %>% 
  rename(n_biovol = n) %>% 
  ungroup() 

biovol_table <- biovol %>% select(-c(Station:biovol)) %>% distinct()
```

| Cruise | phase   | mean\_biovol | sd\_biovol | n\_biovol |
| :----- | :------ | -----------: | ---------: | --------: |
| AT38   | initial |        0.032 |      0.008 |         5 |
| AT38   | later   |        0.040 |      0.027 |        14 |
| AT39   | initial |        0.055 |      0.023 |         6 |
| AT39   | later   |        0.054 |      0.021 |        12 |

Mean Biovolume AT38, AT39; ‘later’ represents any timepoint after 0

It would be good to see if there is a log-linear relationship between
the cell carbon calculated from the GF75 POC data (fg C L<sup>-1</sup>)
and the estimated biovolumes (µm<sup>3</sup>). Unfortunately we don’t
yet have concurrent measurements of GF75 POC data and biovolume at
“stationary” so that will have to wait. We do have concurrent
measurements at initial though, so we can plot those.

``` r
biovol_merge <- biovol %>% 
  filter(Timepoint == 0) %>% 
  select(Cruise:Bottle, biovol) %>% 
  rename(i.biovol = biovol) %>% 
  left_join(fg_cell.qc, .)
```

#### Biovolume v. Cell Carbon (this study and literature)

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-77-1.png" style="display: block; margin: auto;" />

For the samples in which we have concurrent biovolume and cellular C
measurements, we see that their relationship to one another is pretty
widespread relative to published literature relationshops. Circles are
samples corrected using filter 2 of either all POC blanks or only TFF
blanks. Diamonds are samples corrected using the sum of the filters for
TFF blanks.

# Merge Data

``` r
ba.merge <- gcstat.df %>% 
  mutate(key = paste(Cruise, Station, Depth, Treatment, Bottle, Timepoint, sep = "."),
         Depth = as.numeric(Depth)) %>% 
  filter(!Treatment %in% c("GF75", "Niskin", "TFF-Ret", "Volume", "Parallel"))

oc.merge <- oc_p %>% 
  select(-Days) %>% 
  filter(!Treatment == "Parallel")


merge <- full_join(ba.merge, oc.merge) %>% 
  left_join(., fg_cell.qc %>% select(-c(contains(".ug"), contains("gf75.cells")))) %>% 
  arrange(Cruise, Station, Depth, Treatment, Bottle, Hours) %>%
  mutate(Days = round(Hours/24, 1)) %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) %>%
  select(Season:Depth, facet_depth,  Treatment, facet_treatment, Bottle, facet_bottle, Timepoint, Hours, stationary, s.timepoint, Days, cells:sd_p_cells, doc, sd_doc, ptoc, sd_ptoc, doc_from_t0, ptoc_from_t0, gf75.flag, i.poc.cf1:i.ccf.c3, s.poc.cf1:s.ccf.c3, ccf.white:ccf.lee) %>% 
  fill(stationary, facet_treatment, facet_depth, facet_bottle, .direction = "downup") %>% 
  ungroup() %>% 
  drop_na(Hours) 
```

# Bottle v. Vial Incubation Comparisons

## PTOC v. PDOC

### Absolute

``` r
ptoc_pdoc.data <- oc_p %>% 
  filter(Cruise == "AT39", !Treatment == "Parallel") 
ptoc_pdoc.reg <- lmodel2(pdoc ~ ptoc, data = ptoc_pdoc.data , nperm = 99)
```

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = pdoc ~ ptoc, data = ptoc_pdoc.data, nperm = 99)
    ## 
    ## n = 71   r = 0.9732196   r-square = 0.9471563 
    ## Parametric P-values:   2-tailed = 8.622836e-46    1-tailed = 4.311418e-46 
    ## Angle between the two OLS regression lines = 1.554332 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  3.337716 0.9424277        43.30229              0.01
    ## 2     MA  1.863641 0.9675048        44.05379              0.01
    ## 3    SMA  1.813323 0.9683609        44.07912                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS      0.1837596        6.491672  0.8889664   0.9958891
    ## 2     MA     -1.4548804        5.004919  0.9140651   1.0239599
    ## 3    SMA     -1.4159087        4.869193  0.9163741   1.0232969
    ## 
    ## Eigenvalues: 48.72439 0.6605831 
    ## 
    ## H statistic used for computing C.I. of MA: 0.0008036195

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-81-1.png" style="display: block; margin: auto;" />

### Delta

``` r
delta_ptoc_pdoc.reg <- lmodel2(pdoc_from_t0 ~ ptoc_from_t0, data = ptoc_pdoc.data , nperm = 99)
```

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = pdoc_from_t0 ~ ptoc_from_t0, data =
    ## ptoc_pdoc.data, nperm = 99)
    ## 
    ## n = 71   r = 0.9711193   r-square = 0.9430727 
    ## Parametric P-values:   2-tailed = 1.126827e-44    1-tailed = 5.634135e-45 
    ## Angle between the two OLS regression lines = 1.678043 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method   Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  0.09589786 0.9411488        43.26345              0.01
    ## 2     MA -0.02473289 0.9682354        44.07541              0.01
    ## 3    SMA -0.02875343 0.9691382        44.10209                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS     -0.3359442       0.5277399  0.8856156    0.996682
    ## 2     MA     -0.2867641       0.2227450  0.9126664    1.027072
    ## 3    SMA     -0.2831516       0.2114846  0.9151948    1.026261
    ## 
    ## Eigenvalues: 78.77687 1.153065 
    ## 
    ## H statistic used for computing C.I. of MA: 0.0008695124

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-84-1.png" style="display: block; margin: auto;" />

## Bottle TOC v. Vial TOC

Because the same samples are measured for TOC and PTOC for T0, these
data are omitted from the regression analysis.

### Absolute

``` r
toc_ptoc.data <- oc_p %>% 
  filter(!Cruise == "AT34", !Timepoint == 0, !Treatment == "Parallel") 
toc_ptoc.reg <- lmodel2(ptoc ~ toc, data = toc_ptoc.data, nperm = 99)
```

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ptoc ~ toc, data = toc_ptoc.data, nperm = 99)
    ## 
    ## n = 72   r = 0.9607479   r-square = 0.9230365 
    ## Parametric P-values:   2-tailed = 1.034424e-40    1-tailed = 5.172118e-41 
    ## Angle between the two OLS regression lines = 2.287741 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  4.590834 0.8938486        41.79189              0.01
    ## 2     MA  2.378699 0.9276331        42.85002              0.01
    ## 3    SMA  2.199658 0.9303675        42.93411                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS      0.5468351        8.634833  0.8323213   0.9553760
    ## 2     MA     -1.9423592        6.432959  0.8657150   0.9936259
    ## 3    SMA     -1.9620843        6.095265  0.8708724   0.9939271
    ## 
    ## Eigenvalues: 60.6155 1.206903 
    ## 
    ## H statistic used for computing C.I. of MA: 0.001177878

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-87-1.png" style="display: block; margin: auto;" />

### Delta

``` r
delta_toc_ptoc.reg <- lmodel2(ptoc_from_t0 ~ toc_from_t0, data = toc_ptoc.data, nperm = 99)
```

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = ptoc_from_t0 ~ toc_from_t0, data =
    ## toc_ptoc.data, nperm = 99)
    ## 
    ## n = 72   r = 0.9158719   r-square = 0.8388214 
    ## Parametric P-values:   2-tailed = 1.864721e-29    1-tailed = 9.323607e-30 
    ## Angle between the two OLS regression lines = 4.975145 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept     Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS  2.529435 0.7907834        38.33636              0.01
    ## 2     MA  2.490784 0.8519461        40.42921              0.01
    ## 3    SMA  2.483532 0.8634214        40.80803                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS       2.201944        2.856927  0.7081515   0.8734153
    ## 2     MA       2.431858        2.544742  0.7665617   0.9451910
    ## 3    SMA       2.428820        2.533258  0.7847345   0.9499984
    ## 
    ## Eigenvalues: 25.98195 1.114466 
    ## 
    ## H statistic used for computing C.I. of MA: 0.002660834

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-90-1.png" style="display: block; margin: auto;" />

## Bottle v. Vial Cell Abundance

``` r
btl_vial_cell.data <- merge %>% 
 drop_na(p_cells)
btl_vial_cell.reg <- lmodel2(p_cells ~ cells, data = btl_vial_cell.data, nperm = 99)
```

    ## 
    ## Model II regression
    ## 
    ## Call: lmodel2(formula = p_cells ~ cells, data = btl_vial_cell.data,
    ## nperm = 99)
    ## 
    ## n = 17   r = 0.979752   r-square = 0.959914 
    ## Parametric P-values:   2-tailed = 6.869573e-12    1-tailed = 3.434786e-12 
    ## Angle between the two OLS regression lines = 1.169037 degrees
    ## 
    ## Permutation tests of OLS, MA, RMA slopes: 1-tailed, tail corresponding to sign
    ## A permutation test of r is equivalent to a permutation test of the OLS slope
    ## P-perm for SMA = NA because the SMA slope cannot be tested
    ## 
    ## Regression results
    ##   Method Intercept    Slope Angle (degrees) P-perm (1-tailed)
    ## 1    OLS -62025907 1.051394        46.43514              0.01
    ## 2     MA -89040705 1.074686        47.06169              0.01
    ## 3    SMA -87227207 1.073123        47.02008                NA
    ## 
    ## Confidence intervals
    ##   Method 2.5%-Intercept 97.5%-Intercept 2.5%-Slope 97.5%-Slope
    ## 1    OLS     -226182286       102130471  0.9331511    1.169637
    ## 2     MA     -238461187        43365926  0.9605251    1.203517
    ## 3    SMA     -231900733        42380931  0.9613744    1.197860
    ## 
    ## Eigenvalues: 1.3175e+18 1.340652e+16 
    ## 
    ## H statistic used for computing C.I. of MA: 0.003145633

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-93-1.png" style="display: block; margin: auto;" />

# Tidy and Wrangle Merged Data

## Condense DOC Data

Because DOC and TOC samples are essentially the same in this dataset, we
will treat them as such. We’ll use the DOC samples from AT34 and the
parallel TOC samples from AT38 and AT39.

``` r
tidy_merge <- merge %>% 
  mutate(ptoc = ifelse(Cruise == "AT34", doc, ptoc),
         sd_ptoc = ifelse(Cruise == "AT34", sd_doc, sd_ptoc),
         ptoc_from_t0 = ifelse(Cruise == "AT34", doc_from_t0, ptoc_from_t0)) %>% 
  select(-c(p_cells, sd_p_cells, doc, sd_doc, doc_from_t0)) %>% 
  rename(doc = ptoc,
         sd_doc = sd_ptoc,
         doc_from_t0 = ptoc_from_t0) %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) %>% 
  mutate(s.timepoint = ifelse(s.timepoint == Timepoint, Hours, NA)) %>% 
  fill(s.timepoint, .direction = "downup") %>% 
  rename(s.gf75 = s.timepoint) %>% 
  select(-Timepoint) %>% 
  # merge duplicate stationary timepoints
  group_by(Cruise, Station, Depth, Treatment, Bottle, Hours) %>%
  fill(cells:ccf.lee, .direction = "downup") %>% 
  distinct() %>% 
  ungroup() %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) 

tidy_merge_keys <- tidy_merge %>% 
  group_keys() %>%
  mutate(key = paste(Cruise, ", S", Station, ", Z =", Depth, ",", Treatment, ",", Bottle))
tidy_merge_header <- tidy_merge_keys$key

tidy_merge_list <- tidy_merge %>%
  group_split()
names(tidy_merge_list) <- tidy_merge_header
```

## Interpolate Stationary Timepoint

``` r
interp.func <- function(x) {
  y <- zoo(x, order.by = x$Hours)
  interp_cells <- round(as.numeric(na.approx(y$cells, na.rm = F)))
  interp_doc <- round(as.numeric(na.approx(y$doc, na.rm = F)), 1)
  z <- cbind(y, interp_cells,  interp_doc)
  as_tibble(z)
}

interp_st <- lapply(tidy_merge_list, interp.func) %>% 
  plyr::ldply(., as.data.frame) %>% 
  select(-.id) %>% 
  mutate_at(vars(Depth, Hours:interp_doc), as.numeric) %>% 
  select(Season:doc_from_t0, interp_cells, interp_doc, everything())
```

# Calculate Derived Variables

We need to readjust our initial bacterial carbon numbers before
calculating other parameters. Since the initial condition was 1.2-µm
filtrate, the carbon content of that filtrate would be higher than the
30:70 incubation mix at the beginning of the experiment. As a result
we’ll correct the initial carbon number to account for this:
POC<sub>initial</sub> = 0.3(POC~1.2 µm filtrate~)

We’ll calculate and define:

  - Division rates (**r**) between each timepoint (Monod 1949). We
    assume all cells are viable.
  - Duration of lag, log, stationary phases based on calculated division
    rates and timing of stationary determined via GrowthCurver
      - lag phase ends when r \> 0.01 (Yates and Smotzer 2007)
      - we’ve already filtered out the death phases so stationary runs
        from GrowthCurver determined timepoint until cell counts cease
      - we do not determine growth phases where there is no observable
        growth
  - Carrying Capacity (**k**), mean of cell abundance during stationary
  - Exponential growth rate day<sup>-1</sup>
  - Dynamic CCFs: we use a weighted function to determine the CCF from
    the initial condition to the beginning of stationary, end members
    which we have measured:
      - CCF<sub>t</sub> = CCF<sub>initial</sub> +
        w(CCF<sub>stationary</sub> - CCF<sub>initial</sub>), where w =
        cells<sub>t</sub> / cells<sub>st</sub>
      - for this to work, we normalize the cell counts to those at
        stationary
  - Cell carbon in µmol C L<sup>-1</sup>, using the dynamic CCFS
  - Bacterial growth efficiencies (BGE) using multiple approaches, only
    where ∆DOC is resolvable:
      - point-to-point (**BGE\_p**), ∆BC and ∆DOC from T0 to stationary
      - phase-to-phase (**BGE\_ph**), ∆BC and ∆DOC using means of lag
        phase and stationary phase values
      - area under the curve (**BGE\_ac**)

<!-- end list -->

``` r
interp_st %>% 
  select(Cruise, Station, Depth, Treatment, Bottle, stationary, s.gf75) %>% 
  unique() %>% 
  filter(s.gf75 < stationary) 
```

    ##   Cruise Station Depth Treatment Bottle stationary s.gf75
    ## 1   AT34       3   200   Control      C        289    288
    ## 2   AT34       3   200   Control      D        555    288
    ## 3   AT38       1    10    SynExd      J        171    120
    ## 4   AT38       1   200   Control      E        390    240
    ## 5   AT38       1   200   Control      F        379    240
    ## 6   AT38       3    10   Control      A        121    120
    ## 7   AT38       6    10   Control      B        196    192
    ## 8   AT39       4    10   Control      A        273    216
    ## 9   AT39       4    10   Control      B        262    216

There are a few cases where the onset of stationary predicted by the
GrowthCurver model occurs later than the time of sampling (n = 9 of 78).

``` r
calcs <- interp_st %>%
  mutate_at(vars( i.poc.cf1:i.pon.c.um), funs(.*0.3)) %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) %>% 
  #division rates 
  mutate(r = round((log(interp_cells) - lag(log(interp_cells)))/(Hours - lag(Hours)), 2),
         r = ifelse(Hours == 0, 0, r)) %>% 
  #cell growth phase
  mutate(cell_div = ifelse(Hours >= stationary & !is.na(interp_cells), "stationary", NA),
         cell_div = ifelse(Hours > stationary & is.na(cell_div), "out of bounds", cell_div),
         cell_div = ifelse(Hours < stationary & r > 0.01, "exponential", cell_div)) %>% 
  fill(cell_div, .direction = "down") %>% 
  mutate(cell_div = ifelse(Hours < stationary & r == 0 & is.na(cell_div), "lag", cell_div)) %>% 
  fill(cell_div, .direction = "down") %>% 
  mutate(full_curve = ifelse(cell_div == "stationary", cell_div, NA)) %>%
  fill(full_curve, .direction = "downup") %>% 
  mutate(full_curve = ifelse(full_curve == "stationary", T, F),
         full_curve = ifelse(is.na(full_curve), F, full_curve),
         cell_div = ifelse(full_curve == F, "no growth", cell_div),
         full_curve = ifelse(cell_div == "exponential", T, NA)) %>% 
  fill(full_curve, .direction = "downup") %>% 
  mutate(full_curve = ifelse(is.na(full_curve), F, full_curve),
         cell_div = ifelse(full_curve == F, "no growth", cell_div)) %>% 
  # carrying capacity
  mutate(k = ifelse(cell_div == "stationary", interp_cells, NA),
         k = ifelse(!is.na(k), mean(k, na.rm = T), NA)) %>% 
  #exponential growth rate
  mutate(end_lag_t = ifelse(cell_div == "lag", Hours, NA),
         end_lag_t = ifelse(!is.na(end_lag_t), max(end_lag_t, na.rm = T), NA),
         end_lag_cells = ifelse(Hours == end_lag_t, interp_cells, NA),
         beg_stat_cells = ifelse(Hours == stationary, interp_cells, NA)) %>% 
  fill(k:beg_stat_cells, .direction = "downup") %>% 
  mutate(mew = ifelse(!is.na(end_lag_cells), ((log(beg_stat_cells) - log(end_lag_cells))/(stationary - end_lag_t)) * 24, NA )) %>% 
  #weighted CCFs, raw POC BC, and BC based on CCFs
  mutate(i.cells = ifelse(Hours == 0, cells, NA)) %>% 
  fill(i.cells, .direction = "updown") %>% 
  mutate(delta_cells = interp_cells - i.cells,
         s.delta_cells = ifelse(Hours == stationary & stationary != 0 & full_curve == T, delta_cells, NA)) %>% 
  fill(s.delta_cells, .direction = "updown") %>% 
  mutate(norm_cells = ifelse(!is.na(s.delta_cells), round(delta_cells/s.delta_cells,2), NA)) %>% 
  rename(weight = norm_cells) %>% 
  mutate(dyn.ccf.c1 = round(i.ccf.c1 + (weight * (s.ccf.c1 - i.ccf.c1)),1),
         dyn.ccf.c3 = round(i.ccf.c3 + (weight * (s.ccf.c3 - i.ccf.c3)),1),
         dyn.ccf.c1 = round(i.ccf.c1 + (weight * (s.ccf.c1 - i.ccf.c1)),1),
         dyn.ccf.white_fukuda = round(ccf.white + (weight * (ccf.fukuda - ccf.white)),1),
         dyn.ccf.white_lee = round(ccf.white + (weight * (ccf.lee - ccf.white)),1),
         bc.ccf.c1 =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * i.ccf.c1) / (12*10^9), 1), NA), 
         bc.ccf.c3 =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * i.ccf.c3) / (12*10^9), 1), NA), 
         bc.dyn.ccf.c1 =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * dyn.ccf.c1) / (12*10^9), 1), NA), 
         bc.dyn.ccf.c3 =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * dyn.ccf.c3) / (12*10^9), 1), NA), 
         bc.ccf.white =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * ccf.white) / (12*10^9), 1), NA), 
         bc.ccf.fukuda =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * ccf.fukuda) / (12*10^9), 1), NA), 
         bc.ccf.lee =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * ccf.lee) / (12*10^9), 1), NA), 
         bc.dyn.ccf.white_fukuda =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * dyn.ccf.white_fukuda) / (12*10^9), 1), NA), 
         bc.dyn.ccf.white_lee =  ifelse(cell_div %in% c("lag", "exponential", "stationary"), round((interp_cells * dyn.ccf.white_lee) / (12*10^9), 1), NA)) %>% 
  ungroup() %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle, cell_div) %>% 
  mutate(ph.bc.ccf.c1 = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.ccf.c1, na.rm = T), 1), NA),
         ph.bc.ccf.c3 = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.ccf.c3, na.rm = T), 1), NA),
         ph.bc.dyn.ccf.c1 = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.dyn.ccf.c1, na.rm = T), 1), NA),
         ph.bc.dyn.ccf.c3 = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.dyn.ccf.c3, na.rm = T), 1), NA),
         ph.bc.ccf.white = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.ccf.white, na.rm = T), 1), NA),
         ph.bc.ccf.fukuda = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean( bc.ccf.fukuda, na.rm = T), 1), NA),
         ph.bc.ccf.lee = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.ccf.lee, na.rm = T), 1), NA),
         ph.bc.dyn.ccf.white_fukuda =  ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.dyn.ccf.white_fukuda, na.rm = T), 1), NA),
         ph.bc.dyn.ccf.white_lee = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(bc.dyn.ccf.white_lee, na.rm = T), 1), NA),
         ph.doc = ifelse(!cell_div %in% c("no growth", "out of bounds"), round(mean(interp_doc, na.rm = T), 1), NA)
         ) %>% 
  ungroup() %>% 
  group_by(Cruise, Station, Depth, Treatment, Bottle) %>% 
  fill(ph.bc.ccf.c1:ph.doc, .direction = "down") %>%
  mutate(del.bc.poc = s.poc.um - i.poc.um,
         del.bc.poc.c1 = s.poc.c1.um - i.poc.c1.um,
         del.bc.poc.c3 = s.poc.c3.um - i.poc.c3.um,
         del.p.bc.ccf.c1 = ifelse(Hours == stationary, bc.ccf.c1, NA) - first(bc.ccf.c1),
         del.p.bc.ccf.c3 = ifelse(Hours == stationary, bc.ccf.c3, NA) - first(bc.ccf.c3),
         del.p.bc.dyn.ccf.c1 = ifelse(Hours == stationary, bc.dyn.ccf.c1, NA) - first(bc.dyn.ccf.c1),
         del.p.bc.dyn.ccf.c3 = ifelse(Hours == stationary, bc.dyn.ccf.c3, NA) - first(bc.dyn.ccf.c3),
          del.p.bc.ccf.white = ifelse(Hours == stationary, bc.ccf.white, NA) - first(bc.ccf.white),
         del.p.bc.ccf.fukuda = ifelse(Hours == stationary, bc.ccf.fukuda, NA) - first(bc.ccf.fukuda),
         del.p.bc.ccf.lee = ifelse(Hours == stationary, bc.ccf.lee, NA) - first(bc.ccf.lee),
         del.p.bc.dyn.ccf.white_fukuda = ifelse(Hours == stationary, bc.dyn.ccf.white_fukuda, NA) - first(bc.dyn.ccf.white_fukuda),
         del.p.bc.dyn.ccf.white_lee = ifelse(Hours == stationary, bc.dyn.ccf.white_lee, NA) - first(bc.dyn.ccf.white_lee),
         
         
          del.ph.bc.ccf.c1 = ifelse(Hours == stationary, ph.bc.ccf.c1, NA) - first(ph.bc.ccf.c1),
          del.ph.bc.ccf.c3 = ifelse(Hours == stationary, ph.bc.ccf.c3, NA) - first(ph.bc.ccf.c3),
         del.ph.bc.dyn.ccf.c1 = ifelse(Hours == stationary, ph.bc.dyn.ccf.c1, NA) - first(ph.bc.dyn.ccf.c1),
         del.ph.bc.dyn.ccf.c3 = ifelse(Hours == stationary, ph.bc.dyn.ccf.c3, NA) - first(ph.bc.dyn.ccf.c3),
          del.ph.bc.ccf.white = ifelse(Hours == stationary, ph.bc.ccf.white, NA) - first(ph.bc.ccf.white),
         del.ph.bc.ccf.fukuda = ifelse(Hours == stationary, ph.bc.ccf.fukuda, NA) - first(ph.bc.ccf.fukuda),
         del.ph.bc.ccf.lee = ifelse(Hours == stationary, ph.bc.ccf.lee, NA) - first(ph.bc.ccf.lee),
         del.ph.bc.dyn.ccf.white_fukuda = ifelse(Hours == stationary, ph.bc.dyn.ccf.white_fukuda, NA) - first(ph.bc.dyn.ccf.white_fukuda),
         del.ph.bc.dyn.ccf.white_lee = ifelse(Hours == stationary, ph.bc.dyn.ccf.white_lee, NA) - first(ph.bc.dyn.ccf.white_lee),
         del.ph.doc = first(ph.doc) - last(ph.doc),
         del.doc = first(doc) - ifelse(Hours == stationary, interp_doc, NA),
         ) %>% 
  fill(del.p.bc.ccf.c1:del.doc, .direction = "downup") %>% 
  mutate(del.ph.doc.flag = ifelse(del.ph.doc >= 1.5, "Acceptable", "NR"),
         del.doc.flag = ifelse(del.doc >= 1.5, "Acceptable", "NR"),
         bge.bc.poc = ifelse(del.doc.flag != "NR", del.bc.poc/(del.doc + del.bc.poc), NA),
         bge.bc.poc.c1 = ifelse(del.doc.flag != "NR", del.bc.poc.c1/(del.doc + del.bc.poc.c1), NA),
         bge.bc.poc.c3 = ifelse(del.doc.flag != "NR", del.bc.poc.c3/(del.doc + del.bc.poc.c3), NA),
         bge.p.bc.ccf.c1 = ifelse(del.doc.flag != "NR", del.p.bc.ccf.c1/(del.doc + del.p.bc.ccf.c1), NA),
         bge.p.bc.ccf.c3 = ifelse(del.doc.flag != "NR", del.p.bc.ccf.c3/(del.doc + del.p.bc.ccf.c3), NA),
         bge.p.bc.dyn.ccf.c1 = ifelse(del.doc.flag != "NR", del.p.bc.dyn.ccf.c1/(del.doc + del.p.bc.dyn.ccf.c1), NA),
         bge.p.bc.dyn.ccf.c3 = ifelse(del.doc.flag != "NR", del.p.bc.dyn.ccf.c3/(del.doc + del.p.bc.dyn.ccf.c3), NA),
         bge.p.bc.ccf.white = ifelse(del.doc.flag != "NR", del.p.bc.ccf.white/(del.doc + del.p.bc.ccf.white), NA),
         bge.p.bc.ccf.fukuda = ifelse(del.doc.flag != "NR", del.p.bc.ccf.fukuda/(del.doc +  del.p.bc.ccf.fukuda), NA),
         bge.p.bc.ccf.lee = ifelse(del.doc.flag != "NR", del.p.bc.ccf.lee/(del.doc + del.p.bc.ccf.lee), NA),
         bge.p.bc.dyn.ccf.white_fukuda = ifelse(del.doc.flag != "NR", del.p.bc.dyn.ccf.white_fukuda/(del.doc + del.p.bc.dyn.ccf.white_fukuda), NA),
         bge.p.bc.dyn.ccf.white_lee = ifelse(del.doc.flag != "NR", del.p.bc.dyn.ccf.white_lee/(del.doc + del.p.bc.dyn.ccf.white_lee), NA),
         
         bge.ph.bc.ccf.c1 = ifelse(del.ph.doc.flag != "NR", del.ph.bc.ccf.c1/(del.ph.doc + del.ph.bc.ccf.c1), NA),
         bge.ph.bc.ccf.c3 = ifelse(del.ph.doc.flag != "NR", del.ph.bc.ccf.c3/(del.ph.doc + del.ph.bc.ccf.c3), NA),
         bge.ph.bc.dyn.ccf.c1 = ifelse(del.ph.doc.flag != "NR", del.ph.bc.dyn.ccf.c1/(del.ph.doc + del.ph.bc.dyn.ccf.c1), NA),
         bge.ph.bc.dyn.ccf.c3 = ifelse(del.ph.doc.flag != "NR", del.ph.bc.dyn.ccf.c3/(del.ph.doc + del.ph.bc.dyn.ccf.c3), NA),
         bge.ph.bc.ccf.white = ifelse(del.ph.doc.flag != "NR", del.ph.bc.ccf.white/(del.ph.doc + del.ph.bc.ccf.white), NA),
         bge.ph.bc.ccf.fukuda = ifelse(del.ph.doc.flag != "NR", del.ph.bc.ccf.fukuda/(del.ph.doc + del.ph.bc.ccf.fukuda), NA),
         bge.ph.bc.ccf.lee = ifelse(del.ph.doc.flag != "NR", del.ph.bc.ccf.lee/(del.ph.doc + del.ph.bc.ccf.lee), NA),
         bge.ph.bc.dyn.ccf.white_fukuda = ifelse(del.ph.doc.flag != "NR", del.ph.bc.dyn.ccf.white_fukuda/(del.ph.doc + del.ph.bc.dyn.ccf.white_fukuda), NA),
         bge.ph.bc.dyn.ccf.white_lee = ifelse(del.ph.doc.flag != "NR", del.ph.bc.dyn.ccf.white_lee/(del.ph.doc + del.ph.bc.dyn.ccf.white_lee), NA)
         ) %>% 
  ungroup() %>% 
  mutate_at(vars(bge.bc.poc:bge.ph.bc.dyn.ccf.white_lee), round, 2) %>% 
  left_join(., readRDS("~/naames_export_ms/Output/processed_export_n2_n3_n4.rds") %>% 
              mutate(Cruise = gsub("AT39-6", "AT39", Cruise)) %>%
              mutate_at(vars(Station), as.character) %>% 
              select(Cruise:Subregion) %>% 
              distinct()
  )
 
# bge <- calcs %>%
#   group_by(Cruise, Station, Depth, Treatment, Bottle) %>%
#   filter(full_curve == T) %>% 
#   drop_na(interp_doc_from_t0, ccf_bc_from_t0) %>% 
#   mutate(auc_bc = round(integrateTrapezoid(Hours, ccf_bc_from_t0, type = "A"), 1),
#          auc_doc = round(integrateTrapezoid(Hours, interp_doc_from_t0, type = "A"), 1),
#          bge_ac = round(auc_bc/auc_doc, 2),
#          bge_ac = ifelse(ddoc_resolve_p == F, NA, bge_ac)) %>% 
#   select(Cruise, Station, Depth, Treatment, Bottle, auc_bc:bge_ac) %>%  
#   distinct() %>% 
#   drop_na() %>% 
#   left_join(calcs, .) %>% 
#   ungroup() 

bge_summary <- calcs %>% 
  select(Season:Station, ave_lat:Subregion, Depth:facet_bottle, stationary:s.gf75, gf75.flag, del.doc, del.doc.flag, del.ph.doc, del.ph.doc.flag, del.bc.poc:del.ph.bc.dyn.ccf.white_lee, bge.bc.poc:bge.ph.bc.dyn.ccf.white_lee) %>% 
  distinct() %>% 
  mutate(type = ifelse(Treatment == "Control", "Control", "Non-Control")) %>% 
  select(Season:Treatment, type, everything())
```

## BGEs

``` r
bge_lvls <- c("p2p, POC", 
              "p2p, POC 1f corr.", 
              "p2p, POC 2f corr.", 
              "p2p, CCF 1f corr.", 
              "ph2ph, CCF 1f corr.", 
              "p2p, CCF 2f corr.", 
              "ph2ph, CCF 2f corr.", 
              "p2p, Dyn. CCF 1f corr.", 
              "ph2ph, Dyn. CCF 1f corr.", 
              "p2p, Dyn. CCF 2f corr.", 
              "ph2ph, Dyn. CCF 2f corr.", 
              "p2p, White CCF", 
              "ph2ph, White CCF", 
              "p2p, Fukuda CCF", 
              "ph2ph, Fukuda CCF", 
              "p2p, Lee CCF", 
              "ph2ph, Lee CCF", 
              "p2p, Dyn. White_Fukuda CCF", 
              "ph2ph, Dyn. White_Fukuda CCF", 
              "p2p, Dyn. White_Lee CCF", 
              "ph2ph, Dyn. White_Lee CCF" )

bge_summary.plot.data <- bge_summary %>%
  select(Season:Bottle, contains("bge")) %>% 
  distinct() %>% 
  gather(key, bge, contains("bge")) %>% 
  drop_na(bge) %>% 
  mutate(appr = key,
         appr = ifelse(appr == "bge.bc.poc", "p2p, POC", appr),
         appr = ifelse(appr == "bge.bc.poc.c1", "p2p, POC 1f corr.", appr),
         appr = ifelse(appr == "bge.bc.poc.c3", "p2p, POC 2f corr.", appr),
         appr = ifelse(appr == "bge.p.bc.ccf.c1", "p2p, CCF 1f corr.", appr),
         appr = ifelse(appr == "bge.p.bc.ccf.c3", "p2p, CCF 2f corr.", appr),
         appr = ifelse(appr == "bge.p.bc.dyn.ccf.c1", "p2p, Dyn. CCF 1f corr.", appr),
         appr = ifelse(appr == "bge.p.bc.dyn.ccf.c3", "p2p, Dyn. CCF 2f corr.", appr),
         appr = ifelse(appr == "bge.p.bc.ccf.white", "p2p, White CCF", appr),
         appr = ifelse(appr == "bge.p.bc.ccf.fukuda", "p2p, Fukuda CCF", appr),
         appr = ifelse(appr == "bge.p.bc.ccf.lee", "p2p, Lee CCF", appr),
         appr = ifelse(appr == "bge.p.bc.dyn.ccf.white_fukuda", "p2p, Dyn. White_Fukuda CCF", appr),
         appr = ifelse(appr == "bge.p.bc.dyn.ccf.white_lee", "p2p, Dyn. White_Lee CCF", appr),
         appr = ifelse(appr == "bge.ph.bc.ccf.c1", "ph2ph, CCF 1f corr.", appr),
         appr = ifelse(appr == "bge.ph.bc.ccf.c3", "ph2ph, CCF 2f corr.", appr),
         appr = ifelse(appr == "bge.ph.bc.dyn.ccf.c1", "ph2ph, Dyn. CCF 1f corr.", appr),
         appr = ifelse(appr == "bge.ph.bc.dyn.ccf.c3", "ph2ph, Dyn. CCF 2f corr.", appr),
         appr = ifelse(appr == "bge.ph.bc.ccf.white", "ph2ph, White CCF", appr),
         appr = ifelse(appr == "bge.ph.bc.ccf.fukuda", "ph2ph, Fukuda CCF", appr),
         appr = ifelse(appr == "bge.ph.bc.ccf.lee", "ph2ph, Lee CCF", appr),
         appr = ifelse(appr == "bge.ph.bc.dyn.ccf.white_fukuda", "ph2ph, Dyn. White_Fukuda CCF", appr),
         appr = ifelse(appr == "bge.ph.bc.dyn.ccf.white_lee", "ph2ph, Dyn. White_Lee CCF", appr),
         p = appr) %>% 
  separate(., p, into = c("points", "ccf"), sep = "," ) %>% 
  group_by(Cruise, type, appr) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!Station == "U")
```

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-99-1.png" style="display: block; margin: auto;" />

Excluding non-control and 200 m
experiments.

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-100-1.png" style="display: block; margin: auto;" />

## Other Parameters

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-101-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-102-1.png" style="display: block; margin: auto;" />

# Merge Data with Export MS Data

Here, we put our experimental data in the context of the results from
the export MS to explore remineralization of the seasonally accumulated
DOC pool:

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
  - Lastly we calculate the rate of DOC drawdown (nmol C d<sup>-1</sup>,
    **ddoc**)

<!-- end list -->

``` r
export <- readRDS("~/naames_export_ms/Output/processed_export_n2_n3_n4.rds") %>% 
  mutate(Cruise = gsub("AT39-6", "AT39", Cruise)) %>% 
  select(Cruise, Season, degree_bin, Station, redis_DOC_vol, NCP_mol_100, doc_ncp_100) %>% 
  distinct() %>% 
  mutate_at(vars(redis_DOC_vol), round, 1) %>% 
  mutate_at(vars(NCP_mol_100:doc_ncp_100), round, 2) %>%
  mutate_at(vars(Station), as.character) 

export_bioav <- left_join(calcs %>% filter(!Station == "U"), export) %>% 
  mutate(redis_DOC_vol = ifelse(Station %in% c("S2RD", "S2RF"), 55.0, redis_DOC_vol),
         degree_bin = ifelse(is.na(degree_bin), 39, degree_bin)) %>%  
  group_by(Cruise, Station, Depth, Treatment, Hours) %>% 
  mutate(trt_ave_doc = ifelse(!is.na(doc), round(mean(doc, na.rm = T),1), NA),
         trt_sd_doc = ifelse(!is.na(doc), round(sd(doc, na.rm = T),1), NA)) %>% 
  ungroup() %>% 
  group_by(Cruise, Station, Depth, Treatment) %>% 
  mutate(norm_doc = ifelse(Depth != 200 | Treatment == "MixDS", trt_ave_doc - redis_DOC_vol, NA),
         accm_doc = ifelse(Depth == 10 & Treatment == "Control", first(trt_ave_doc) - redis_DOC_vol, NA),
         bioav_accm_doc = ifelse(Depth == 10 & Treatment == "Control", last(trt_ave_doc) - redis_DOC_vol, NA),
         norm_doc_from_t0 = ifelse(!is.na(norm_doc), first(norm_doc) - norm_doc, NA),
         bioav_doc = first(norm_doc) - last(norm_doc),
         persis_doc = accm_doc - bioav_doc,
         per_bioav = round((bioav_doc/accm_doc * 100)),
         per_persis = round((persis_doc/accm_doc * 100)),
         time = last(Days),
         ddoc = round((bioav_doc/time) * 1000)) %>%   
  ungroup() %>% 
  group_by(Cruise, Station, Depth) %>% 
  fill(accm_doc:bioav_accm_doc, .direction = "downup") %>% 
  ungroup()
```

# Remineralization of Seasonally Accumulated DOC

### NAAMES 2

#### No Addition

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-104-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-105-1.png" style="display: block; margin: auto;" />

#### Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-106-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-107-1.png" style="display: block; margin: auto;" />

## NAAMES 3

#### No Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-108-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-109-1.png" style="display: block; margin: auto;" />

#### Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-110-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-111-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-112-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-113-1.png" style="display: block; margin: auto;" />

### NAAMES 4

#### No Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-114-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-115-1.png" style="display: block; margin: auto;" />

#### Additions

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-116-1.png" style="display: block; margin: auto;" />

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-117-1.png" style="display: block; margin: auto;" />

## Seasonal Comparison

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-118-1.png" style="display: block; margin: auto;" />

| Season       | min\_bioav\_doc | mean\_bioav\_doc | med\_bioav\_doc | max\_bioav\_doc | min\_persis\_doc | mean\_persis\_doc | med\_persis\_doc | max\_persis\_doc | min\_ddoc | mean\_ddoc | med\_ddoc | max\_ddoc | min\_bioav | mean\_bioav | med\_bioav | max\_bioav | min\_persis | mean\_persis | med\_persis | max\_persis |
| :----------- | --------------: | ---------------: | --------------: | --------------: | ---------------: | ----------------: | ---------------: | ---------------: | --------: | ---------: | --------: | --------: | ---------: | ----------: | ---------: | ---------: | ----------: | -----------: | ----------: | ----------: |
| Early Autumn |             3.8 |              4.9 |             4.9 |             5.9 |              3.9 |               9.1 |              9.3 |             13.0 |        59 |         75 |        76 |        87 |         25 |          37 |         34 |         56 |          44 |           63 |          66 |          75 |
| Early Spring |             1.8 |              3.2 |             3.1 |             4.7 |            \-0.5 |               0.2 |              0.2 |              0.8 |        18 |         36 |        38 |        56 |         82 |          94 |         92 |        112 |        \-12 |            6 |           8 |          18 |
| Late Spring  |             2.1 |              3.0 |             3.0 |             4.0 |            \-0.5 |               4.3 |              5.3 |              7.0 |        38 |        145 |       160 |       224 |         27 |          52 |         34 |        115 |        \-15 |           48 |          66 |          73 |

Seasonal Accumulated DOC Bioavailability and Persistance

This seasonal comparison is best exemplifies at 44˚N where we have
experiments from all cruises. It shows:

  - an increase in the accumulated DOC pool from the early spring to the
    early autumn
  - an increase in the persistent fraction of the bulk DOC as the
    seasons progress
      - DOC in the early spring does not accumulate (also apparent at
        39˚N)
      - 64% of DOC pool in the late spring was persistent while 69% of
        the DOC pool in the early autumn was persistent. This would
        imply an increase of the persistent DOC fraction by 5%. I think
        this suggests that most of the net production of persistent DOC
        occurs in the late spring, closer to the timing of the bloom
        peak.

Overall, there are elevated persistent fractions in the late spring and
early autumn, with greater proportion of the early autumn pool beeing
more consistently persistant.

Our addition experiments are consistent with our definition of DOC
persistence. In every cruise, the addition of SPE-DOM substrates
(different quality and substrates) stimulated drawdown of accumulated
DOC and bacterioplankton growth, but did not result in drawdown into the
persistent
pool.

<img src="NAAMES_DOC_Remin_Bioassays_files/figure-gfm/unnamed-chunk-121-1.png" style="display: block; margin: auto;" />
