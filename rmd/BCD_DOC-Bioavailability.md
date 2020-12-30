BCD and DOC Bioavailability
================
Nicholas Baetge
5/26/2020

# Intro

This document shows plots and tables of the merged field-experiment
NAAMES DOC data.

``` r
library(tidyverse) 
library(patchwork)
library(zoo)
#stat tests
library(lmtest)
library(lmodel2)
library(rstatix)
library(ggpubr)
library(blandr)
```

# Import Data

``` r
export <- readRDS("~/GITHUB/naames_bioav_ms/Input/master/processed_export_for_bioavMS.6.7.20.rds") %>% 
  select(Season, Cruise, degree_bin, Station, int_delta_DOC_100) %>% 
  distinct()

export_bioav <- readRDS("~/GITHUB/naames_bioav_ms/Output/processed_DOC_bioavailability.rds") 

doc_og <- read_rds("~/GITHUB/naames_bioav_ms/Input/master/DOC_Input") %>% 
  filter(Cruise == "AT34", Treatment == "Control", Depth == 10) %>% 
  drop_na(doc) %>% 
  left_join(., export_bioav %>% select(Season, Cruise, Station, Bottle, stationary.harvest)) %>% 
  select(Season, Cruise, Station, Bottle, stationary.harvest, Hours, doc) %>% 
  mutate(Days = Hours/24) %>% 
  filter(Station == 4) %>% 
  group_by(Season, Cruise, Station, Days, stationary.harvest) %>% 
  summarise_at(vars(doc), list(mean = mean, sd = sd)) %>%
  ungroup() %>% 
  mutate(increase = 56.2 - 53.4,
         per_increase = increase / 53.4 * 100)
```

    ## Joining, by = c("Season", "Cruise", "Station", "Bottle")

``` r
bge <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_bge.rds") %>% 
  group_by(Cruise, Station, Treatment, Bottle) %>% 
  mutate(remin.end = last(Days),
         stationary.harvest = stationary.harvest/24) %>%
  ungroup()

bge_summary <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_bge_summary.rds") %>% 
  ungroup()

bcd <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_BCD.rds") 

lat44 <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_lat44_remins.rds") %>% 
  drop_na(sd_combined_doc) %>% 
  mutate(station_line = ifelse(Cruise == "AT39" & Station == 4, "ab", "a"))
```

Units for imported data frames are currently:

  - BP, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - BCD, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - BA, cells m<sup>-3</sup>
  - BC, µmol C m<sup>-3</sup>
  - mew, d<sup>-1</sup>
  - NPP, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - ∆DOC and NCP (from export MS, integrated to Ez), mol C
    m<sup>-2</sup>

NPP, NCP and BCD are converted to: mmol C m<sup>-3</sup> d<sup>-1</sup>

## Bottle v. Vial Incubation Comparisons

We will run the Bland-Altman method (aka Tukey mean-difference) for
assessing agreement between two methods (in this case, vial and bottle
measurements/doc star and filtered doc). The Bland-Altman plot is a
scatter plot in which the difference between the paired measurements
(A-B) is plotted against their mean value \[(A+B)/2\]. The graph
provides two main pieces of information, namely the average of all the
differences, which is also provided by the t test and which is often
called the bias, and the 95% limits of agreement.

We’ll also do a lmodel2 regression.

### DOC filtered v . DOC\*

``` r
vessel_comparisons <- bge %>% 
  select(Season, Cruise, Station, Bottle, Hours, Days, stationary.harvest, cells, p_cells, ptoc,  toc, doc, pdoc, i.ccf, s.ccf, boc_i.ccf, boc_s.ccf) %>% 
  group_by(Cruise, Station, Bottle) %>% 
  mutate(doc_star = ifelse(Hours < stationary.harvest, ptoc - boc_i.ccf, ptoc - boc_s.ccf),
         doc_star = ifelse(is.na(doc_star) & Hours < stationary.harvest, toc - boc_i.ccf, toc - boc_s.ccf),
         combined_doc = doc,
         combined_doc = ifelse(is.na(combined_doc), pdoc, combined_doc))  


doc_doc.star.reg <- lmodel2(doc_star ~ doc, data = vessel_comparisons, nperm = 99) #filtered DOC (from bottles) v DOC star
```

    ## RMA was not requested: it will not be computed.

``` r
doc_doc.star.blandr <- blandr.statistics(vessel_comparisons$combined_doc, vessel_comparisons$doc_star, sig.level = 0.95)

doc_doc.star.blandr$bias
```

    ## [1] 1.606

``` r
doc_doc.star.blandr$upperLOA
```

    ## [1] 4.123638

``` r
doc_doc.star.blandr$lowerLOA
```

    ## [1] -0.9116376

``` r
doc_doc.star.reg.plot + doc_doc.star.blandr.plot +
  plot_annotation(tag_levels = "a") &
  plot_layout(guides = "collect") +
  theme(plot.tag = element_text(size = 14))
```

![](BCD_DOC-Bioavailability_files/figure-gfm/combine%20doc%20star%20plots-1.png)<!-- -->

### PDOC v DOC bottle

``` r
vessel_doc.reg <- lmodel2(pdoc ~ doc, data = vessel_comparisons, nperm = 99) #vial v bottle DOC (both filtered)
```

    ## RMA was not requested: it will not be computed.

``` r
vessel_doc.blandr <- blandr.statistics(vessel_comparisons$doc, vessel_comparisons$pdoc, sig.level = 0.95)

vessel_doc.blandr$bias
```

    ## [1] 0.3833333

``` r
vessel_doc.blandr$upperLOA
```

    ## [1] 2.578183

``` r
vessel_doc.blandr$lowerLOA
```

    ## [1] -1.811517

``` r
vl_btl_doc_blandr <- vessel_comparisons %>% 
  drop_na(pdoc) %>% 
  mutate(diff = doc - pdoc,
         mean = (doc + pdoc) / 2) %>% 
  ggplot(., aes(x = mean, y = diff, group = Bottle)) +
  geom_hline(yintercept = 0, size = 1, linetype = 1) +
    geom_hline(yintercept = vessel_doc.blandr$bias, size = 1, linetype = 2) +
  geom_hline(yintercept = vessel_doc.blandr$upperLOA, size = 1, linetype = 3) +
  geom_hline(yintercept = vessel_doc.blandr$lowerLOA, size = 1, linetype = 3) +
  geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  labs(y = expression(italic(paste("Difference: Bottle DOC - Vial DOC, µmol C L"^-1))), x = expression(italic(paste("Mean: (Bottle DOC + Vial DOC) / 2, µmol C L"^-1)))) +
  custom_theme() +
  guides(fill = F)
```

``` r
vl_btl_doc <- vessel_comparisons %>% 
  drop_na(pdoc) %>% 
  ggplot(aes(x = doc, y = pdoc, group = Bottle)) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_abline(intercept = vessel_doc.reg$regression.results[3,2],
              slope = vessel_doc.reg$regression.results[3,3],colour = "black", linetype = 2, size = 1) +
  geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  labs(x = expression(italic(paste("Bottle DOC, µmol C L"^-1))), y = expression(italic(paste("Vial DOC, µmol C L"^-1))), fill = "Early Spring Station") +
  custom_theme() +
   annotate( geom = "text", label = expression(atop("y = 1.02x - 1.44", paste("r"^2,"= 0.96, ", italic("p "), "<< 0.01"))), x = 80, y = 51, size = 3) +
  xlim(50, 85) +
  ylim(50, 85) +
  guides(fill = F)
```

### Bottle v. Vial Cell Abundance

``` r
vessel_cells.reg <- lmodel2(p_cells ~ cells, data = vessel_comparisons, nperm = 99) #vial v bottle DOC (both filtered)
```

    ## RMA was not requested: it will not be computed.

``` r
vessel_cells.blandr <- blandr.statistics(vessel_comparisons$cells, vessel_comparisons$p_cells, sig.level = 0.95)

vessel_cells.blandr$bias
```

    ## [1] 2417934

``` r
vessel_cells.blandr$upperLOA
```

    ## [1] 342965561

``` r
vessel_cells.blandr$lowerLOA
```

    ## [1] -338129693

``` r
vl_btl_cell_blandr <- vessel_comparisons %>% 
  drop_na(p_cells) %>% 
  mutate(diff = cells - p_cells,
         mean = (cells + p_cells) / 2) %>% 
  ggplot(., aes(x = mean, y = diff, group = Bottle)) +
  geom_hline(yintercept = 0, size = 1, linetype = 1) +
    geom_hline(yintercept = vessel_cells.blandr$bias, size = 1, linetype = 2) +
  geom_hline(yintercept = vessel_cells.blandr$upperLOA, size = 1, linetype = 3) +
  geom_hline(yintercept = vessel_cells.blandr$lowerLOA, size = 1, linetype = 3) +
  geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  labs(y = expression(italic(paste("Difference: Bottle Cells - Vial Cells, L"^-1))), x = expression(italic(paste("Mean: (Bottle Cells + Vial Cells) / 2, L"^-1)))) +
  custom_theme() +
  guides(fill = F)
```

``` r
vl_btl_cell <- vessel_comparisons %>% 
  drop_na(p_cells) %>% 
  ggplot(aes(x = cells, y = p_cells)) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_abline(intercept = vessel_cells.reg$regression.results[3,2],
              slope = vessel_cells.reg$regression.results[3,3],colour = "black", linetype = 2, size = 1) +
 geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  labs(x = expression(italic(paste("Bottle Cells L"^-1))), y = expression(italic(paste("Vial Cells L"^-1))), fill = "NAAMES 4 Station") +
  theme_test(base_size = 16) +
  theme(legend.position = "top") +
  annotate( geom = "text", label = expression(atop("y = 1.07x - 8.72 * 10"^7, paste("r"^2,"= 0.96, ", italic("p "), "<< 0.01"))), x = 2.6E9, y = 5.0E8, size = 3) +
  #theme(plot.caption = element_text(face = "italic")) +
  ylim(4.5*10^8, 3.0*10^9) +
  xlim(4.5*10^8, 3.0*10^9)
```

``` r
vl_btl_cell + vl_btl_cell_blandr + vl_btl_doc + vl_btl_doc_blandr +
  plot_annotation(tag_levels = "a") &
  plot_layout(guides = "collect") +
  theme(plot.tag = element_text(size = 14))
```

![](BCD_DOC-Bioavailability_files/figure-gfm/combine%20regression%20plots-1.png)<!-- -->

# N2 Post-Cruise contamination

``` r
unique(doc_og$increase)
```

    ## [1] 2.8

``` r
unique(doc_og$per_increase)
```

    ## [1] 5.243446

![](BCD_DOC-Bioavailability_files/figure-gfm/N2%20remin%20plot-1.png)<!-- -->

# Experiments at 44˚N

![](BCD_DOC-Bioavailability_files/figure-gfm/Remin%20plots-1.png)<!-- -->

# Tables: Remin Experiments

Seasonally accumulated DOC is expressed as mmol C m<sup>-3</sup> or µmol
C L<sup>-1</sup>

``` r
bge.table.data <- bge %>%
  left_join(., export_bioav %>% select(Season, Cruise, Station, degree_bin) %>% distinct()) %>% 
  mutate(degree_bin = ifelse(is.na(degree_bin), 39, degree_bin)) %>% 
  filter(Treatment == "Control") %>% 
  select(Season, Cruise, Station, degree_bin, Bottle, stationary.harvest, remin.end) %>% 
  distinct() %>% 
  left_join(bge_summary, .) %>% 
  select(Season, Cruise, Station, degree_bin, Bottle, stationary.harvest, remin.end, everything()) %>% 
  select(-n) %>% 
  arrange(factor(Cruise, levels = c("AT39", "AT34", "AT38")), degree_bin) %>% 
  group_by(Season, Cruise, Station, degree_bin) %>% 
  add_tally() %>% 
  ungroup()
```

    ## Joining, by = c("Season", "Cruise", "Station")

    ## Joining, by = c("Season", "Station", "Bottle")

``` r
subset <- bge.table.data %>% 
  select(Season:degree_bin, ave_del.poc_station:sd_del.doc.star_station, ave_bge_station, sd_bge_station) %>% distinct() %>% 
  group_by(Season, Cruise, Station, degree_bin) %>% 
  mutate(combined.deldoc = ifelse(!is.na(ave_del.doc_station) & !is.na(ave_del.doc.star_station), (ave_del.doc_station + ave_del.doc.star_station) / 2, NA),
         sd_combined.deldoc = ifelse(!is.na(ave_del.doc_station) & !is.na(ave_del.doc.star_station), sqrt(mean(c((ave_del.doc_station - combined.deldoc)^2, (ave_del.doc.star_station - combined.deldoc)^2))), NA))

bge.table.season <- subset %>% 
  mutate(combined.deldoc = ifelse(is.na(combined.deldoc), ave_del.doc_station, combined.deldoc),
         combined.deldoc = ifelse(is.na(combined.deldoc), ave_del.doc.star_station, combined.deldoc)) %>% 
  group_by(Season) %>% 
  summarise_at(vars(combined.deldoc, contains("ave")), list(mean = mean, sd = sd), na.rm = T) %>% 
  arrange(factor(Season, levels = levels))


max(bge.table.data$remin.end)
```

    ## [1] 110

``` r
min(bge.table.data$remin.end)
```

    ## [1] 8

``` r
bioav.table.data <- export_bioav %>% 
  select(Cruise:Station, degree_bin, end.remin, stationary.harvest, initial.doc:longterm.ddoc) %>% 
  distinct() %>% 
  group_by(Season, Cruise, Station, degree_bin)

bioav.table <- bioav.table.data %>% 
  summarise_at(vars(initial.doc:longterm.ddoc), list(mean = mean, sd = sd), na.rm = T) %>% 
  arrange(factor(Season, levels = levels), degree_bin) %>% 
  select(Season:degree_bin, contains(c("accm.doc", "shortterm.bioav.doc", "shortterm.ddoc", "longterm.del.doc",  "longterm.bioav.doc", "longterm.ddoc",  "persis") ))

bioav.table.summary <- bioav.table %>% 
  ungroup() %>% 
  group_by(Season) %>% 
  summarise_at(vars(contains("mean")), list(season.ave = mean, sd.season = sd), na.rm = T) %>% 
  arrange(factor(Season, levels = levels))
```

# Box plots: NPP, BP, BA, µ

<img src="BCD_DOC-Bioavailability_files/figure-gfm/patchwork rate & abundance plots-1.png" style="display: block; margin: auto;" />

# Bar plots: BCD and BCD:NPP

Convert BCD and NPP to mmol C m<sup>-3</sup> d<sup>-1</sup>

``` r
bcd.data <- bcd %>% 
  select(Season, Cruise:degree_bin, int.bcd_station, int.bp, int.NPP, bp.npp, bcd.npp_station) %>% 
  rename(int.bcd = int.bcd_station,
         bcd.npp = bcd.npp_station) %>% 
  distinct() %>% 
  mutate_at(vars(int.bcd, int.bp, int.NPP), funs(./10^3)) %>% 
  mutate_at(vars(bp.npp, bcd.npp), funs(./10^2)) %>% 
  mutate(degree_bin = as.character(degree_bin)) 

bcd.summary <- bcd.data %>% 
  group_by(Season, Cruise, degree_bin) %>% 
  summarise_at(vars(int.bcd, int.bp, int.NPP, bp.npp, bcd.npp), list(mean = mean, sd = sd)) %>% 
  arrange(factor(Season, levels = levels), degree_bin) %>% 
  select(Season:degree_bin, contains(c("int.NPP", "int.bp", "bp.npp", "int.bcd", "bcd.npp")))
  
bcd.cruise.summary <- bcd.data %>% 
  group_by(Season, Cruise) %>% 
  summarise_at(vars(int.bcd, int.bp, int.NPP, bp.npp, bcd.npp), list(cruise_mean = mean, cruise_sd = sd, cruise_max = max, cruise_min = min)) %>% 
  arrange(factor(Season, levels = levels)) %>% 
  select(Season, contains(c("int.NPP", "int.bp", "bp.npp", "int.bcd", "bcd.npp")))

bcd.overall.summary <- bcd.data %>% 
  summarise_at(vars(int.bcd, int.bp, int.NPP, bp.npp, bcd.npp), list(overall_mean = mean, overall_sd = sd, max = max, min = min))

bcd.npp.summary <- bcd.summary %>% 
  select(Season, Cruise, degree_bin, int.NPP_mean, int.bcd_mean) %>% 
  rename(NPP = int.NPP_mean,
         BCD = int.bcd_mean) %>% 
  pivot_longer(c(NPP, BCD), names_to = "rate", values_to = "ave") %>% 
  left_join(., bcd.summary %>% 
  select(Season, Cruise, degree_bin, int.NPP_sd, int.bcd_sd) %>% 
    rename(NPP = int.NPP_sd,
         BCD = int.bcd_sd) %>% 
  pivot_longer(c(NPP, BCD), names_to = "rate", values_to = "sd"))
```

    ## Joining, by = c("Season", "Cruise", "degree_bin", "rate")

``` r
bcd.npp_summary_table <- left_join(bcd.summary, bcd.cruise.summary)
```

    ## Joining, by = "Season"

``` r
field.doc_table <- export_bioav %>% 
  select(Cruise, Season, accm.doc) %>% 
  drop_na(accm.doc) %>% 
  group_by(Cruise, Season) %>% 
  summarise_at(vars(accm.doc), list(mean = mean, sd = sd))

field.doc_table
```

    ## # A tibble: 3 x 4
    ## # Groups:   Cruise [3]
    ##   Cruise Season        mean    sd
    ##   <chr>  <chr>        <dbl> <dbl>
    ## 1 AT34   Late Spring   6.51 3.42 
    ## 2 AT38   Early Autumn 12.7  3.33 
    ## 3 AT39   Early Spring  3.58 0.755

![](BCD_DOC-Bioavailability_files/figure-gfm/BCD%20and%20NPP%20bar%20plots-1.png)<!-- -->

# Regressions: Property-Property

## BCD v NPP

    ## RMA was not requested: it will not be computed.

## BCD:NPP v NPP

### exponential model using lm function

``` r
reg2 <- lm(log(bcd.npp) ~ int.NPP, data = bcd.data)

reg2.df <- data.frame(x = bcd.data$int.NPP, y = exp(fitted(reg2)))
```

### manual exponential model

y = alpha \* (e^beta \* x) + theta

    ## $alpha
    ## (Intercept) 
    ##     0.22086 
    ## 
    ## $beta
    ##    int.NPP 
    ## -0.3066193 
    ## 
    ## $theta
    ## [1] 0.05297107

``` r
model <- nls(bcd.npp ~ alpha * exp(beta * int.NPP) + theta , data = bcd.data, start = start)

model.df <- data.frame(x = bcd.data$int.NPP, y = predict(model, list(x = bcd.data$int.NPP)) )

summary(model)
```

    ## 
    ## Formula: bcd.npp ~ alpha * exp(beta * int.NPP) + theta
    ## 
    ## Parameters:
    ##        Estimate Std. Error t value Pr(>|t|)    
    ## alpha   2.15650    1.63651   1.318  0.19725    
    ## beta  -15.81365    5.73445  -2.758  0.00967 ** 
    ## theta   0.21637    0.02043  10.591    8e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09965 on 31 degrees of freedom
    ## 
    ## Number of iterations to convergence: 11 
    ## Achieved convergence tolerance: 6.218e-07

## Plots

![](BCD_DOC-Bioavailability_files/figure-gfm/BCD%20and%20NPP%20regression%20plots-1.png)<!-- -->
