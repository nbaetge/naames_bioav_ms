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
  select(Season, Cruise, degree_bin, Station, int_delta_DOC_100, NCP_mol_100, doc_ncp_100) %>% 
  distinct() %>% 
  mutate(Season = ifelse(Season == "Late Autumn", "Early Winter", Season))

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
  mutate(doc_star = ifelse(Cruise != "AT34" & Hours < stationary.harvest, round(interp_ptoc - i.poc,1), round(interp_ptoc - s.poc,1)),
         combined_doc = ifelse(Cruise == "AT34", interp_doc, doc_star),
          sd_combined_doc = ifelse(Cruise == "AT34", sd_doc, sd_ptoc)) %>% 
  group_by(Cruise, Station, Treatment, Bottle) %>% 
  mutate(remin.end = last(Days),
         stationary.harvest = stationary.harvest/24) %>%
  ungroup() %>% 
  group_by(Cruise, Station, Treatment, Hours) %>% 
  mutate(trt_combined_doc = ifelse(!is.na(sd_combined_doc), round(mean(combined_doc, na.rm = T), 1), NA),
         sd_trt_combined_doc = ifelse(!is.na(sd_combined_doc), round(sd(combined_doc, na.rm = T), 1), NA)) %>% 
  ungroup() %>% 
  group_by(Season, Station, Treatment) %>% 
  mutate(norm_doc =  round(trt_combined_doc - first(trt_combined_doc), 1),
         sd_norm_doc = ifelse(Hours != 0,  sd_trt_combined_doc, NA),
         norm_cells = station_cells - first(station_cells),
         sd_norm_cells = ifelse(Hours != 0, sd_station_cells, NA)) %>% 
  ungroup() %>% 
  select(Season:stationary.harvest, remin.end, Hours:sd_station_cells, norm_cells, sd_norm_cells, everything()) %>% 
  left_join(., export_bioav %>% select(Season, Cruise, Station, degree_bin) %>% distinct()) %>% 
  mutate(degree_bin = ifelse(is.na(degree_bin), 39, degree_bin)) %>% 
  select(Season:Station, degree_bin, everything()) %>% 
  mutate(Season = ifelse(Season == "Late Autumn", "Early Winter", Season))
```

    ## Joining, by = c("Season", "Cruise", "Station")

``` r
bge_summary <- read_rds("~/GITHUB/naames_bioav_ms/Output/bge_summary.rds") %>% 
  left_join(., read_rds("~/GITHUB/naames_bioav_ms/Output/bge_cruise_means.rds") %>% select(Season, Cruise, bge_mean) %>% rename(cruise_bge = bge_mean)) %>% 
  left_join(., read_rds("~/GITHUB/naames_bioav_ms/Output/bge_station_means.rds") %>% select(Season, Cruise, Station, bge_mean) %>% rename(station_bge = bge_mean)) %>% 
  left_join(., export_bioav %>% select(Season, Cruise, Station, degree_bin) %>% distinct()) %>% 
  mutate(degree_bin = ifelse(is.na(degree_bin), 39, degree_bin)) %>% 
  select(Season:Station, degree_bin, everything()) %>% 
  mutate(Season = ifelse(Season == "Late Autumn", "Early Winter", Season)) %>% 
  mutate(ccf = (del.poc * (12 * 10^9)/del.cells))
```

    ## Joining, by = c("Season", "Cruise")

    ## Joining, by = c("Season", "Cruise", "Station")
    ## Joining, by = c("Season", "Cruise", "Station")

``` r
bcd <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_BCD.rds") %>% 
  mutate(Season = ifelse(Season == "Late Autumn", "Early Winter", Season))

lat44 <- read_rds("~/GITHUB/naames_bioav_ms/Output/processed_lat44_remins.rds") %>% 
  drop_na(sd_combined_doc) %>% 
  mutate(station_line = ifelse(Cruise == "AT39" & Station == 4, "ab", "a")) %>% 
  mutate(Season = ifelse(Season == "Late Autumn", "Early Winter", Season))
```

Units for imported data frames are currently:

  - BP, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - BCD, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - BA, cells m<sup>-3</sup>
  - NPP, µmol C m<sup>-3</sup> d<sup>-1</sup>
  - ∆DOC and NCP (from export MS, integrated to Ez), mol C
    m<sup>-2</sup>

NPP, NCP and BCD are converted to: mmol C m<sup>-3</sup> d<sup>-1</sup>

# Bottle v. Vial Incubation Comparisons

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
  select(Season, Cruise, Station, degree_bin, Bottle, Hours, Days,  cells, p_cells, ptoc,  toc, doc, pdoc, doc_star) 

doc_doc.star <- vessel_comparisons %>% 
  drop_na(doc) %>% 
  drop_na(doc_star)

doc_doc.star.reg <- lmodel2(doc_star ~ doc, data = doc_doc.star, nperm = 99) #filtered DOC (from bottles) v DOC star
```

    ## RMA was not requested: it will not be computed.

``` r
doc_doc.star.blandr <- blandr.statistics(doc_doc.star$doc, doc_doc.star$doc_star, sig.level = 0.95)

doc_doc.star.blandr$bias
```

    ## [1] 2.128571

``` r
doc_doc.star.blandr$upperLOA
```

    ## [1] 6.747419

``` r
doc_doc.star.blandr$lowerLOA
```

    ## [1] -2.490276

``` r
doc_doc.star.reg.plot + doc_doc.star.blandr.plot +
  plot_annotation(tag_levels = "a") &
  plot_layout(guides = "collect") +
  theme(plot.tag = element_text(size = 14))
```

![](BCD_DOC-Bioavailability_files/figure-gfm/combine%20doc%20star%20plots-1.png)<!-- -->

### PDOC v DOC bottle

``` r
vessel_doc <- vessel_comparisons %>% 
  drop_na(pdoc)

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
vl_btl_doc_blandr <- vessel_doc %>% 
  mutate(diff = doc - pdoc,
         mean = (doc + pdoc) / 2) %>% 
  ggplot(., aes(x = mean, y = diff, group = Bottle)) +
  geom_hline(yintercept = 0, size = 1, linetype = 1) +
    geom_hline(yintercept = vessel_doc.blandr$bias, size = 1, linetype = 2) +
  geom_hline(yintercept = vessel_doc.blandr$upperLOA, size = 1, linetype = 3) +
  geom_hline(yintercept = vessel_doc.blandr$lowerLOA, size = 1, linetype = 3) +
  geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = expression(paste("Difference: Bottle DOC - Vial DOC, µmol C L"^-1)), x = expression(paste("Mean: (Bottle DOC + Vial DOC) / 2, µmol C L"^-1))) +
  theme_classic2(base_size = 16) +
  theme(legend.title = element_blank()) +
  guides(fill = F)
```

``` r
vl_btl_doc <- vessel_doc %>% 
  ggplot(aes(x = doc, y = pdoc, group = Bottle)) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_abline(intercept = vessel_doc.reg$regression.results[3,2],
              slope = vessel_doc.reg$regression.results[3,3],colour = "black", linetype = 2, size = 1) +
  geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = expression(paste("Bottle DOC, µmol C L"^-1)), y = expression(paste("Vial DOC, µmol C L"^-1)), fill = "Early Spring Station") +
  theme_classic2(base_size = 16) +
  theme(legend.title = element_blank()) +
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
  scale_fill_brewer(palette = "Dark2") +
  labs(y = expression(paste("Difference: Bottle Cells - Vial Cells, L"^-1)), x = expression(paste("Mean: (Bottle Cells + Vial Cells) / 2, L"^-1))) +
  theme_classic2(base_size = 16) 
```

``` r
vl_btl_cell <- vessel_comparisons %>% 
  drop_na(p_cells) %>% 
  ggplot(aes(x = cells, y = p_cells)) + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_abline(intercept = vessel_cells.reg$regression.results[3,2],
              slope = vessel_cells.reg$regression.results[3,3],colour = "black", linetype = 2, size = 1) +
 geom_point(aes(fill = Station), color = "black", shape = 21, size = 4, alpha = 0.7 ) +
  labs(x = expression(paste("Bottle Cells L"^-1)), y = expression(paste("Vial Cells L"^-1)), fill = "NAAMES 4 Station") +
  guides(fill = F) +
  theme_classic2(base_size = 16) +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.title = element_blank()) +
  annotate( geom = "text", label = expression(atop("y = 1.07x - 8.72 * 10"^7, paste("r"^2,"= 0.96, ", italic("p "), "<< 0.01"))), x = 2.6E9, y = 5.0E8, size = 3) +
  #theme(plot.caption = element_text(face = "")) +
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

# Curves

## Cells

``` r
ba_curves <- bge %>% 
  filter(Treatment == "Control") %>% 
  drop_na(norm_cells) %>% 
  select(Season, Station, degree_bin, Days, norm_cells, sd_norm_cells) %>% 
  distinct() %>% 
  # filter(!Cruise == "AT34" | !Station == 3) %>% #did not calc bge or cell resp due to data coming from post stationary
  ggplot(aes(x = Days, y = norm_cells, group = interaction(Season, Station))) +
  geom_errorbar(aes(ymin = norm_cells - sd_norm_cells, ymax = norm_cells + sd_norm_cells, color = factor(degree_bin)), width = 0.3, alpha = 0.5) +
  geom_line(aes(color = factor(degree_bin)), alpha = 0.7) +
  geom_point(aes(fill = factor(degree_bin)), size = 3, shape = 21, alpha = 0.7) +
  scale_color_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu") +
 labs(x = expression(""), y = expression(paste("Bacterial Abundance, Cells L"^-1)), fill = "Latitude, ˚N") +
  theme_classic2(base_size = 16) +
  guides(color = F, linetype = F) + facet_grid(~factor(Season, levels = levels))
```

## DOC

``` r
doc_curves <-  bge %>% 
  filter(Treatment == "Control") %>% 
  drop_na(sd_combined_doc) %>% 
  select(Season, Station, degree_bin, Days, norm_doc, sd_norm_doc) %>% 
  distinct() %>% 
  # filter(!Cruise == "AT34" | !S,tation == 3) %>% #did not calc bge or cell resp due to data coming from post stationary
  ggplot(aes(x = Days, y = norm_doc, group = interaction(Season, Station))) +
  geom_errorbar(aes(ymin = norm_doc - sd_norm_doc, ymax = norm_doc + sd_norm_doc, color = factor(degree_bin)), width = 1, alpha = 0.5) +
  geom_line(aes(color = factor(degree_bin)), alpha = 0.7) +
  geom_point(aes(fill = factor(degree_bin)), size = 3, shape = 21, alpha = 0.7) +
  scale_color_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu") +
 # labs(x = expression("Days"), y = expression(paste("DOC"^"(*)",", µmol C L"^-1)), fill = "") +
  labs(x = expression("Days"), y = expression(paste("DOC(*), µmol C L"^-1)), fill = "") +
  theme_classic2(base_size = 16) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  guides(color = F, fill = F) +
  facet_grid(~factor(Season, levels = levels))
```

``` r
doc_curves2 <-  bge %>% 
  filter(Treatment == "Control", Days <= 20) %>% 
  drop_na(sd_combined_doc) %>% 
  select(Season, Station, degree_bin, Days, norm_doc, sd_norm_doc) %>% 
  distinct() %>% 
  # filter(!Cruise == "AT34" | !S,tation == 3) %>% #did not calc bge or cell resp due to data coming from post stationary
  ggplot(aes(x = Days, y = norm_doc, group = interaction(Season, Station))) +
  geom_errorbar(aes(ymin = norm_doc - sd_norm_doc, ymax = norm_doc + sd_norm_doc, color = factor(degree_bin)), width = 0.5, alpha = 0.3) +
  geom_line(aes(color = factor(degree_bin)), alpha = 0.7) +
  geom_point(aes(fill = factor(degree_bin)), size = 3, shape = 21, alpha = 0.7) +
  scale_color_brewer(palette = "RdBu") +
  scale_fill_brewer(palette = "RdBu") +
 labs(x = expression("Days"), y = expression(paste("DOC, µmol C L"^-1)), fill = "") +
  theme_classic2(base_size = 16) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  guides(color = F, fill = F) +
  facet_grid(~factor(Season, levels = levels))
```

``` r
ba_curves / doc_curves + 
   plot_annotation(tag_levels = "a") &
  plot_layout(guides = "collect") +
  theme(plot.tag = element_text(size = 14))
```

![](BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ba_curves / doc_curves2 + 
   plot_annotation(tag_levels = "a") &
  plot_layout(guides = "collect") +
  theme(plot.tag = element_text(size = 14))
```

![](BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Delta Stats

stats for ∆POC and ∆DOC

``` r
delta_stats <- bge_summary %>% 
  select(Cruise, degree_bin,  del.poc, del.doc, del.doc.star) %>% 
  mutate(del.doc = ifelse(is.na(del.doc), del.doc.star, del.doc)) %>% 
  select(-del.doc.star)

longterm_delta_stats <- export_bioav %>% 
  select(Cruise, degree_bin, longterm.del.doc) %>% 
  drop_na(longterm.del.doc) %>% 
  mutate(season = ifelse(Cruise == "AT38", "Autumn", "Spring"))
```

``` r
# delta_stats %>% 
  # kruskal_test(del.poc ~ Cruise)

compare_means(del.poc ~ Cruise, delta_stats, method = "kruskal.test")
```

    ## # A tibble: 1 x 6
    ##   .y.         p p.adj p.format p.signif method        
    ##   <chr>   <dbl> <dbl> <chr>    <chr>    <chr>         
    ## 1 del.poc 0.443  0.44 0.44     ns       Kruskal-Wallis

``` r
# delta_stats %>%
#   kruskal_test(del.doc ~ Cruise)

compare_means(del.doc ~ Cruise, delta_stats, method = "kruskal.test")
```

    ## # A tibble: 1 x 6
    ##   .y.         p p.adj p.format p.signif method        
    ##   <chr>   <dbl> <dbl> <chr>    <chr>    <chr>         
    ## 1 del.doc 0.132  0.13 0.13     ns       Kruskal-Wallis

``` r
delta_stats %>%
  group_by(Cruise) %>% 
  kruskal_test(del.doc ~ degree_bin)
```

    ## # A tibble: 3 x 7
    ##   Cruise .y.         n statistic    df     p method        
    ## * <chr>  <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 AT34   del.doc    10      1.5      3 0.682 Kruskal-Wallis
    ## 2 AT38   del.doc    12      3.32     2 0.191 Kruskal-Wallis
    ## 3 AT39   del.doc    12      1.5      1 0.221 Kruskal-Wallis

``` r
compare_means(longterm.del.doc ~ Cruise, longterm_delta_stats, method = "kruskal.test")
```

    ## # A tibble: 1 x 6
    ##   .y.                    p  p.adj p.format p.signif method        
    ##   <chr>              <dbl>  <dbl> <chr>    <chr>    <chr>         
    ## 1 longterm.del.doc 0.00167 0.0017 0.0017   **       Kruskal-Wallis

``` r
compare_means(longterm.del.doc ~ Cruise, longterm_delta_stats)
```

    ## # A tibble: 3 x 8
    ##   .y.              group1 group2       p  p.adj p.format p.signif method  
    ##   <chr>            <chr>  <chr>    <dbl>  <dbl> <chr>    <chr>    <chr>   
    ## 1 longterm.del.doc AT34   AT38   0.00130 0.0039 0.0013   **       Wilcoxon
    ## 2 longterm.del.doc AT34   AT39   0.120   0.12   0.1198   ns       Wilcoxon
    ## 3 longterm.del.doc AT38   AT39   0.0165  0.033  0.0165   *        Wilcoxon

``` r
longterm_delta_stats %>% 
  group_by(Cruise) %>% 
  kruskal_test(longterm.del.doc ~ degree_bin)
```

    ## # A tibble: 3 x 7
    ##   Cruise .y.                  n statistic    df      p method        
    ## * <chr>  <chr>            <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 AT34   longterm.del.doc     7      5.14     4 0.273  Kruskal-Wallis
    ## 2 AT38   longterm.del.doc    12      9.22     5 0.101  Kruskal-Wallis
    ## 3 AT39   longterm.del.doc    12      4.15     1 0.0415 Kruskal-Wallis

``` r
longterm_delta_stats %>% 
  group_by(season) %>% 
  summarize_at(vars(longterm.del.doc), list(mean = mean, sd = sd), na.rm = T)  %>% 
  mutate(foldchange = first(mean)/last(mean))
```

    ## # A tibble: 2 x 4
    ##   season  mean    sd foldchange
    ##   <chr>  <dbl> <dbl>      <dbl>
    ## 1 Autumn  5.78  1.35       1.54
    ## 2 Spring  3.76  1.24       1.54

``` r
longterm_delta_stats %>% 
  summarize_at(vars(longterm.del.doc), list(mean = mean, sd = sd), na.rm = T) %>% 
  mutate(cv = sd/mean)
```

    ## # A tibble: 1 x 3
    ##    mean    sd    cv
    ##   <dbl> <dbl> <dbl>
    ## 1  4.54  1.61 0.355

\#CCFs

``` r
bge_summary %>% 
  group_by(Cruise) %>% 
  summarize_at(vars(ccf), list(mean = mean, sd = sd), na.rm = T)
```

    ## # A tibble: 3 x 3
    ##   Cruise  mean    sd
    ##   <chr>  <dbl> <dbl>
    ## 1 AT34    29.4 63.1 
    ## 2 AT38    18.8 17.4 
    ## 3 AT39    20.5  8.25

# BGE

``` r
bge_box <- bge_summary %>% 
  drop_na(bge) %>% 
  ggplot(aes(x = factor(Season, levels = levels), y = bge)) +
  geom_boxplot(width = 0.3) +
  geom_point(aes(fill = factor(degree_bin)), shape = 21, size = 3, alpha = 0.7) +
  theme_classic2(base_size = 16) +
  # scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = matlab.colors) +
  labs(x = "", y = "BGE", fill = "Latitude, ˚N") + 
  stat_compare_means()  
  
  bge_box
```

![](BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
ba_curves / doc_curves2 + bge_box +
   plot_annotation(tag_levels = "a") &
  plot_layout(guides = "collect") +
  theme(plot.tag = element_text(size = 20))
```

![](BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# Experiments at 44˚N

![](BCD_DOC-Bioavailability_files/figure-gfm/Remin%20plots-1.png)<!-- -->

![](BCD_DOC-Bioavailability_files/figure-gfm/Remin%20plots%20initial-1.png)<!-- -->

![](BCD_DOC-Bioavailability_files/figure-gfm/Remin%20plots2-1.png)<!-- -->

# Bioavailability & BGE Tables

``` r
bge.table.data <- bge_summary %>% 
  group_by(Cruise, Season, Station, degree_bin) %>% 
  summarize_at(vars(del.cells:cell_resp), list(mean = mean, sd = sd), na.rm = T)  %>% 
  arrange(factor(Season, levels = levels), degree_bin)

bge.table.data2 <- bge_summary %>% 
  group_by(Cruise, Season) %>% 
  summarize_at(vars(del.cells:cell_resp), list(mean = mean, sd = sd), na.rm = T)  %>% 
  arrange(factor(Season, levels = levels))

bioav.table.data <- export_bioav %>% 
  select(Cruise:Station, degree_bin, end.remin, stationary.harvest, redis_DOC_vol, initial.doc:longterm.ddoc) %>% 
  mutate(shortterm.del.doc = ifelse(shortterm.del.doc < 2, NA, shortterm.del.doc),
         shortterm.bioav.doc = ifelse(is.na(shortterm.del.doc), NA, shortterm.bioav.doc),
         shortterm.ddoc = ifelse(is.na(shortterm.del.doc), NA, shortterm.ddoc),
         
         longterm.del.doc = ifelse(end.remin < 10, NA, longterm.del.doc),
          longterm.bioav.doc = ifelse(is.na( longterm.del.doc), NA,  longterm.bioav.doc),
          longterm.ddoc = ifelse(is.na( longterm.del.doc), NA,longterm.ddoc),
         
         
         persis.doc = ifelse(persis.doc < 2, NA, persis.doc),
         
         percent_persis.doc = ifelse(is.na(persis.doc), NA, percent_persis.doc)) 

station_bioav <- bioav.table.data %>% 
 group_by(Cruise, Season, Station, degree_bin) %>% 
  summarize_at(vars(redis_DOC_vol:longterm.ddoc), list(mean = mean, sd = sd), na.rm = T)  %>% 
  arrange(factor(Season, levels = levels), degree_bin) %>% 
  ungroup()

cruise_bioav <- bioav.table.data %>% 
 group_by(Cruise, Season) %>% 
  summarize_at(vars(redis_DOC_vol:longterm.ddoc), list(mean = mean, sd = sd), na.rm = T)  %>% 
  arrange(factor(Season, levels = levels)) %>% 
  ungroup() 

time <- export_bioav %>% 
  select(Season, Cruise, degree_bin, Station, Bottle, stationary.harvest, end.remin) %>% 
  arrange(factor(Season, levels = levels), degree_bin) %>% 
  ungroup()
```

# Box plots: NPP, BP, BA

<img src="BCD_DOC-Bioavailability_files/figure-gfm/patchwork rate & abundance plots-1.png" style="display: block; margin: auto;" />

# Bar plots: BCD and BCD:NPP

Convert BCD and NPP to mmol C m<sup>-3</sup> d<sup>-1</sup>

``` r
bcd.data <- bcd %>% 
  select(Season, Cruise:degree_bin, Station, int.bcd, int.bcd_global, int.bp, int.NPP, bp.npp, bcd.npp, bcd.npp_global) %>% 
  distinct() %>% 
  mutate_at(vars(int.bcd, int.bcd_global, int.bp, int.NPP), funs(./10^3)) %>% 
  mutate_at(vars(bp.npp, bcd.npp, bcd.npp_global), funs(./10^2)) %>% 
  mutate(degree_bin = as.character(degree_bin)) 

bcd.summary <- bcd.data %>% 
  group_by(Season, Cruise, Station, degree_bin) %>% 
  summarise_at(vars(int.bcd_global, int.bp, int.NPP, bp.npp, bcd.npp_global), list(mean = mean, sd = sd)) %>% 
  arrange(factor(Season, levels = levels), degree_bin) %>% 
  select(Season:degree_bin, contains(c("int.NPP", "int.bp", "bp.npp", "int.bcd", "bcd.npp")))
  
bcd.cruise.summary <- bcd.data %>% 
  group_by(Season, Cruise) %>% 
  summarise_at(vars(int.bcd_global, int.bp, int.NPP, bp.npp, bcd.npp_global), list(cruise_mean = mean, cruise_sd = sd, cruise_max = max, cruise_min = min)) %>% 
  arrange(factor(Season, levels = levels)) %>% 
  select(Season, contains(c("int.NPP", "int.bp", "bp.npp", "int.bcd", "bcd.npp")))

bcd.overall.summary <- bcd.data %>% 
  summarise_at(vars(int.bcd_global, int.bp, int.NPP, bp.npp, bcd.npp_global), list(overall_mean = mean, overall_sd = sd, max = max, min = min))

bcd.npp.summary <- bcd.summary %>% 
  select(Season, Cruise,  degree_bin, int.NPP_mean, int.bcd_global_mean) %>% 
  rename(NPP = int.NPP_mean,
         BCD = int.bcd_global_mean) %>% 
  pivot_longer(c(NPP, BCD), names_to = "rate", values_to = "ave") %>% 
  left_join(., bcd.summary %>% 
  select(Season, Cruise, degree_bin,  int.NPP_sd, int.bcd_global_sd) %>% 
    rename(NPP = int.NPP_sd,
         BCD = int.bcd_global_sd) %>% 
  pivot_longer(c(NPP, BCD), names_to = "rate", values_to = "sd"))
```

    ## Adding missing grouping variables: `Station`
    ## Adding missing grouping variables: `Station`

    ## Joining, by = c("Station", "Season", "Cruise", "degree_bin", "rate")

``` r
bcd.npp_summary_table <- left_join(bcd.summary, bcd.cruise.summary) %>% 
  mutate_at(vars(contains(c("int.NPP", "int.bp", "int.bcd"))), round, 3) %>% 
  mutate_at(vars(contains(c("bp.npp", "bcd.npp"))), round, 2) 
```

    ## Joining, by = "Season"

``` r
export_bar <- export %>% 
  mutate(degree_bin = as.character(degree_bin)) %>% 
  group_by(Season, Cruise, degree_bin) %>% 
  summarise_at(vars(int_delta_DOC_100, doc_ncp_100, NCP_mol_100), list(mean = mean, sd = sd)) %>% 
  arrange(factor(Season, levels = levels), degree_bin)
```

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
    ## 1 AT34   Late Spring   6.13  2.83
    ## 2 AT38   Early Autumn 13.1   3.34
    ## 3 AT39   Early Spring  2.8   1.37

![](BCD_DOC-Bioavailability_files/figure-gfm/BCD%20and%20NPP%20bar%20plots-1.png)<!-- -->

![](BCD_DOC-Bioavailability_files/figure-gfm/BCD%20and%20NPP%20bar%20plots2-1.png)<!-- -->

# Regressions: Property-Property

``` r
summary(bcd.data$int.NPP)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1142  0.2328  0.6400  0.7681  0.9524  3.3406

## BCD v NPP

    ## RMA was not requested: it will not be computed.

## BCD:NPP v NPP

### exponential model using lm function

``` r
reg2 <- lm(log(bcd.npp_global) ~ int.NPP, data = bcd.data)

reg2.df <- data.frame(x = bcd.data$int.NPP, y = exp(fitted(reg2)))
```

### manual exponential model

y = alpha \* (e^beta \* x) + theta

    ## $alpha
    ## (Intercept) 
    ##    0.312114 
    ## 
    ## $beta
    ##    int.NPP 
    ## -0.4936019 
    ## 
    ## $theta
    ## [1] 0.04601972

``` r
model <- nls(bcd.npp_global ~ alpha * exp(beta * int.NPP) + theta , data = bcd.data, start = start)

model.df <- data.frame(x = bcd.data$int.NPP, y = predict(model, list(x = bcd.data$int.NPP)) )

summary(model)
```

    ## 
    ## Formula: bcd.npp_global ~ alpha * exp(beta * int.NPP) + theta
    ## 
    ## Parameters:
    ##       Estimate Std. Error t value Pr(>|t|)    
    ## alpha  1.94713    0.54626   3.564   0.0012 ** 
    ## beta  -8.89252    1.90546  -4.667 5.57e-05 ***
    ## theta  0.20132    0.02245   8.969 4.03e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1012 on 31 degrees of freedom
    ## 
    ## Number of iterations to convergence: 11 
    ## Achieved convergence tolerance: 3.152e-06

## Plots

![](BCD_DOC-Bioavailability_files/figure-gfm/BCD%20and%20NPP%20regression%20plots-1.png)<!-- -->

``` r
compare_means(bcd.npp_global_mean ~ Cruise, bcd.summary, method = "kruskal.test")
```

    ## Adding missing grouping variables: `Season`, `Station`

    ## # A tibble: 1 x 6
    ##   .y.                       p  p.adj p.format p.signif method        
    ##   <chr>                 <dbl>  <dbl> <chr>    <chr>    <chr>         
    ## 1 bcd.npp_global_mean 0.00411 0.0041 0.0041   **       Kruskal-Wallis

``` r
compare_means(bcd.npp_global_mean ~ Season, bcd.summary)
```

    ## Adding missing grouping variables: `Cruise`, `Station`

    ## # A tibble: 6 x 8
    ##   .y.             group1      group2           p p.adj p.format p.signif method 
    ##   <chr>           <chr>       <chr>        <dbl> <dbl> <chr>    <chr>    <chr>  
    ## 1 bcd.npp_global… Early Spri… Late Spri… 0.762   0.76  0.7619   ns       Wilcox…
    ## 2 bcd.npp_global… Early Spri… Early Aut… 0.0667  0.2   0.0667   ns       Wilcox…
    ## 3 bcd.npp_global… Early Spri… Early Win… 0.00952 0.048 0.0095   **       Wilcox…
    ## 4 bcd.npp_global… Late Spring Early Aut… 0.132   0.26  0.1320   ns       Wilcox…
    ## 5 bcd.npp_global… Late Spring Early Win… 0.00216 0.013 0.0022   **       Wilcox…
    ## 6 bcd.npp_global… Early Autu… Early Win… 0.0411  0.16  0.0411   *        Wilcox…

## bp/bcd profile

![](BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

![](BCD_DOC-Bioavailability_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->
