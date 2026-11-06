# =============================================================================
# Cheng et al. (2024) Watkins landraces: build a master table and shortlist
# one or two candidate landraces per ancestral group (AG1-AG7)
#
# Paper:  Cheng et al. 2024, Nature 632:823. doi:10.1038/s41586-024-07682-9
# Input:  Raw_data/Tables_1-10.xlsx   (Supplementary Tables 1-10)
#         Raw_data/Tables_11-20.xlsx  (Supplementary Tables 11-20)
#         Raw_data/Tables_22-33.xlsx  (Supplementary Tables 21-33)
#
# Sections 1-6: master table, exploration plots, candidate shortlist
# Sections 7-9: take-all - GWAS loci, haplotypes by AG, overlapping Paragon NILs
# Section 10:   each line's own take-all haplotypes (downloaded from wwwg2b.com;
#               needs an internet connection the first time it runs)
# Section 11:   is the take-all score confounded by flowering time / height?
#               (JIC field phenotype file, downloaded from Earlham Open Data)
# Output: Outputs/ (CSV tables and PNG plots)
#
# =============================================================================


# ---- 0. Setup ---------------------------------------------------------------

# Run this line to install the packages if needed
# install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "stringr", "readr", "maps", "jsonlite"))

library(readxl)   # read Excel files
library(dplyr)    # filter / select / mutate / join / summarise
library(tidyr)    # pivot_longer (wide -> long)
library(ggplot2)  # plots
library(stringr)  # text cleaning
library(readr)    # write_csv
library(jsonlite) # read data downloaded from the WWWG2B portal (section 10)

file_1_10  <- "Raw_data/Tables_1-10.xlsx"
file_11_20 <- "Raw_data/Tables_11-20.xlsx"
file_22_33 <- "Raw_data/Tables_22-33.xlsx"
dir.create("Outputs", showWarnings = FALSE)

# --- Shortlist settings: change these and re-run section 6 -------------------
min_membership <- 0.80        # min ADMIXTURE proportion in the line's own AG
keep_habits    <- c("spring") # growth habits allowed; c("spring", "winter") for both
n_per_ag       <- 3           # candidates to keep per AG

# --- Take-all settings (sections 8-9) -----------------------------------------
min_hap_carriers <- 5         # ignore haplotypes carried by fewer landraces than this
max_nil_mb       <- 100       # "tight" NIL = marker-selected segment shorter than this (Mb)

# One fixed colour per group (Okabe-Ito colour-blind-safe palette)
ag_colours <- c(AG1 = "#E69F00", AG2 = "#56B4E9", AG3 = "#009E73",
                AG4 = "#F0E442", AG5 = "#0072B2", AG6 = "#D55E00",
                AG7 = "#CC79A7", Modern = "#999999")

# Significance stars for BH-corrected q-values (used in the stats sections)
stars <- function(q) case_when(is.na(q) ~ "", q < 0.001 ~ "***", q < 0.01 ~ "**",
                               q < 0.05 ~ "*", TRUE ~ "ns")
# Format a p-value for plot labels ("p < 0.001" or "p = 0.03")
fmt_p <- function(p, label = "p") ifelse(p < 0.001, paste(label, "< 0.001"),
                                          paste(label, "=", signif(p, 2)))
# NB on p-values in this script: landraces within an AG are related, and the
# take-all haplotype effects (Table 22) were estimated by GWAS on these same
# landraces. Tests of take-all differences between groups are therefore not
# fully independent evidence - treat them as supporting, not as proof.


# ---- 1. Table 1: germplasm description (all 1,051 lines) --------------------
# Row 1 of the sheet is the table title, so skip it; row 2 holds the headers.
# "/" is used for missing values in this table, so treat it as NA.

table1_raw <- read_excel(file_1_10, sheet = "S1", skip = 1,
                         na = c("", "/", "NA"))

table1 <- table1_raw %>%
  select(accession     = `Accession name`,
         wwwg2b_id     = `Accession name in China (WWWG2B)`,
         gru_storcode  = `GRU StorCode`,
         country       = `Country of origin`,
         latitude      = Latitude,
         longitude     = Longitude,
         growth_habit  = `Growth habit`,
         group         = `Group in this study`,
         seedstor_url  = starts_with("Accession URL")) %>%
  mutate(
    accession    = str_remove(accession, "\\*$") %>% str_trim(), # "Paragon*" -> "Paragon"
    latitude     = as.numeric(latitude),
    longitude    = as.numeric(longitude),
    growth_habit = str_to_lower(growth_habit),                   # "Spring" -> "spring"
    country      = na_if(country, "unknown"),
    is_watkins   = str_starts(accession, "WATDE")                 # Watkins IDs start with WATDE
  )

# Quick checks
glimpse(table1)
count(table1, group)          # AG1-AG7 + Modern
count(table1, is_watkins)     # expect 827 Watkins (Chinese Spring is placed in AG1 but isn't Watkins)


# ---- 2. Table 4: ADMIXTURE ancestry proportions -----------------------------
# Row 1 = title, row 2 = headers ("K=2", "K=3", ...). Each block of columns is
# one value of K, separated by an empty column. The paper uses K = 7, which is
# columns 28-34. These are in AG1 ... AG7 order (checked below).

table4_raw <- read_excel(file_1_10, sheet = "S4", skip = 2, col_names = FALSE,
                         .name_repair = "unique_quiet")
table4_hdr <- read_excel(file_1_10, sheet = "S4", range = cell_rows(2), col_names = FALSE,
                         .name_repair = "unique_quiet")

stopifnot(all(unlist(table4_hdr[28:34]) == "K=7"))  # stops if the columns ever move

admixture <- table4_raw %>%
  select(1, 2, 28:34) %>%
  setNames(c("accession", "ag_admixture", paste0("AG", 1:7))) %>%
  filter(!is.na(accession)) %>%
  mutate(accession = str_replace(accession, "_", " "),  # "Chinese_Spring" -> "Chinese Spring"
         across(AG1:AG7, as.numeric))

# Long format: one row per accession x AG
admixture_long <- admixture %>%
  pivot_longer(AG1:AG7, names_to = "component", values_to = "proportion")

# Sanity check: each group's average membership should be highest in its own
# component (AG1 lines highest in AG1, etc.). Expect a strong diagonal.
admixture_long %>%
  group_by(ag_admixture, component) %>%
  summarise(mean_prop = round(mean(proportion), 2), .groups = "drop") %>%
  pivot_wider(names_from = component, values_from = mean_prop)

# Per accession: proportion in its own AG, plus its next-largest AG
membership <- admixture_long %>%
  group_by(accession) %>%
  summarise(
    own_membership = proportion[component == first(ag_admixture)],
    second_ag      = component[component != first(ag_admixture)][
                       which.max(proportion[component != first(ag_admixture)])],
    second_prop    = max(proportion[component != first(ag_admixture)]),
    .groups = "drop"
  )


# ---- 3. Table 12: Paragon x Watkins RIL (mapping) populations ---------------
# 109 populations; the 73 marked "Y" were used in the paper. A few Watkins
# lines are parents of more than one population, so summarise per parent.

table12 <- read_excel(file_11_20, sheet = "S12", skip = 1) %>%
  select(used_in_study  = `Used in this study`,
         ril_population = `RIL population`,
         common_parent  = `Common NAM parent`,
         accession      = starts_with("Non-common Watkins parent"),
         core_set       = `Watkins Core Set Line`,
         n_rils         = starts_with("Number of"))

ril_parents <- table12 %>%
  group_by(accession) %>%
  summarise(
    paragon_ril       = any(common_parent == "Paragon"),
    ril_used_in_paper = any(used_in_study == "Y", na.rm = TRUE),
    ril_populations   = paste(ril_population, collapse = "; "),
    ril_crosses       = paste(common_parent, "x", accession, collapse = "; "),
    n_rils            = sum(n_rils, na.rm = TRUE),
    in_core_set       = any(!is.na(core_set)),
    .groups = "drop"
  )


# ---- 3b. Table 31: Paragon near-isogenic lines (NILs) -----------------------
# 738 NILs in 143 families. Each family comes from one Watkins donor and has
# sibling lines carrying either the Watkins (W) or the Paragon (P) allele of a
# targeted segment - so the P siblings are perfect controls for the W siblings.
#   MASstart/MASend     = the marker-selected (targeted) segment
#   startIntro/endIntro = the whole introgressed region (often much larger)

table31 <- read_excel(file_22_33, sheet = "S31", skip = 1) %>%
  filter(!is.na(NILFamily)) %>%                       # drops the Paragon row
  mutate(across(c(MASstart, MASend, startIntro, endIntro), as.numeric),
         allele_type = case_when(str_starts(Allele, "W")      ~ "Watkins",
                                 str_starts(Allele, "[Pp]")   ~ "Paragon",
                                 TRUE                         ~ "Other"))

nil_families <- table31 %>%
  group_by(NILFamily) %>%
  summarise(donor        = first(Parent),
            toolkit_set  = first(TK),
            chr          = first(Chr),
            MASstart     = first(MASstart),
            MASend       = first(MASend),
            mas_size_mb  = (first(MASend) - first(MASstart)) / 1e6,
            target_trait = first(Crop.Ontology.Variable.name),
            n_W_lines    = sum(allele_type == "Watkins"),
            n_P_lines    = sum(allele_type == "Paragon"),
            W_gru_codes  = paste(GRU[allele_type == "Watkins"], collapse = "; "),
            P_gru_codes  = paste(GRU[allele_type == "Paragon"], collapse = "; "),
            .groups = "drop")

# Number of NIL families per Watkins donor (joined to the master table below)
nil_donors <- nil_families %>%
  count(donor, name = "n_nil_families") %>%
  rename(accession = donor)


# ---- 4. Master table: everything joined by accession ------------------------

master <- table1 %>%
  left_join(membership,  by = "accession") %>%
  left_join(ril_parents, by = "accession") %>%
  left_join(nil_donors,  by = "accession") %>%
  mutate(paragon_ril       = coalesce(paragon_ril, FALSE),
         ril_used_in_paper = coalesce(ril_used_in_paper, FALSE),
         n_nil_families    = coalesce(n_nil_families, 0L),
         admixed           = own_membership < min_membership)

# Every landrace should have found its ancestry values, and every RIL parent
# should be a Watkins line in Table 1 (both should print 0)
sum(is.na(master$own_membership[master$group != "Modern"]))
sum(!ril_parents$accession %in% table1$accession)

write_csv(master, "Outputs/watkins_master_table.csv")

# Your elite lines and the reference genome line
master %>%
  filter(accession %in% c("Paragon", "Hereward", "Soissons", "Chinese Spring")) %>%
  select(accession, group, country, growth_habit)


# ---- 5. Explore the collection ----------------------------------------------

watkins <- master %>% filter(is_watkins)

# 5a. Lines per AG, split by growth habit
ag_summary <- watkins %>%
  group_by(group) %>%
  summarise(n_lines         = n(),
            n_spring        = sum(growth_habit == "spring", na.rm = TRUE),
            n_winter        = sum(growth_habit == "winter", na.rm = TRUE),
            n_clear_member  = sum(own_membership >= min_membership),
            n_paragon_ril   = sum(paragon_ril),
            n_countries     = n_distinct(country, na.rm = TRUE),
            .groups = "drop")
ag_summary
write_csv(ag_summary, "Outputs/ag_summary.csv")

# 5-stats. Do the AGs differ in growth habit and in ancestry "purity"?
habit_table <- table(watkins$group, coalesce(watkins$growth_habit, "unknown"))
habit_table
# Fisher's exact test (simulated p, because many cells are small)
habit_test <- fisher.test(habit_table, simulate.p.value = TRUE, B = 1e4)
habit_test
# Own-AG ancestry proportion by AG (Kruskal-Wallis, non-parametric)
membership_kw <- kruskal.test(own_membership ~ factor(group), data = watkins)
membership_kw

p_counts <- ggplot(watkins, aes(x = group, fill = growth_habit)) +
  geom_bar() +
  labs(x = "Ancestral group", y = "Number of landraces", fill = "Growth habit",
       title = "Watkins landraces per ancestral group",
       subtitle = paste0("Growth habit differs between AGs (Fisher's exact test, ",
                         fmt_p(habit_test$p.value), ")")) +
  theme_bw()
p_counts
ggsave("Outputs/01_landraces_per_AG.png", p_counts, width = 7, height = 4.5, dpi = 300)

# 5b. Top countries of origin in each AG
top_countries <- watkins %>%
  count(group, country, sort = TRUE) %>%
  group_by(group) %>%
  slice_max(n, n = 5, with_ties = FALSE) %>%
  arrange(group, desc(n))
print(top_countries, n = Inf)
write_csv(top_countries, "Outputs/top_countries_per_AG.csv")

# 5c. Map of collection sites. Many lines share a country-level coordinate
# (e.g. every Australian line plots at Canberra), so points are jittered a little.
world <- map_data("world")

p_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_jitter(data = filter(watkins, !is.na(latitude)),
              aes(longitude, latitude, colour = group),
              width = 1.5, height = 1.5, size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = ag_colours) +
  coord_quickmap(xlim = c(-20, 150), ylim = c(-45, 65)) +
  labs(x = NULL, y = NULL, colour = "AG",
       title = "Collection sites of Watkins landraces, coloured by ancestral group") +
  theme_bw()
p_map
ggsave("Outputs/02_collection_map.png", p_map, width = 10, height = 6, dpi = 300)

# 5d. How "pure" is each line? Own-AG membership, with the shortlist threshold
p_membership <- ggplot(watkins, aes(x = group, y = own_membership, colour = group)) +
  geom_boxplot(outlier.shape = NA, colour = "grey40") +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  geom_hline(yintercept = min_membership, linetype = "dashed") +
  scale_colour_manual(values = ag_colours, guide = "none") +
  labs(x = "Ancestral group", y = "Proportion of ancestry in own AG (K = 7)",
       title = "ADMIXTURE membership of each landrace in its assigned group",
       subtitle = paste0("No significant difference between AGs (Kruskal-Wallis, ",
                         fmt_p(membership_kw$p.value), ")")) +
  theme_bw()
p_membership
ggsave("Outputs/03_membership_by_AG.png", p_membership, width = 7, height = 4.5, dpi = 300)

# 5e. Classic ADMIXTURE bar plot: one thin bar per landrace, sorted by AG
#     then by membership (similar to Cheng et al. Extended Data Fig. 2)
bar_order <- watkins %>% arrange(group, desc(own_membership)) %>% pull(accession)

p_admix <- admixture_long %>%
  filter(accession %in% bar_order) %>%
  mutate(accession = factor(accession, levels = bar_order)) %>%
  ggplot(aes(x = accession, y = proportion, fill = component)) +
  geom_col(width = 1) +
  facet_grid(~ ag_admixture, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = ag_colours) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Component",
       title = "ADMIXTURE (K = 7) for 827 Watkins landraces") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.spacing = unit(0.1, "lines"))
p_admix
ggsave("Outputs/04_admixture_barplot.png", p_admix, width = 12, height = 4, dpi = 300)


# ---- 6. Shortlist candidates per AG -----------------------------------------
# Priority order within each AG:
#   1. Parent of a Paragon x Watkins RIL population (allows later gene mapping)
#   2. That population was used in the paper (phenotyped, genotyped)
#   3. Highest membership of its own AG (a clean representative of the group)
# Filters: growth habit in keep_habits, own membership >= min_membership.

# How many lines survive each filter? If an AG drops to 0, relax the settings.
filter_funnel <- watkins %>%
  group_by(group) %>%
  summarise(all_lines        = n(),
            right_habit      = sum(growth_habit %in% keep_habits),
            plus_membership  = sum(growth_habit %in% keep_habits &
                                     own_membership >= min_membership),
            plus_paragon_ril = sum(growth_habit %in% keep_habits &
                                     own_membership >= min_membership & paragon_ril),
            .groups = "drop")
filter_funnel

shortlist <- watkins %>%
  filter(growth_habit %in% keep_habits,
         own_membership >= min_membership) %>%
  arrange(group, desc(paragon_ril), desc(ril_used_in_paper), desc(own_membership)) %>%
  group_by(group) %>%
  slice_head(n = n_per_ag) %>%
  ungroup() %>%
  select(group, accession, country, growth_habit,
         own_membership, second_ag, second_prop,
         paragon_ril, ril_used_in_paper, ril_populations, n_rils,
         n_nil_families, in_core_set, gru_storcode, seedstor_url)

print(shortlist, n = Inf, width = Inf)
write_csv(shortlist, "Outputs/candidate_shortlist.csv")


# ---- 7. Take-all GWAS loci (Table 18) ---------------------------------------
# Take-all incidence was scored on the Watkins landraces at Rothamsted
# 2008-2012 (trait code "GgtriInc"). GWAS found 58 marker-trait associations
# (MTAs). Overlapping MTAs are grouped into regions; a region found in more
# than one YEAR is more convincing than one found in a single year.

takeall_mtas <- read_excel(file_11_20, sheet = "S18", skip = 1) %>%
  filter(`New Trait abbreviation` == "GgtriInc") %>%
  transmute(environment = `Environments and Years`,
            year        = as.integer(str_extract(environment, "20[0-9]{2}")),
            chr         = str_remove(Chromosome, "chr"),
            start       = as.numeric(Start),
            end         = as.numeric(End),
            n_genes     = `# of gene models annotated`,
            prioritised = `priotized (y/n)` == "y",
            top_snp     = `The most significant SNP`,
            p_value     = as.numeric(`p-value`),
            pve_pct     = as.numeric(`pve (%)`))

nrow(takeall_mtas)   # expect 58

# Group overlapping MTAs on the same chromosome into regions
takeall_regions <- takeall_mtas %>%
  arrange(chr, start) %>%
  group_by(chr) %>%
  mutate(region_no = cumsum(coalesce(start > lag(cummax(end)), TRUE))) %>%
  group_by(chr, region_no) %>%
  summarise(start_mb    = round(min(start) / 1e6, 1),
            end_mb      = round(max(end) / 1e6, 1),
            n_mtas      = n(),
            years       = paste(sort(unique(year)), collapse = ", "),
            n_years     = n_distinct(year),
            prioritised = any(prioritised),
            max_pve_pct = max(pve_pct),
            .groups = "drop") %>%
  mutate(recurring = n_years > 1,
         region    = paste0(chr, ":", start_mb, "-", end_mb, " Mb"))

takeall_regions %>% filter(n_mtas > 1)   # recurring = found in >1 year
write_csv(takeall_regions, "Outputs/takeall_regions.csv")


# ---- 8. Take-all haplotypes by ancestral group (Table 22) -------------------
# At each take-all MTA, Table 22 lists every haplotype, its estimated effect on
# take-all (as % of the trial mean; NEGATIVE = LESS take-all = protective), and
# how many landraces in each AG (and modern cultivars) carry it.
# NB: it gives COUNTS per AG only - not which individual landrace carries
# which haplotype - so this is a group-level comparison.

takeall_haps <- read_excel(file_22_33, sheet = "S22", skip = 1) %>%
  filter(`New Trait abbreviation` == "GgtriInc") %>%
  transmute(environment = MTA,
            chr         = str_remove(chromosome, "chr"),
            start       = as.numeric(`peak start`),
            end         = as.numeric(`peak end`),
            block       = sig_block_id,
            haplotype   = haplotype_id,
            effect_pct  = 100 * as.numeric(`effect in percentage of the mean`),
            across(c(AG1:AG7, landrance, modern), as.numeric)) %>%
  mutate(locus = paste0(chr, ":", round(start / 1e6), "-", round(end / 1e6), " Mb (",
                        str_extract(environment, "20[0-9]{2}"), ")"))

n_distinct(takeall_haps$locus)   # expect 58

# Long format: one row per haplotype x group (AG1-AG7, modern)
haps_long <- takeall_haps %>%
  mutate(n_modern_carriers = modern) %>%      # keep a copy for the filter in 8d
  pivot_longer(c(AG1:AG7, modern), names_to = "group", values_to = "n_lines") %>%
  mutate(group = recode(group, modern = "Modern"))

# 8a. Average effect of the haplotypes each group carries at each locus
#     (weighted by how many lines carry each haplotype)
locus_by_group <- haps_long %>%
  group_by(locus, chr, start, group) %>%
  summarise(mean_effect_pct = sum(n_lines * effect_pct) / sum(n_lines),
            .groups = "drop")

# 8b. One number per group: average across all 58 loci
takeall_ag_score <- locus_by_group %>%
  group_by(group) %>%
  summarise(loci_protective      = sum(mean_effect_pct < 0, na.rm = TRUE),
            loci_total           = sum(!is.na(mean_effect_pct)),
            avg_effect_pct       = round(mean(mean_effect_pct, na.rm = TRUE), 1),
            .groups = "drop")
takeall_ag_score
write_csv(takeall_ag_score, "Outputs/takeall_score_by_AG.csv")

# 8b-stats. Is each AG's take-all haplotype profile different from modern wheat?
#   The 58 loci are used as paired blocks: at each locus, compare the AG's
#   average effect with modern wheat's (paired Wilcoxon, BH across 7 AGs).
#   Friedman test = do the 8 groups differ overall?
#   NB: nearby loci (e.g. 3 x 6A in 2012) aren't fully independent.
locus_wide <- locus_by_group %>%
  select(locus, group, mean_effect_pct) %>%
  pivot_wider(names_from = group, values_from = mean_effect_pct) %>%
  drop_na()
friedman_res <- friedman.test(as.matrix(select(locus_wide, -locus)))
friedman_res

ag_vs_modern <- tibble(group = paste0("AG", 1:7)) %>%
  rowwise() %>%
  mutate(loci                   = nrow(locus_wide),
         median_diff_vs_modern  = round(median(locus_wide[[group]] - locus_wide$Modern), 1),
         loci_higher_than_modern = sum(locus_wide[[group]] > locus_wide$Modern),
         p = wilcox.test(locus_wide[[group]], locus_wide$Modern, paired = TRUE)$p.value) %>%
  ungroup() %>%
  mutate(q = p.adjust(p, "BH"), sig = stars(q))
ag_vs_modern     # positive diff = more take-all-increasing haplotypes than modern wheat
write_csv(ag_vs_modern, "Outputs/stats_takeall_AG_vs_modern.csv")

# 8c. Heatmap: rows = take-all loci (in genome order), columns = groups
locus_order <- locus_by_group %>% distinct(locus, chr, start) %>%
  arrange(chr, start) %>% pull(locus)

p_takeall <- locus_by_group %>%
  mutate(locus = factor(locus, levels = rev(locus_order)),
         group = factor(group, levels = c(paste0("AG", 1:7), "Modern"))) %>%
  ggplot(aes(x = group, y = locus, fill = mean_effect_pct)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-60, 60), oob = scales::squish,
                       name = "Take-all effect\n(% of trial mean)") +
  labs(x = NULL, y = "Take-all GWAS locus (year)",
       title = "Take-all haplotypes carried by each group",
       subtitle = paste0("Blue = reduces take-all, red = increases it. Groups differ (Friedman ",
                         fmt_p(friedman_res$p.value), ").\nStars under each AG = different from modern ",
                         "wheat across the 58 loci (paired Wilcoxon, BH)")) +
  scale_x_discrete(labels = function(g) ifelse(
    g == "Modern", "Modern\n(reference)",
    paste0(g, "\n", ag_vs_modern$sig[match(g, ag_vs_modern$group)]))) +
  theme_minimal(base_size = 9)
p_takeall
ggsave("Outputs/05_takeall_haplotypes_heatmap.png", p_takeall, width = 7, height = 10, dpi = 300)

# 8d. Protective haplotypes found ONLY in Watkins (absent from modern wheat):
#     a source of take-all variation that breeders don't already have.
#     Divided by group size to give "per landrace".
ag_sizes <- watkins %>% count(group, name = "n_landraces")

watkins_only_protective <- haps_long %>%
  filter(group != "Modern",
         effect_pct < 0, n_modern_carriers == 0, landrance >= min_hap_carriers) %>%
  group_by(group) %>%
  summarise(carriers = sum(n_lines), .groups = "drop") %>%
  left_join(ag_sizes, by = "group") %>%
  mutate(per_landrace = round(carriers / n_landraces, 1))
watkins_only_protective
write_csv(watkins_only_protective, "Outputs/takeall_watkins_only_protective.csv")

p_protective <- ggplot(watkins_only_protective,
                       aes(x = group, y = per_landrace, fill = group)) +
  geom_col() +
  scale_fill_manual(values = ag_colours, guide = "none") +
  labs(x = "Ancestral group", y = "Watkins-only protective haplotypes per landrace",
       title = "Take-all-reducing haplotypes that are absent from modern wheat") +
  theme_bw()
p_protective
ggsave("Outputs/06_takeall_watkins_only_protective.png", p_protective, width = 7, height = 4.5, dpi = 300)


# ---- 9. Paragon NILs whose targeted segment covers a take-all locus ---------
# NILs were made for other traits (height, heading date, yield...) and were
# NOT tested for take-all. But if a NIL's Watkins segment covers a take-all
# locus, the W vs P siblings could be compared for take-all / root microbiome
# with the rest of the genome identical.

nil_takeall <- nil_families %>%
  inner_join(takeall_mtas, by = "chr", relationship = "many-to-many") %>%
  filter(MASstart <= end, MASend >= start) %>%                # segments overlap
  left_join(select(master, accession, group), by = c("donor" = "accession")) %>%
  group_by(NILFamily, donor, group, chr, MASstart, MASend, mas_size_mb,
           target_trait, n_W_lines, n_P_lines, W_gru_codes, P_gru_codes) %>%
  summarise(takeall_mtas  = n(),
            takeall_years = paste(sort(unique(year)), collapse = ", "),
            .groups = "drop") %>%
  # flag NILs that cover a recurring (multi-year) take-all region
  left_join(takeall_regions %>% filter(recurring) %>%
              transmute(chr, rec_start = start_mb * 1e6, rec_end = end_mb * 1e6, region),
            by = "chr", relationship = "many-to-many") %>%
  mutate(covers_recurring = !is.na(region) & MASstart <= rec_end & MASend >= rec_start) %>%
  group_by(across(NILFamily:takeall_years)) %>%
  summarise(recurring_region = paste(region[covers_recurring], collapse = "; "),
            .groups = "drop") %>%
  arrange(mas_size_mb)

nrow(nil_takeall)   # families overlapping any take-all MTA (most segments are large)
write_csv(nil_takeall, "Outputs/NILs_overlapping_takeall_all.csv")

# Tight segments only: the most useful for testing a specific locus
nil_takeall_tight <- nil_takeall %>% filter(mas_size_mb < max_nil_mb)
print(nil_takeall_tight, n = Inf, width = Inf)
write_csv(nil_takeall_tight, "Outputs/NILs_overlapping_takeall_tight.csv")

# Which shortlisted candidates also donated NILs?
shortlist %>% filter(n_nil_families > 0) %>% select(group, accession, n_nil_families)


# ---- 10. Each line's own take-all haplotypes (WWWG2B portal) ----------------
# Table 22 only gives haplotype COUNTS per AG. The WWWG2B portal lists which
# lines carry each haplotype in each block, so here we:
#   1. download the haplotype carriers for every take-all block (saved in
#      Raw_data/wwwg2b_cache/ so they're only downloaded once)
#   2. find the haplotype each line carries at each of the 58 take-all loci
#   3. add the haplotype effects from Table 22
#   4. average them into a take-all haplotype score per line
#      (NEGATIVE = carries haplotypes linked to LESS take-all)
# NB: this is a GENETIC PREDICTION from GWAS effects (one site, small effects,
# assumes effects add up) - not a measured take-all score.

cache_dir  <- "Raw_data/wwwg2b_cache"
dir.create(cache_dir, showWarnings = FALSE)
portal_url <- paste0("https://wwwg2b.com/api/tools/breeding/haploblock/",
                     "get_block_sample_cluster_info?type=block&block=")

get_block_carriers <- function(block) {
  cache_file <- file.path(cache_dir, paste0(block, ".json"))
  if (!file.exists(cache_file)) {                 # download only if not cached
    download.file(paste0(portal_url, block), cache_file, quiet = TRUE, mode = "wb")
    Sys.sleep(0.5)                                # be polite to the server
  }
  res <- fromJSON(cache_file)
  if (res$code != 0) stop("Portal returned an error for block ", block)
  carriers <- res$clusterData                     # list: haplotype -> line names
  tibble(block     = block,
         haplotype = rep(names(carriers), lengths(carriers)),
         line      = unlist(carriers, use.names = FALSE))
}

take_blocks <- unique(takeall_haps$block)
length(take_blocks)                               # 56 blocks (2 blocks have 2 loci)

hap_carriers <- bind_rows(lapply(take_blocks, get_block_carriers)) %>%
  mutate(line_key = str_to_lower(str_remove_all(line, "[^A-Za-z0-9]")))  # "Maris Huntsman" = "MarisHuntsman"

# Check names match Table 1 (should be few or none unmatched)
master <- master %>% mutate(line_key = str_to_lower(str_remove_all(accession, "[^A-Za-z0-9]")))
unmatched <- setdiff(unique(hap_carriers$line_key), master$line_key)
length(unmatched); head(unmatched, 20)

# Haplotype effects from Table 22 (block + haplotype number identify a haplotype)
hap_effects <- takeall_haps %>%
  transmute(locus, block, haplotype = as.character(haplotype), effect_pct,
            watkins_only_protective = effect_pct < 0 & modern == 0)

# One row per line x locus: the effect of the haplotype that line carries
line_locus <- hap_carriers %>%
  inner_join(hap_effects, by = c("block", "haplotype"), relationship = "many-to-many")

line_scores <- line_locus %>%
  group_by(line_key) %>%
  summarise(takeall_loci                    = n(),
            takeall_score                   = round(mean(effect_pct), 1),
            takeall_loci_protective         = sum(effect_pct < 0),
            takeall_watkins_only_protective = sum(watkins_only_protective),
            .groups = "drop")

master_takeall <- master %>%
  left_join(line_scores, by = "line_key") %>%
  select(accession, group, country, growth_habit, own_membership,
         paragon_ril, n_nil_families, starts_with("takeall_"))
write_csv(master_takeall, "Outputs/takeall_line_scores.csv")

# 10a. Sanity check: per-line scores averaged by group should match section 8
master_takeall %>%
  group_by(group) %>%
  summarise(n_scored = sum(!is.na(takeall_score)),
            mean_line_score = round(mean(takeall_score, na.rm = TRUE), 1),
            .groups = "drop")

# 10b. Your elite lines
master_takeall %>%
  filter(accession %in% c("Paragon", "Hereward", "Soissons")) %>%
  select(accession, starts_with("takeall_"))

# 10b-stats. Do groups differ in per-line take-all score?
scored   <- filter(master_takeall, !is.na(takeall_score))
kw_score <- kruskal.test(takeall_score ~ factor(group), data = scored)
kw_score
epsilon2 <- unname(kw_score$statistic) / (nrow(scored) - 1)   # effect size, 0-1
round(epsilon2, 2)                                            # > 0.26 = large
pairwise_score <- pairwise.wilcox.test(scored$takeall_score, scored$group,
                                       p.adjust.method = "BH")
pairwise_score
write_csv(as_tibble(pairwise_score$p.value, rownames = "group"),
          "Outputs/stats_takeall_score_pairwise.csv")

# Watkins-only protective haplotypes per line: do AGs differ?
protective_kw <- kruskal.test(takeall_watkins_only_protective ~ factor(group),
                              data = filter(scored, group != "Modern"))
protective_kw
# Re-save plot 06 (section 8d) with this result added
p_protective <- p_protective +
  labs(subtitle = paste0("Per-landrace counts differ between AGs (Kruskal-Wallis, ",
                         fmt_p(protective_kw$p.value), ")"))
ggsave("Outputs/06_takeall_watkins_only_protective.png", p_protective, width = 7, height = 4.5, dpi = 300)

# Where does each candidate sit within its own AG? (0 = most protective line)
candidate_percentiles <- scored %>%
  group_by(group) %>%
  mutate(percentile_in_group = round(100 * (rank(takeall_score) - 1) / (n() - 1))) %>%
  ungroup() %>%
  filter(accession %in% c(shortlist$accession, "Paragon", "Hereward", "Soissons")) %>%
  select(group, accession, takeall_score, percentile_in_group) %>%
  arrange(group, takeall_score)
print(candidate_percentiles, n = Inf)
write_csv(candidate_percentiles, "Outputs/stats_candidate_percentiles.csv")

# 10c. The shortlist with each candidate's take-all haplotype score
shortlist_takeall <- shortlist %>%
  left_join(select(master_takeall, accession, starts_with("takeall_")), by = "accession")
shortlist_takeall %>%
  select(group, accession, paragon_ril, n_nil_families, starts_with("takeall_")) %>%
  print(n = Inf, width = Inf)
write_csv(shortlist_takeall, "Outputs/candidate_shortlist_takeall.csv")

# 10d. Within each AG: the lines with the LOWEST and HIGHEST scores that pass
#      the same filters as the shortlist (habit + membership). A low/high pair
#      from one AG gives a take-all contrast with similar ancestry.
takeall_extremes <- master_takeall %>%
  filter(str_starts(accession, "WATDE"),
         growth_habit %in% keep_habits,
         own_membership >= min_membership,
         !is.na(takeall_score)) %>%
  group_by(group) %>%
  arrange(takeall_score, .by_group = TRUE) %>%
  mutate(rank_low = row_number(), rank_high = rev(row_number())) %>%
  filter(rank_low <= n_per_ag | rank_high <= n_per_ag) %>%
  mutate(end = if_else(rank_low <= n_per_ag, "lowest (protective)", "highest (susceptible)")) %>%
  ungroup() %>%
  select(group, end, accession, country, takeall_score, takeall_loci_protective,
         takeall_watkins_only_protective, paragon_ril, n_nil_families)
print(takeall_extremes, n = Inf, width = Inf)
write_csv(takeall_extremes, "Outputs/takeall_lowest_highest_per_AG.csv")

# 10e. Plot: every line's score by group, with shortlisted candidates labelled
p_line_scores <- ggplot(filter(master_takeall, !is.na(takeall_score)),
                        aes(x = group, y = takeall_score)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, colour = "grey40") +
  geom_jitter(aes(colour = group), width = 0.2, alpha = 0.4, size = 1) +
  geom_point(data = filter(master_takeall, accession %in% c(shortlist$accession,
                                                            "Paragon", "Hereward", "Soissons")),
             shape = 21, fill = "black", colour = "white", size = 2.5) +
  geom_text(data = filter(master_takeall, accession %in% c(shortlist$accession,
                                                           "Paragon", "Hereward", "Soissons")),
            aes(label = accession), size = 2.3, hjust = -0.15, check_overlap = TRUE) +
  scale_colour_manual(values = ag_colours, guide = "none") +
  labs(x = NULL, y = "Take-all haplotype score (% of trial mean)",
       title = "Predicted take-all from each line's haplotypes at 58 GWAS loci",
       subtitle = paste0("Below 0 = haplotypes linked to less take-all; black = shortlist & elite lines. ",
                         "Kruskal-Wallis across groups: p = ", format.pval(kw_score$p.value, digits = 2),
                         ", epsilon-squared = ", round(epsilon2, 2))) +
  theme_bw()
p_line_scores
ggsave("Outputs/07_takeall_line_scores.png", p_line_scores, width = 9, height = 5.5, dpi = 300)

# 10f. Heatmap: each shortlisted candidate's haplotype effect at each locus
profile_lines <- c(shortlist$accession, "Paragon", "Hereward", "Soissons")
p_profiles <- line_locus %>%
  inner_join(select(master, line_key, accession, group), by = "line_key") %>%
  filter(accession %in% profile_lines) %>%
  mutate(label = paste0(accession, " (", group, ")"),
         label = factor(label, levels = rev(unique(label[order(group, accession)]))),  # AG1 at top, elites at bottom
         locus = factor(locus, levels = locus_order)) %>%
  ggplot(aes(x = locus, y = label, fill = effect_pct)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-60, 60), oob = scales::squish,
                       name = "Take-all effect\n(% of trial mean)") +
  labs(x = "Take-all GWAS locus", y = NULL,
       title = "Take-all haplotype profile of shortlisted candidates and elite lines") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p_profiles
ggsave("Outputs/08_takeall_candidate_profiles.png", p_profiles, width = 13, height = 6, dpi = 300)


# ---- 11. Is the take-all score confounded by flowering time / height? -------
# Worry: if AG1/AG3 lines simply flower earlier (or are shorter), their take-all
# haplotypes might really be flowering-time haplotypes. Test it with heading
# date and plant height from the JIC field trials (Earlham Open Data):
#   heading date: 2006, 2010, 2014, 2020 trials;  height: 2006, 2020 trials
# Each trial counts days from a different reference date, so values are
# converted to z-scores within each trial (0 = trial average, +1 = one
# standard deviation later/taller) and then averaged per landrace.

jic_file <- "Raw_data/Watkins_JIC_phenotype_data.xlsx"
if (!file.exists(jic_file)) {
  download.file(paste0("https://opendata.earlham.ac.uk/wheat/under_license/toronto/",
                       "WatSeq_2023-09-15_landrace_modern_Variation_Data/Watseq_phenotype_data/",
                       "Natural_Populations/Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx"),
                jic_file, mode = "wb", quiet = TRUE)
}

# Read one trait from one trial sheet, average sub-lines, convert to z-scores
read_trial <- function(sheet, column, trial) {
  read_excel(jic_file, sheet = sheet) %>%
    transmute(accession = StoreCode,
              value     = suppressWarnings(as.numeric(.data[[column]]))) %>%   # "*" = missing
    group_by(accession) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(value)) %>%
    mutate(z = (value - mean(value)) / sd(value), trial = trial)
}
max_abs_z <- 5   # values more than 5 SD from the trial average are treated as errors

heading_trials <- bind_rows(
  read_trial("WGIN_Watkins_JIC_CFLN06", "Hd_dto_days-CFLN06",        "2006"),
  read_trial("WISP_Watkins_JIC_CFLN10", "Hd_dto_day-CFLN10",         "2010"),
  read_trial("DFW_Watkins_JI_CFLN14",   "Hd_dto_days-CFLN14",        "2014"),
  read_trial("DFW_Watkins_JI_CFLN20",   "Hd_dto_days-CFLN20-CFLN20", "2020"))

height_trials <- bind_rows(
  read_trial("WGIN_Watkins_JIC_CFLN06", "PH_M_cm-CFLN06",      "2006"),
  read_trial("DFW_Watkins_JI_CFLN20",   "PH_M_cm.GS91-CFLN20", "2020"))

# Remove implausible values (e.g. WATDE0714: 130 days in 2010 vs a median of 38)
filter(heading_trials, abs(z) > max_abs_z)    # prints what is removed
filter(height_trials,  abs(z) > max_abs_z)
heading_trials <- filter(heading_trials, abs(z) <= max_abs_z)
height_trials  <- filter(height_trials,  abs(z) <= max_abs_z)

# Do the heading-date trials agree with each other? (they should, r ~0.5-0.7)
heading_trials %>%
  select(accession, trial, z) %>%
  pivot_wider(names_from = trial, values_from = z) %>%
  select(-accession) %>%
  cor(use = "pairwise.complete.obs") %>% round(2)

phenology <- full_join(
  heading_trials %>% group_by(accession) %>% summarise(heading_z = mean(z), .groups = "drop"),
  height_trials  %>% group_by(accession) %>% summarise(height_z  = mean(z), .groups = "drop"),
  by = "accession")

pheno_takeall <- master_takeall %>%
  filter(str_starts(accession, "WATDE"), !is.na(takeall_score)) %>%
  inner_join(phenology, by = "accession") %>%
  mutate(growth_habit = coalesce(growth_habit, "unknown"))
nrow(pheno_takeall)

# 11a. Correlations: overall, and within each AG
# Pearson correlation with p-value and 95% confidence interval
cor_row <- function(d, label) {
  h <- cor.test(d$takeall_score, d$heading_z)
  t <- cor.test(d$takeall_score, d$height_z)
  tibble(group = label, n = nrow(d),
         r_heading = round(unname(h$estimate), 2),
         ci_heading = paste0(round(h$conf.int[1], 2), " to ", round(h$conf.int[2], 2)),
         p_heading = h$p.value,
         r_height = round(unname(t$estimate), 2),
         p_height = t$p.value)
}
confound_cor <- bind_rows(
  cor_row(pheno_takeall, "All landraces"),
  pheno_takeall %>% group_split(group) %>%
    lapply(function(d) cor_row(d, d$group[1])) %>% bind_rows() %>%
    mutate(q_heading = p.adjust(p_heading, "BH"),      # BH across the 7 AGs
           q_height  = p.adjust(p_height,  "BH"))) %>%
  mutate(sig_heading = stars(coalesce(q_heading, p_heading)),
         sig_height  = stars(coalesce(q_height,  p_height)))
confound_cor
write_csv(confound_cor, "Outputs/takeall_vs_phenology_correlations.csv")

# How early/late and tall/short is each AG on average?
pheno_takeall %>%
  group_by(group) %>%
  summarise(across(c(heading_z, height_z, takeall_score), ~ round(mean(.x, na.rm = TRUE), 2)),
            .groups = "drop")

# 11b. Do the AG differences survive adjusting for phenology?
#      Compare the AG effects (relative to AG1) with and without covariates.
m_group      <- lm(takeall_score ~ group, data = pheno_takeall)
m_covariates <- lm(takeall_score ~ heading_z + height_z + growth_habit, data = pheno_takeall)
m_both       <- lm(takeall_score ~ group + heading_z + height_z + growth_habit, data = pheno_takeall)

c(R2_group_only      = summary(m_group)$r.squared,
  R2_covariates_only = summary(m_covariates)$r.squared,
  R2_both            = summary(m_both)$r.squared) %>% round(3)
summary(m_both)
anova(m_covariates, m_both)    # does AG still matter once phenology is accounted for?

ag_effects <- tibble(term    = names(coef(m_group))[-1],
                     without_covariates = round(coef(m_group)[-1], 1)) %>%
  left_join(tibble(term = names(coef(m_both)),
                   with_covariates = round(coef(m_both), 1)), by = "term")
ag_effects
write_csv(ag_effects, "Outputs/takeall_AG_effects_adjusted.csv")

# 11c. Phenology-adjusted score for each line: remove the part of the score
#      predicted by heading date and height (slopes estimated WITHIN groups,
#      from m_both), so AG differences are kept but flowering time is not.
b_heading <- coef(m_both)[["heading_z"]]
b_height  <- coef(m_both)[["height_z"]]

pheno_takeall <- pheno_takeall %>%
  mutate(takeall_score_adj = round(takeall_score
                                   - b_heading * coalesce(heading_z, 0)
                                   - b_height  * coalesce(height_z, 0), 1))

shortlist_pheno <- pheno_takeall %>%
  filter(accession %in% c(shortlist$accession)) %>%
  select(group, accession, heading_z, height_z, takeall_score, takeall_score_adj) %>%
  mutate(across(c(heading_z, height_z), ~ round(.x, 2))) %>%
  arrange(group, takeall_score_adj)
print(shortlist_pheno, n = Inf)
write_csv(shortlist_pheno, "Outputs/candidate_takeall_phenology_adjusted.csv")
write_csv(pheno_takeall, "Outputs/takeall_line_scores_with_phenology.csv")

# 11d. Plot: take-all score vs heading date, with a trend line per AG
p_confound <- ggplot(pheno_takeall, aes(x = heading_z, y = takeall_score, colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_point(data = filter(pheno_takeall, accession %in% shortlist$accession),
             shape = 21, fill = "black", colour = "white", size = 2.5) +
  geom_text(data = filter(pheno_takeall, accession %in% shortlist$accession),
            aes(label = accession), size = 2.2, colour = "black",
            hjust = -0.15, check_overlap = TRUE) +
  scale_colour_manual(values = ag_colours,
                      labels = function(g) {
                        cc <- confound_cor[match(g, confound_cor$group), ]
                        paste0(g, "  r = ", cc$r_heading, " ", cc$sig_heading)
                      }) +
  labs(x = "Heading date (z-score; positive = later than average)",
       y = "Take-all haplotype score (% of trial mean)",
       colour = "AG (r within AG,\nBH-corrected)",
       title = "Is the take-all score linked to flowering time?",
       subtitle = with(filter(confound_cor, group == "All landraces"),
                       paste0("All landraces: r = ", r_heading, " (95% CI ", ci_heading, "), ",
                              fmt_p(p_heading), ". Lines = trend within each AG; black = candidates"))) +
  theme_bw()
p_confound
ggsave("Outputs/09_takeall_vs_heading.png", p_confound, width = 9, height = 6, dpi = 300)


# Save the R session info so you know which package versions made these results
writeLines(capture.output(sessionInfo()), "Outputs/sessionInfo.txt")
