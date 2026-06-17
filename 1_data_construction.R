# Data construction for dementia algorithms
# Created by: Gina Nam
# Date: 11/23/2025
# Update log:
# 12/11/2025: perform a complete-case analysis

#---- Package loading ----

rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load(
  # Data manipulation 
  "tidyverse", "here", "haven", "stringr", "dplyr", "labelled",
  "gtsummary", "gt", "survey", "flextable", "officer", "purrr", "rlang",
  
  # Visualization
  "ggplot2", "gganimate", "gifski", "UpSetR"
)

#---- Load the datasets ----
source(here::here("0_paths_jl.R"))

## RAND longitudinal file ----
rand <- read_sas(here(
  path_to_hrs,
  "RANDLongitudinal",
  "RAND HRS Longitudinal File 2022",
  "randhrs1992_2022v1_SAS",
  "randhrs1992_2022v1.sas7bdat"
))

## Tracker file ----
trkr <- read_dta(here(
  path_to_hrs,
  "HRSRawData",
  "Tracker",
  "trk2022v2",
  "trk2022tr_r.dta"
))

## Apoe file ----
apoe_data <- read_sas(here(
  path_to_hrs,
  "HRSRawData",
  "Sensitive",
  "APOE",
  "apoe_serotonin_release.sas7bdat"
))

## Langa-Weir ----
langaweir <- read_sas(here(
  path_to_hrs,
  "HRSRawData",
  "Contributed Projects",
  "Cognition",
  "Langa-Weir Classification of Cognitive Function (1995-2020)",
  "cogfinalimp_9520wide.sas7bdat"
))

## Wu ----
memimp <- read_sas(here(
  path_to_hrs,
  "HRSCleanData",
  "ImputedMemoryScores",
  "dpmemimp_jan2024.sas7bdat"
))

## Expert Model dementia probability scores Gianatassio ----
dp <- read_sas(here(
  path_to_hrs,
  "HRSRawData",
  "Contributed Projects",
  "Cognition",
  "Gianattasio-Power Predicted Dementia Probability Scores and Dementia Classifications",
  "hrsdementia_2024_1217.sas7bdat"
))


## Hudomiet model ----
hudomiet <- read_dta(here(
  path_to_hrs,
  "HRSRawData",
  "Contributed Projects",
  "Cognition",
  "Predicted Cognition and Dementia Measures/PredictedCognitionDementiaMeasures",
  "Dementia_HRS_2000-2016_Basic_Release1_2m.dta"
))

#---- Data cleaning ----

## Convert all variable names to lowercase ----
colnames(apoe_data) <- tolower(colnames(apoe_data))
apoe_data <-apoe_data %>% 
  mutate(hhidpn = str_pad(paste(hhid, pn, sep=""),9,pad="0"))

names(rand) <- tolower(names(rand))
rand <- rand %>% 
  mutate(hhidpn = str_pad(hhidpn,9,pad="0"))

names(trkr) <- tolower(names(trkr))
trkr <- trkr %>% 
  mutate(hhidpn = str_pad(paste(hhid, pn, sep = ""), 9, pad = "0"))

names(langaweir) <- tolower(names(langaweir))
langaweir <- langaweir %>% 
  mutate(hhidpn = str_pad(paste(hhid, pn, sep = ""), 9, pad = "0"))

names(memimp) <- tolower(names(memimp))
memimp <- memimp %>% 
  mutate(hhidpn = str_pad(paste(hhid, pn, sep = ""), 9, pad = "0"))

names(dp) <- tolower(names(dp))
dp <- dp %>% 
  mutate(hhidpn = str_pad(paste(hhid, pn, sep = ""), 9, pad = "0"))

names(hudomiet) <- tolower(names(hudomiet))
hudomiet <- hudomiet %>%
  mutate(hhidpn = str_pad(hhidpn,9,pad="0"))

## Reshape and select variables ----
## For each algorithm, either extract or generate the cognitive status variable 
## for the year 2010.

### Langa-Weir reshape and create variables ----

langa_long <- langaweir %>%
  pivot_longer(
    cols = starts_with("cogfunction"),
    names_to = "year",
    values_to = "dem_3lw",
    names_pattern = "cogfunction(\\d{4})"   # pulls the number after w
  ) %>%
  mutate(year = as.integer(year),
         ## Langa–Weir: 1/2 = Normal/CIND (no dementia), 3 = Dementia
         dem_lw = case_when(
           dem_3lw %in% c(1, 2) ~ 0L,
           dem_3lw == 3         ~ 1L,
           TRUE                ~ NA_integer_
         )) %>%
  filter(year %in% c(2010, 2012, 2014, 2016, 2018, 2020)) %>%
  select(hhidpn, year, dem_3lw, dem_lw)

n_distinct(langa_long$hhidpn) #40130

# convert back to wide
langa_wide <- langa_long %>%
  pivot_wider(
    id_cols = hhidpn,
    names_from = year,
    values_from = c(dem_lw, dem_3lw),
    names_glue = "{.value}{year}"
  )

n_distinct(langa_wide$hhidpn) #40130

### Wu reshape and create variables ----
# For the Wu algorithm, we create "cog_memimpXX" using a cutoff of 0.5 to determine dementia
memimp_filter <- memimp %>%
  mutate(
    across(
      starts_with("dementpimp"),
      ~ if_else(. >= 0.5, 1L, 0L),
      .names = "cog_{.col}"
    )
  )

# check
table(memimp_filter$cog_dementpimp10,exclude=NULL)

# convert wide to long  
wu_long <- memimp_filter %>%
  pivot_longer(
    cols = starts_with("cog_dementpimp"),
    names_to = "year2",
    values_to = "dem_wu",
    names_pattern = "cog_dementpimp(\\d+)"
  ) %>%
  mutate(
    year2 = as.integer(year2),
    year = if_else(year2 >= 90, 1900L + year2, 2000L + year2)
  ) %>%
  filter(year %in% c(2010, 2012, 2014, 2016, 2018, 2020)) %>%
  select(hhidpn, year, dem_wu) 

n_distinct(wu_long$hhidpn) #38354

# convert back to wide
wu_wide <- wu_long %>%
  pivot_wider(
    id_cols     = hhidpn,
    names_from  = year,
    values_from = dem_wu,
    names_glue  = "dem_wu{year}"
  )

n_distinct(wu_wide$hhidpn) #38354

### Expert rename variables ----
expert_long <- dp %>%
  rename(
    dem_expert = expert_dem,  
    year       = hrs_year    
  ) %>%
  mutate(year = as.integer(year)) %>%
  select(hhidpn, year, dem_expert) %>%
  filter(year %in% c(2010, 2012, 2014, 2016, 2018, 2020))

# count unique ID
n_distinct(expert_long$hhidpn) #12903

# convert long to wide
expert_wide <- expert_long %>%
  pivot_wider(
    id_cols = hhidpn,
    names_from = year,
    values_from = dem_expert,
    names_prefix = "dem_expert"
  )

n_distinct(expert_wide$hhidpn) #12903

### Hudomiet create cog variables ----
wave_year_map <- tibble(
  wave = 5:13,
  year = c(2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014, 2016) 
)

hudomiet_long <- hudomiet %>%
  rename(hudomiet_cog = cog) %>%
  left_join(wave_year_map, by = "wave") %>%
  mutate(
    dem_3hudomiet = case_when(
      hudomiet_cog < 0       ~ 3,
      hudomiet_cog >= 0 & hudomiet_cog < 1 ~ 2,
      hudomiet_cog >= 1      ~ 1
  ),
  ## Hudomiet: 1/2 = Normal/CIND (no dementia), 3 = Dementia
    dem_hudomiet = case_when(
      dem_3hudomiet %in% c(1, 2) ~ 0L,
      dem_3hudomiet == 3         ~ 1L,
      TRUE                      ~ NA_integer_
    )) %>%
  select(hhidpn, year, dem_3hudomiet, dem_hudomiet, hudomiet_cog) %>%
  filter(year %in% c(2010, 2012, 2014, 2016))

n_distinct(hudomiet_long$hhidpn) #14643

# convert long to wide
hudomiet_wide <- hudomiet_long %>%
  pivot_wider(
    id_cols     = hhidpn,
    names_from  = year,
    values_from = c(dem_hudomiet, dem_3hudomiet),
    names_glue  = "{.value}{year}"
  )

n_distinct(hudomiet_wide$hhidpn) #14643

## Select variables from RAND and Tracker files ----
trkr_filter<-trkr%>%
  select(hhidpn,mage,yrenter,mwgtr,stratum,firstiw,secu) %>%
  mutate(mage = ifelse(mage == 999, NA, mage))

rand_filter<- rand %>%
  select(hhidpn,r10iwstat:r15iwstat,
         ragender,raracem,rahispan,rabplace,rafeduc,rameduc,raedyrs,r10mstat,r10shlt,raedegrm,
         r10diabe,r10stroke,r10smokev)

### convert rand to long since it has repeated measures included (iwstat)
### rand_long will be used for creating event variable
rand_long <- rand_filter %>%
  
  # pivot from wide to long
  pivot_longer(
    cols = r10iwstat:r15iwstat,
    names_to = "wave_raw",
    values_to = "iwstat"
  ) %>%
  
  # extract the wave number (10, 11, ..., 15)
  mutate(
    wave_num = as.integer(sub("r(\\d+)iwstat", "\\1", wave_raw)),
    # map wave number to calendar year
    year = 1992 + (wave_num - 1) * 2
  ) %>%
  
  # keep what you need 
  select(-wave_raw) %>%
  arrange(hhidpn, year)

### Combine RAND + Tracker files ----
covars_wide <- rand_filter %>%
  left_join(trkr_filter, by = "hhidpn")

n_distinct(covars_wide$hhidpn) #45234 (matched with RAND obs #)

# Create algorithm-specific wide datasets with generic demYYYY columns ----
## Function to ensure all dem-year column exist
dem_years <- c(2010, 2012, 2014, 2016, 2018, 2020)
dem_vars  <- paste0("dem", dem_years)

add_missing_dem_cols <- function(df, dem_vars) {
  missing <- setdiff(dem_vars, names(df))
  if (length(missing) > 0) {
    # add missing dem columns as NA_integer_
    df[missing] <- NA_integer_
  }
  df %>%
    select(hhidpn, algorithm, all_of(dem_vars), everything())
}

## Langa-weir analytic dataset ----
lw_wide_alg <- langa_wide %>%
  # keep only ID + the binary lw dementia variables
  select(hhidpn, matches("^dem_lw\\d{4}$|^dem_3lw\\d{4}$")) %>%
  # rename dem_lwYYYY -> demYYYY
  rename_with(
    ~ sub("^dem_lw(\\d{4})$", "dem\\1", .x),
    matches("^dem_lw\\d{4}$")
  ) %>%
  # merge covariates
  left_join(covars_wide, by = "hhidpn") %>%
  mutate(algorithm = "lw") %>%
  relocate(algorithm, .after = hhidpn) %>%
  add_missing_dem_cols(dem_vars = dem_vars)

n_distinct(lw_wide_alg$hhidpn) #40130

## Wu analytic dataset ----
wu_wide_alg <- wu_wide %>%
  # keep ID + Wu dementia indicators
  select(hhidpn, matches("^dem_wu\\d{4}$")) %>%
  # rename dem_wuYYYY -> demYYYY
  rename_with(
    ~ sub("^dem_wu(\\d{4})$", "dem\\1", .x),
    matches("^dem_wu\\d{4}$")
  ) %>%
  left_join(covars_wide, by = "hhidpn") %>%
  mutate(algorithm = "wu") %>%
  relocate(algorithm, .after = hhidpn) %>%
  add_missing_dem_cols(dem_vars = dem_vars)

## Expert analytic dataset ----
expert_wide_alg <- expert_wide %>%
  # ID + expert dementia indicators
  select(hhidpn, matches("^dem_expert\\d{4}$")) %>%
  # rename dem_expertYYYY -> demYYYY
  rename_with(
    ~ sub("^dem_expert(\\d{4})$", "dem\\1", .x),
    matches("^dem_expert\\d{4}$")
  ) %>%
  left_join(covars_wide, by = "hhidpn") %>%
  mutate(algorithm = "expert") %>%
  relocate(algorithm, .after = hhidpn) %>%
  add_missing_dem_cols(dem_vars = dem_vars)

## Hudomiet analytic dataset ----
hudomiet_wide_alg <- hudomiet_wide %>%
  # ID + Hudomiet binary dementia
  select(hhidpn, matches("^dem_hudomiet\\d{4}$|^dem_3hudomiet\\d{4}$"))%>%
  # rename dem_hudomiet2YYYY -> demYYYY
  rename_with(
    ~ sub("^dem_hudomiet(\\d{4})$", "dem\\1", .x),
    matches("^dem_hudomiet\\d{4}$")
  ) %>%
  left_join(covars_wide, by = "hhidpn") %>%
  mutate(algorithm = "hudomiet") %>%
  relocate(algorithm, .after = hhidpn) %>%
  add_missing_dem_cols(dem_vars = dem_vars)

# Merge ----
# stack all 4 algorithm datasets
analysis_ready <- bind_rows(
  lw_wide_alg,
  wu_wide_alg,
  expert_wide_alg,
  hudomiet_wide_alg
) %>%
  mutate(
    algorithm = factor(
      algorithm,
      levels = c("lw", "wu", "expert", "hudomiet")
    )
  )

table(analysis_ready$algorithm, useNA = "ifany") #good! each algorithm matches with observation # in each dataset

n_distinct(analysis_ready$hhidpn) #43371

table(analysis_ready$algorithm, analysis_ready$dem2010, useNA = "ifany")

## Recode and create variables ----
### wave enter ----
analysis_ready <- analysis_ready %>%
  mutate(wave_enter = case_when(
    yrenter == 1992 ~ "1992",
    yrenter %in% c(1993, 1994) ~ "1993–94",
    yrenter %in% c(1995, 1996) ~ "1995–96",
    yrenter == 1998 ~ "1998",
    yrenter == 2000 ~ "2000",
    yrenter == 2002 ~ "2002",
    yrenter == 2004 ~ "2004",
    yrenter == 2006 ~ "2006",
    yrenter == 2008 ~ "2008",
    yrenter == 2010 ~ "2010",
    yrenter == 2012 ~ "2012",
    yrenter == 2014 ~ "2014",
    yrenter == 2016 ~ "2016",
    yrenter == 2018 ~ "2018",
    yrenter == 2020 ~ "2020",
    yrenter == 2022 ~ "2022"
  )) %>%
  mutate(wave_enter = factor(wave_enter,
                             levels = c("1992", "1993–94", "1995–96", "1998", "2000", "2002", "2004", "2006", "2008", "2010", "2012", "2014", "2016", "2018", "2020","2022")))

### education_degree ----
# 0 = <HS, 1 = >=HS
analysis_ready <- analysis_ready %>%
  mutate(
    education_degree = case_when(
      raedegrm == 0 ~ 0L,
      raedegrm %in% c(1, 2, 3, 4, 5, 6, 7, 8) ~ 1L,
      TRUE ~ NA_integer_    
    )
  )

### race_ethnicity ----
# 1: NH White, 2: NH Black, 3: Hispanic, 4: NH Other
analysis_ready <- analysis_ready %>%
  mutate(
    race_ethnicity = case_when(
      rahispan == 1 ~ 3L,  # Hispanic
      rahispan == 0 & raracem == 1 ~ 1L,  # NH White
      rahispan == 0 & raracem == 2 ~ 2L,  # NH Black
      rahispan == 0 & raracem == 3 ~ 4L,  # NH Other
      TRUE ~ NA_integer_
    )
  )

### bornsouth ----
# 0 = No, 1 = Yes
analysis_ready <- analysis_ready %>%
  mutate(
    bornsouth = case_when(
      rabplace %in% c(5, 6, 7) ~ 1L,
      rabplace %in% c(1, 2, 3, 4, 8, 9, 10, 11, 12, 13) ~ 0L,
      TRUE ~ NA_integer_               
    )
  )

### highest parental education ----
analysis_ready <- analysis_ready %>%
  mutate(
    highestpedu = pmax(rafeduc, rameduc, na.rm = TRUE),
    pedu_cat = case_when(
      highestpedu < 8                      ~ 1L,  # <8 years
      highestpedu >= 8  & highestpedu <12 ~ 2L,  # 8–<12
      highestpedu == 12                    ~ 3L,  # 12
      highestpedu > 12                     ~ 4L,  # >12
      TRUE                                 ~ NA_integer_
    )
  )

### marital status ----
# 1 = Married/partnered, 0 = Not married/partnered
analysis_ready <- analysis_ready %>%
  mutate(
    mstat = case_when(
      r10mstat %in% c(1, 2, 3) ~ 1L,   
      r10mstat %in% c(4, 5, 6, 7, 8) ~ 0L,
      TRUE ~ NA_integer_               
    )
  )

### self-rated health ----
# 0 = Excellent/very good/good, 1 = Fair/poor
analysis_ready <- analysis_ready %>%
  mutate(
    shlt = case_when(
      r10shlt %in% c(1, 2, 3) ~ 0L,
      r10shlt %in% c(4, 5)    ~ 1L,
      TRUE                    ~ NA_integer_
    )
  )

### centered baseline age ----
analysis_ready <- analysis_ready %>%
  mutate(centered_baseline_age = mage - 70)

## Arrange ----
analysis_ready <- analysis_ready %>%
  arrange(hhidpn, algorithm)

# Event setup ----
## Dementia incidence ----
dem_years       <- c(2010, 2012, 2014, 2016, 2018, 2020)
followup_years  <- dem_years[dem_years > 2010]   # 2012–2020

analysis_ready <- analysis_ready %>%
  rowwise() %>%   # work within each person × algorithm row
  mutate(
    # 1. Baseline dementia at 2010 for this algorithm
    dem_base = as.integer(dem2010 == 1L),
    
    # 2. First dementia year after 2010 (or NA)
    first_dem_year = {
      # if baseline missing or demented at baseline → no incident
      if (is.na(dem2010) || dem2010 == 1L) {
        NA_real_
      } else {
        # collect post-2010 dementia values into a vector
        vals <- c_across(all_of(paste0("dem", followup_years)))
        has_dem <- any(vals == 1L, na.rm = TRUE)
        
        if (!has_dem) {
          NA_real_
        } else {
          # minimum follow-up year where dem == 1
          min(followup_years[vals == 1L], na.rm = TRUE)
        }
      }
    },
    
    # 3. Incident dementia indicator (ever gets dementia after 2010)
    dem_incident = as.integer(!is.na(first_dem_year)),
    
    # 4. Interval year of dementia onset (same idea as year_dem_* = first_dem_year - 1)
    year_dem = ifelse(!is.na(first_dem_year), first_dem_year - 1, NA_real_)
  ) %>%
  ungroup()

## Death ----
death_per_wave <- rand_long %>%    # long: hhidpn, year, iwstat
  group_by(hhidpn) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    first_death_wave = {
      cand <- year[year > 2010 & iwstat == 5L]
      if (length(cand) == 0) NA_real_ else min(cand)
    },
    .groups = "drop"
  ) %>%
  mutate(
    # as before: death assumed to occur in the interval before the death interview
    year_death = ifelse(is.na(first_death_wave),
                        NA_real_,
                        first_death_wave - 2)
  )

# join to stacked wide dataset
analysis_ready <- analysis_ready %>%
  left_join(death_per_wave, by = "hhidpn")


## dem score missing censoring
dem_year_cols <- paste0("dem", c(2010, 2012, 2014, 2016, 2018, 2020))

analysis_ready <- analysis_ready %>%
  rowwise() %>%
  mutate(
    year_dem_score_miss = {
      vals  <- c_across(all_of(dem_year_cols))
      years <- c(2010, 2012, 2014, 2016, 2018, 2020)
  
      if (algorithm == "hudomiet") {
        # Hudomiet:2016
        fu_idx <- which(years > 2010 & years <= 2016)
      } else {
        fu_idx <- which(years > 2010 & years <= 2020)
      }
      
      first_na_idx <- fu_idx[which(is.na(vals[fu_idx]))[1]]
      
      if (is.na(first_na_idx)) {
        NA_real_     # no missing score after baseline
      } else {
        # last non-NA wave BEFORE the first NA
        prior_idx <- which(!is.na(vals[1:(first_na_idx - 1)]))
        if (length(prior_idx) == 0) NA_real_
        else years[max(prior_idx)]
      }
    }
  ) %>%
  ungroup()
# If dem2010 is NA, these individuals will be excluded later during the eligibility step


## LTFU ----
gap_base <- rand_long %>%
  mutate(
    is_gap  = iwstat %in% c(4L, 7L),
    is_resp = iwstat == 1L
  ) %>%
  group_by(hhidpn) %>%
  arrange(year, .by_group = TRUE) %>%
  summarise(
    # first gap after baseline
    first_gap_year = {
      cand <- year[is_gap & year > 2010]
      if (length(cand) == 0) NA_real_ else min(cand)
    },
    # last response before that gap (including 2010)
    last_resp_before_gap = {
      if (is.na(first_gap_year)) {
        NA_real_
      } else {
        resp <- year[is_resp & year < first_gap_year]
        if (length(resp) == 0) NA_real_ else max(resp)
      }
    },
    year_ltfu = last_resp_before_gap,   # same definition as before
    .groups = "drop"
  )

# join to stacked wide dataset
analysis_ready <- analysis_ready %>%
  left_join(gap_base, by = "hhidpn")

## Final event variables per algorithm
analysis_ready <- analysis_ready %>%
  mutate(
    # admin end year depends on algorithm: Hudomiet ends in 2016
    admin_end_year = case_when(
      algorithm == "hudomiet"  ~ 2016,
      TRUE                     ~ 2020
    ),
    
    # treat missing years as Inf so they won't be the minimum
    T_dem   = ifelse(is.na(year_dem),   Inf, year_dem),
    T_death = ifelse(is.na(year_death), Inf, year_death),
    T_ltfu  = ifelse(is.na(year_ltfu),  Inf, year_ltfu),
    T_dem_score_miss = ifelse(is.na(year_dem_score_miss), Inf, year_dem_score_miss),
    T_admin = admin_end_year,
    
    # earliest event year
    end_year = pmin(T_dem, T_death, T_ltfu, T_dem_score_miss, T_admin, na.rm = TRUE),
    
    # event type for that algorithm
    event_type = case_when(
      end_year == T_dem   ~ "dementia",
      end_year == T_death ~ "death",
      end_year == T_ltfu  ~ "ltfu",
      end_year == T_dem_score_miss ~ "dem_score_miss",
      end_year == T_admin ~ "admin",
      TRUE                ~ NA_character_
    ),
    
    event_dementia = as.integer(event_type == "dementia"),
    event_death    = as.integer(event_type == "death"),
    event_ltfu     = as.integer(event_type == "ltfu"),
    event_dem_score_miss = as.integer(event_type == "dem_score_miss"),
    event_admin    = as.integer(event_type == "admin"),
    
    # analysis time in years since 2010
    time_from_2010 = end_year - 2010
  )

### Zero follow-up time flags (per algorithm row) ----

analysis_ready <- analysis_ready %>%
  mutate(
    drop_zero_time = time_from_2010 == 0 & event_type %in% c("death", "ltfu","dem_score_miss")
  )

# Counts of how many will be dropped per algorithm
drop_counts <- analysis_ready %>%
  group_by(algorithm) %>%
  summarise(
    n_drop = sum(drop_zero_time, na.rm = TRUE),
    .groups = "drop"
  )

drop_counts # lw 2369; wu 2099; expert 1250; hudomiet 1456

analysis_ready <- analysis_ready %>%
  filter(!drop_zero_time)


# join apoe dataset
apoe_filter <- apoe_data%>%
  select(hhidpn, apoe)

# recode apoe 1.34 and 44 together vs. none 2.34 and 44 together vs. none (dropping 24 as missing)
table(apoe_filter$apoe,exclude=NULL)

# 1. 24, 34, and 44 vs. all others
apoe_filter <- apoe_filter %>%
  mutate(
    apoe_24_34_44 = ifelse(apoe %in% c(24, 34, 44), 1, 0)
  )

table(apoe_filter$apoe_24_34_44, exclude = NULL)


# 2. 34 & 44 vs. all others, dropping 24 as missing
apoe_filter <- apoe_filter %>%
  mutate(
    apoe_34_44_drop24 = case_when(
      apoe %in% c(34, 44) ~ 1,
      apoe %in% c(22, 23, 33) ~ 0,
      apoe == 24 ~ NA
    )
  )

table(apoe_filter$apoe_34_44_drop24, exclude = NULL)

analysis_ready <-analysis_ready%>%
  left_join(apoe_filter,by="hhidpn")

# Eligibility indicator ----
baseline_year <- 2010

analysis_ready <- analysis_ready %>%
  mutate(
    # Row-level (algorithm-specific) eligibility for Table 1:
    # applies all criteria EXCEPT dementia-free at 2010
    elig_table1 = case_when(
      is.na(mage)                            ~ 0L,
      algorithm == "lw" & mage <50           ~ 0L,
      algorithm == "wu" & mage <70           ~ 0L,
      algorithm == "expert" & mage <70       ~ 0L,
      algorithm == "hudomiet" & mage <65     ~ 0L,
      is.na(race_ethnicity)                  ~ 0L,
      is.na(dem2010)                         ~ 0L,
      algorithm == "wu" & race_ethnicity == "3" ~ 0L,  # Wu-specific race restriction
      algorithm == "expert" & race_ethnicity == "4" ~ 0L,  # expert-specific race restriction
      TRUE                                   ~ 1L
    ),
    
    # Row-level eligibility for main analysis:
    # same as above, but ALSO require dementia-free at baseline (dem2010 == 0)
    elig_analysis = if_else(
      elig_table1 == 1L & dem2010 == 0L,
      1L, 0L
    )
  ) 


analysis_ready_age70 <- analysis_ready %>%
  mutate(
    # Table 1–style eligibility with uniform age >= 70
    # (no dementia-free requirement here)
    elig_table1_70 = case_when(
      is.na(mage)            ~ 0L,
      mage < 70              ~ 0L,                # uniform age restriction
      is.na(race_ethnicity)  ~ 0L,
      is.na(dem2010)         ~ 0L,
      algorithm == "wu" & race_ethnicity == "3" ~ 0L,  # still keep Wu race restriction
      algorithm == "expert" & race_ethnicity == "4" ~ 0L,  # expert-specific race restriction
      TRUE                   ~ 1L
    ),
    
    # dementia-free at baseline
    elig_analysis_70 = if_else(
      elig_table1_70 == 1L & dem2010 == 0L,
      1L, 0L
    )
  ) 


# Datasets ----
## age >=50 ----
### Table 1 dataset ----
dt_tab1_age50 <- analysis_ready %>%
  filter(elig_table1 == 1L)

### Analytic dataset ----
## Must be dementia-free at 2010 

dt_analysis_age50 <- analysis_ready %>%
  filter(
    elig_analysis == 1L   # person eligible somewhere
  )

n_distinct(dt_analysis_age50$hhidpn) #18204

## age >= 70 ----
### Table 1 dataset ----
dt_tab1_age70 <- analysis_ready_age70 %>%
  filter(elig_table1_70 == 1L)

### Analytic dataset ----
dt_analysis_age70 <- analysis_ready_age70 %>%
  filter(
    elig_analysis_70 == 1L
  )

# perform a complete-case analysis by removing individuals with missing values in variables that will be used in the model
model_vars <- c(
  "ragender", "mage", "education_degree", "race_ethnicity",
  "r10diabe", "r10stroke", "r10smokev"
)

dt_tab1_age50 <- dt_tab1_age50 %>%
  filter(complete.cases(across(all_of(model_vars))))

dt_analysis_age50 <- dt_analysis_age50 %>%
  filter(complete.cases(across(all_of(model_vars))))%>%
  select(-starts_with("dem_3lw"), -starts_with("dem_3hudomiet"))

dt_tab1_age70 <- dt_tab1_age70 %>%
  filter(complete.cases(across(all_of(model_vars))))

dt_analysis_age70_2cat <- dt_analysis_age70 %>%
  filter(complete.cases(across(all_of(model_vars))))%>%
  select(-starts_with("dem_3lw"), -starts_with("dem_3hudomiet"))

dt_analysis_age70_3cat <- dt_analysis_age70 %>%
  filter(complete.cases(across(all_of(model_vars))))

# Save merged files ----

save_path <- "D:/OneDrive - University of Southern California/R00/Manuscripts/HRS/Dementia algorithms/Data"

# saveRDS(analysis_ready,      file = file.path(save_path, "dt_before_eligibility.rds"))
# saveRDS(dt_tab1_age50,       file = file.path(save_path, "dt_tab1_age50.rds"))
# saveRDS(dt_analysis_age50,   file = file.path(save_path, "dt_analysis_age50.rds"))
# saveRDS(dt_tab1_age70,       file = file.path(save_path, "dt_tab1_age70.rds"))
# saveRDS(dt_analysis_age70_2cat,   file = file.path(save_path, "dt_analysis_age70.rds"))
# saveRDS(dt_analysis_age70_3cat,   file = file.path(save_path, "dt_analysis_age70_3cat.rds"))

# Sanity check ----
## Labels ----
recode_for_gtsum <- function(df) {
  df %>%
    mutate(
      algorithm = factor(
        algorithm,
        levels = c("lw", "wu", "expert", "hudomiet"),
        labels = c("Langa-Weir", "Wu", "Expert", "Hudomiet")
      ),
      dem2010 = factor(
        dem2010,
        levels = c(0, 1),
        labels = c("No dementia in 2010", "Dementia in 2010")
      ),
      event_type = factor(
        event_type,
        levels = c("dementia", "death", "ltfu", "admin"),
        labels = c("Dementia", "Death", "LTFU", "Administrative censoring")
      ),
      event_dementia = factor(
        event_dementia,
        levels = c(0, 1),
        labels = c("No incident dementia", "Incident dementia")
      )
    )
}

## Event type by dementia status in 2010, stratified by algorithm ----
### age >=50 ----

tbl_event_by_dem2010_main <- dt_tab1_age50 %>%
  recode_for_gtsum() %>%
  tbl_strata(
    strata = algorithm,
    .tbl = ~ .x %>%
      tbl_summary(
        by = dem2010,                    # columns = baseline dementia status
        include = event_type,            # row = event type
        type = event_type ~ "categorical",
        missing = "no"
      )
  ) %>%
  modify_caption(
    "**Event Type by Dementia Status in 2010, Stratified by Algorithm (Age ≥ 50 before dropping 2010 dementia)**"
  ) %>%
  modify_footnote(
    update = everything() ~
      "Rows: event types; columns: dementia status in 2010. One panel per algorithm."
  )

tbl_event_by_dem2010_main

### age >= 70 ----
tbl_event_by_dem2010_age70 <- dt_tab1_age70 %>%
  recode_for_gtsum() %>%
  tbl_strata(
    strata = algorithm,
    .tbl = ~ .x %>%
      tbl_summary(
        by = dem2010,
        include = event_type,
        type = event_type ~ "categorical",
        missing = "no"
      )
  ) %>%
  modify_caption(
    "**Event Type by Dementia Status in 2010, Stratified by Algorithm (Age ≥ 70 before dropping 2010 dementia)**"
  ) %>%
  modify_footnote(
    update = everything() ~
      "Rows: event types; columns: dementia status in 2010. One panel per algorithm; restricted to age ≥ 70 definition."
  )

tbl_event_by_dem2010_age70


## Dementia incidence by algorithm ----
### age >=50 ----

tbl_incidence_main <- dt_analysis_age50 %>%
  recode_for_gtsum() %>%
  tbl_summary(
    by = algorithm,                     # columns = algorithm
    include = event_dementia,           # row = dementia incidence
    type = event_dementia ~ "categorical",
    missing = "no"
  ) %>%
  modify_caption(
    "**Dementia Incidence by Algorithm (Age ≥ 50 Analytic Sample)**"
  ) %>%
  modify_footnote(
    update = everything() ~
      "Denominator: participants dementia-free in 2010 and eligible under main rules."
  )

tbl_incidence_main

### age >= 70 ----

tbl_incidence_age70 <- dt_analysis_age70 %>%
  recode_for_gtsum() %>%
  tbl_summary(
    by = algorithm,
    include = event_dementia,
    type = event_dementia ~ "categorical",
    missing = "no"
  ) %>%
  modify_caption(
    "**Dementia Incidence by Algorithm (Age ≥ 70 Analytic Sample)**"
  ) %>%
  modify_footnote(
    update = everything() ~
      "Denominator: participants dementia-free in 2010 and eligible under age ≥ 70 rules."
  )

tbl_incidence_age70

## Event type by algorithm ----
### age >= 50 ----
tbl_eventtype_main <- dt_analysis_age50 %>%
  recode_for_gtsum() %>%
  tbl_summary(
    by = algorithm,
    include = event_type,
    type = event_type ~ "categorical",
    missing = "no"
  ) %>%
  modify_caption(
    "**Event Type by Algorithm (Age ≥ 50 Analytic Sample)**"
  ) %>%
  modify_footnote(
    update = everything() ~
      "Event type is defined by the earliest of dementia, death, LTFU, or administrative censoring."
  )

tbl_eventtype_main

### age >= 70 ----
tbl_eventtype_age70 <- dt_analysis_age70 %>%
  recode_for_gtsum() %>%
  tbl_summary(
    by = algorithm,
    include = event_type,
    type = event_type ~ "categorical",
    missing = "no"
  ) %>%
  modify_caption(
    "**Event Type by Algorithm (Age ≥ 70 Analytic Sample)**"
  ) %>%
  modify_footnote(
    update = everything() ~
      "Restricted to age ≥ 70 analysis eligibility."
  )

tbl_eventtype_age70

## Save each result ----
# save_path_output <- "D:/OneDrive - University of Southern California/R00/Manuscripts/HRS/Dementia algorithms/Output"
# 
# saveRDS(tbl_event_by_dem2010_main, file = file.path(save_path_output, "tbl_event_by_dem2010_age50.rds"))
# saveRDS(tbl_event_by_dem2010_age70, file = file.path(save_path_output, "tbl_event_by_dem2010_age70.rds"))
# 
# saveRDS(tbl_incidence_main,        file = file.path(save_path_output, "tbl_incidence_age50.rds"))
# saveRDS(tbl_incidence_age70,        file = file.path(save_path_output, "tbl_incidence_age70.rds"))
# 
# saveRDS(tbl_eventtype_main,        file = file.path(save_path_output, "tbl_eventtype_age50.rds"))
# saveRDS(tbl_eventtype_age70,        file = file.path(save_path_output, "tbl_eventtype_age70.rds"))



