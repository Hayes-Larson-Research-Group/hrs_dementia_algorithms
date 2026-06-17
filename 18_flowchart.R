library(dplyr)
library(tidyr)
library(purrr)


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

# analysis_ready <- analysis_ready %>%
#   filter(!drop_zero_time)


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


id_var <- "hhidpn"

gn_check <- rand %>%
  left_join(trkr_filter, by = "hhidpn")

gn_check %>%
  filter(r10iwstat == "1", mage >= 70) %>%
  summarise(count = n())

gn_check <- gn_check %>%
  mutate(
    race_ethnicity = case_when(
      rahispan == 1 ~ 3L,  # Hispanic
      rahispan == 0 & raracem == 1 ~ 1L,  # NH White
      rahispan == 0 & raracem == 2 ~ 2L,  # NH Black
      rahispan == 0 & raracem == 3 ~ 4L,  # NH Other
      TRUE ~ NA_integer_
    )
  )

gn_check <- gn_check %>%
  filter(r10iwstat == "1", mage >= 70)

analysis_ready_age70 <- analysis_ready_age70 %>%
  filter(mage >= 70)

n_distinct(gn_check$hhidpn) #8523
n_distinct(analysis_ready_age70$hhidpn) #7323

# gn_check should already be filtered to r10iwstat == 1 and mage >= 70
gn_check_u <- gn_check %>%
  distinct(.data[[id_var]], .keep_all = TRUE)


# run_algo_flow <- function(analysis_ready_age70, gn_check, algo_code, algo_label = algo_code, id_var = "id") {
#   
#   gn_check_u <- gn_check %>%
#     distinct(.data[[id_var]], .keep_all = TRUE)%>%
#     select(-any_of("race_ethnicity"))
#   
#   algo_dat <- analysis_ready_age70 %>%
#     filter(algorithm == algo_code) %>%
#     distinct(.data[[id_var]], .keep_all = TRUE) %>%
#     select(all_of(id_var), algorithm, race_ethnicity, dem2010, drop_zero_time)
#   
#   merged <- algo_dat %>%
#     right_join(gn_check_u, by = c(id_var))
#   
#   # 1. total available in algo
#   step1 <- merged %>%
#     filter(!is.na(algorithm))
#   
#   # 2. keep those without zero time
#   # If drop_zero_time is logical, this is correct:
#   step2 <- step1 %>%
#     filter(!drop_zero_time)
#   
#   # 3. included if dem2010 == 0
#   step3 <- step2 %>%
#     filter(dem2010 == 0)
#   
#   # 4. no missing in race
#   step4 <- step3 %>%
#     filter(!is.na(race_ethnicity))
#   
#   # 5. what's left
#   final_data <- step4
#   
#   flow_table <- tibble(
#     algorithm = algo_label,
#     step = c(
#       paste0("1. total available in ", algo_label),
#       "2. keep those without zero time",
#       "3. included if dem2010 == 0",
#       "4. no missing in race",
#       "5. what's left"
#     ),
#     n = c(
#       n_distinct(step1[[id_var]]),
#       n_distinct(step2[[id_var]]),
#       n_distinct(step3[[id_var]]),
#       n_distinct(step4[[id_var]]),
#       n_distinct(final_data[[id_var]])
#     )
#   )
#   
#   list(
#     flow_table = flow_table,
#     final_data = final_data
#   )
# }



table(gn_check_u$race_ethnicity,exclude=NULL)


run_algo_flow <- function(analysis_ready_age70, gn_check, algo_code, algo_label = algo_code, id_var = "id") {
  
  gn_check_u <- gn_check %>%
    distinct(.data[[id_var]], .keep_all = TRUE) %>%
    filter(!is.na(race_ethnicity))
  
  algo_dat <- analysis_ready_age70 %>%
    filter(algorithm == algo_code) %>%
    distinct(.data[[id_var]], .keep_all = TRUE) %>%
    select(all_of(id_var), algorithm, dem2010, drop_zero_time)
  
  merged <- algo_dat %>%
    right_join(gn_check_u, by = id_var)
  
  # ---- Step 2: # of race/ethnicity group dropped ----
  if (algo_code == "expert") {
    n_race_dropped <- merged %>% filter(race_ethnicity == 4) %>% nrow()
    after_race_grp <- merged %>% filter(race_ethnicity != 4)
    race_drop_label <- "2. # NH Other dropped (Expert)"
  } else if (algo_code == "wu") {
    n_race_dropped <- merged %>% filter(race_ethnicity == 3) %>% nrow()
    after_race_grp <- merged %>% filter(race_ethnicity != 3)
    race_drop_label <- "2. # Hispanic dropped (Wu)"
  } else {
    n_race_dropped <- 0
    after_race_grp <- merged
    race_drop_label <- "2. # race/ethnicity group dropped (none)"
  }
  
  n_after_race_grp <- nrow(after_race_grp)
  
  # ---- Step 3: # missing dem2010 variable ----
  n_missing_dem2010 <- after_race_grp %>%
    filter(is.na(dem2010)) %>%
    nrow()
  
  step3_remain <- after_race_grp %>%
    filter(!is.na(dem2010))
  n_step3_remain <- nrow(step3_remain)
  
  # ---- Step 4: # of prevalent dementia in 2010 ----
  n_prev_dem <- step3_remain %>%
    filter(dem2010 == 1) %>%
    nrow()
  
  step4_remain <- step3_remain %>%
    filter(dem2010 == 0)
  n_step4_remain <- nrow(step4_remain)
  
  # ---- Step 5: # with zero time follow-up ----
  n_zero_time <- step4_remain %>%
    filter(drop_zero_time) %>%
    nrow()
  
  step5_remain <- step4_remain %>%
    filter(!drop_zero_time)
  n_step5_remain <- nrow(step5_remain)
  
  flow_table <- tibble(
    algorithm = algo_label,
    step = c(
      race_drop_label,
      "2a. Remaining after dropping race group",
      "3. # missing dem2010 variable",
      "3a. Remaining after dropping missing dem2010",
      "4. # prevalent dementia in 2010",
      "4a. Remaining after dropping prevalent dementia",
      "5. # with zero time follow-up",
      "5a. Remaining after dropping zero time FU"
    ),
    n = c(
      n_race_dropped,
      n_after_race_grp,
      n_missing_dem2010,
      n_step3_remain,
      n_prev_dem,
      n_step4_remain,
      n_zero_time,
      n_step5_remain
    )
  )
  
  list(
    flow_table = flow_table,
    final_data = step5_remain
  )
}


algo_map <- tibble::tibble(
  algo_code = c("lw", "wu", "expert", "hudomiet"),
  algo_label = c("LW", "Wu", "Expert", "Hudomiet")
)

results <- purrr::pmap(
  algo_map,
  function(algo_code, algo_label) {
    run_algo_flow(
      analysis_ready_age70 = analysis_ready_age70,
      gn_check = gn_check,
      algo_code = algo_code,
      algo_label = algo_label,
      id_var = "hhidpn"
    )
  }
)

combined_flow_table <- purrr::map_dfr(results, "flow_table")

final_datasets <- setNames(
  purrr::map(results, "final_data"),
  algo_map$algo_label
)

print(combined_flow_table)

# write.csv(
#   combined_flow_table,
#   file = "D:/OneDrive - University of Southern California/R00/Manuscripts/HRS/Dementia algorithms/Output/combined_flow_table.csv",
#   row.names = FALSE
# )

list2env(
  setNames(
    final_datasets,
    paste0("final_", names(final_datasets))
  ),
  envir = .GlobalEnv
)



#complete case
model_vars <- c(
  "ragender", "mage", "education_degree", "race_ethnicity",
  "r10diabe", "r10stroke", "r10smokev"
)

purrr::map(
  final_datasets,
  ~ setdiff(model_vars, names(.x))
)


final_datasets <- purrr::map(
  final_datasets,
  ~ .x %>%
    mutate(
      education_degree = case_when(
        raedegrm == 0 ~ 0L,
        raedegrm %in% c(1, 2, 3, 4, 5, 6, 7, 8) ~ 1L,
        TRUE ~ NA_integer_
      )
    )
)

final_datasets_cc <- purrr::map(
  final_datasets,
  ~ .x %>% filter(complete.cases(across(all_of(model_vars))))
)

tibble(
  algorithm = names(final_datasets),
  n_before = purrr::map_int(final_datasets, nrow),
  n_after  = purrr::map_int(final_datasets_cc, nrow)
)

purrr::map_lgl(final_datasets_cc, ~ "apoe_24_34_44" %in% names(.x))

# add apoe
final_datasets_apoe <- purrr::map(
  final_datasets_cc,
  ~ .x %>% left_join(apoe_filter, by = "hhidpn")
)

# remove apoe_24_34_44 na
final_datasets_apoe_cc <- purrr::map(
  final_datasets_apoe,
  ~ .x %>% filter(!is.na(apoe_24_34_44))
)

apoe_summary <- tibble(
  algorithm   = names(final_datasets_apoe),
  n_before    = purrr::map_int(final_datasets_apoe, nrow),
  n_after     = purrr::map_int(final_datasets_apoe_cc, nrow),
  n_dropped   = n_before - n_after
)
print(apoe_summary)






