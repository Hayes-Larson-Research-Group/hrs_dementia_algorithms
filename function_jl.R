# Labeling for Table 1

recode_for_gtsum <- function(df, dataset_name = NULL) {
  df <- df %>%
    mutate(
      algorithm_f = factor(
        algorithm,
        levels = c("lw", "wu", "expert", "hudomiet"),
        labels = c("Langa-Weir", "Wu", "Expert", "Hudomiet")
      ),
      dem2010_f = factor(
        dem2010,
        levels = c(0, 1),
        labels = c("No dementia in 2010", "Dementia in 2010")
      ),
      event_type_f = factor(
        event_type,
        levels = c("dementia", "death", "ltfu", "admin"),
        labels = c("Dementia", "Death", "LTFU", "Administrative censoring")
      ),
      event_dementia_f = factor(
        event_dementia,
        levels = c(0, 1),
        labels = c("No incident dementia", "Incident dementia")
      ),
      ragender_f = factor(
        ragender,
        levels = c(1, 2),
        labels = c("Male", "Female")
      ),
      race_ethnicity_f = factor(
        race_ethnicity,
        levels = intersect(c(1, 2, 3, 4), unique(race_ethnicity)),
        labels = c("NH White", "NH Black", "Hispanic", "NH Other")[intersect(c(1, 2, 3, 4), unique(race_ethnicity))]
      ),
      bornsouth_f = factor(
        bornsouth,
        levels = c(0, 1),
        labels = c("No", "Yes")
      ),
      education_degree_f = factor(
        education_degree,
        levels = c(1, 0),
        labels = c("≥ HS","< HS")
      ),
      pedu_cat_f = factor(
        pedu_cat,
        levels = c(1, 2, 3, 4),
        labels = c("<8 years", "8-12 years", "12 years", ">12 years")
      ),
      shlt_f = factor(
        shlt,
        levels = c(0, 1),
        labels = c("Excellent/very good/good", "Fair/poor")
      ),
      mstat_f = factor(
        mstat,
        levels = c(1, 0),
        labels = c("Married/partnered", "Not married/partnered")
      ),
      r10smokev_f = factor(
        r10smokev,
        levels = c(0, 1),
        labels = c("No", "Yes")
      ),
      r10diabe_f = factor(
        r10diabe,
        levels = c(0, 1),
        labels = c("No", "Yes")
      ),
      r10stroke_f = factor(
        r10stroke,
        levels = c(0, 1),
        labels = c("No", "Yes")
      ),
      apoe_24_34_44_f = factor(
        apoe_24_34_44,
        levels = c(0, 1),
        labels = c("All others", "24/34/44")
      ),
      apoe_34_44_drop24_f = factor(
        apoe_34_44_drop24,
        levels = c(0, 1),
        labels = c("All others (excluding 24)", "34/44")
      ),
     #  dem_3lw2010_f = factor(
     #    dem_3lw2010,
     #    levels = c(1, 2, 3),
     #    labels = c("Normal", "CIND", "Dementia")
     #  ),
     # dem_3hudomiet2010_f = factor(
     #      dem_3hudomiet2010,
     #      levels = c(1, 2, 3),
     #      labels = c("Normal", "CIND", "Dementia")
     # )
    )
  
  # Only include dem_3* variables for Table 1 datasets
  if (!is.null(dataset_name) && dataset_name %in% c("dt_tab1_age50", "dt_tab1_age70")) {
    df <- df %>%
      mutate(
        dem_3lw2010_f = factor(
          dem_3lw2010,
          levels = c(1, 2, 3),
          labels = c("Normal", "CIND", "Dementia")
        ),
        dem_3hudomiet2010_f = factor(
          dem_3hudomiet2010,
          levels = c(1, 2, 3),
          labels = c("Normal", "CIND", "Dementia")
        )
      )
  }
  
  return(df)
}





## Example of how to use it. See line 32:

# tbl_event_by_dem2010_main <- dt_tab1_age50 %>%
#   recode_for_gtsum() %>%
#   tbl_strata(
#     strata = algorithm,
#     .tbl = ~ .x %>%
#       tbl_summary(
#         by = dem2010,                    # columns = baseline dementia status
#         include = event_type,            # row = event type
#         type = event_type ~ "categorical",
#         missing = "no"
#       )
#   ) ....