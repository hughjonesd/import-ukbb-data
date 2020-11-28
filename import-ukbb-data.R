
suppressPackageStartupMessages({
  requireNamespace("tibble")
  requireNamespace("tidyr")
  library(drake)
  requireNamespace("matrixStats")
  library(raster)
  library(dplyr)
  library(readr)
  library(readxl)
  library(magrittr)
  library(purrr)
  library(forcats)
  library(santoku)
  library(tidync)
  library(abind)
})

data_dir       <- "../negative-selection-data"
pgs_dir        <- file.path(data_dir, "polygenic_scores")

data_dir       <- "../negative-selection-data"
pgs_dir        <- file.path(data_dir, "polygenic_scores")
sun_dir        <- file.path(data_dir, "sunshine-records")

pcs_file       <- file.path(data_dir, "UKB.HM3.100PCs.40310.txt")
pcs_40_file    <- file.path(data_dir, "ukb30545.40PCs.csv")
famhist_files  <- file.path(data_dir, c(
                    "UKB.EA_pheno.coordinates.QC.csv",
                    "david.family_history.traits.out.csv",
                    "david.family_history.traits.20042020.out.csv",
                    "david.family_history.traits.05052020.out.csv",
                    "david.family_history.traits.16052020.out.csv",
                    "david.family_history.traits.18052020.out.csv",
                    "david.family_history.traits.17062020.out.csv",
                    "david.birthinfo.traits.14072020.out.csv",
                    "david.traits.03112020.out.csv",
                    "david_traits.18112020.I.csv.out.csv",
                    "david_traits.18112020.II.csv.out.csv",
                    "david_traits.18112020.III.out.csv"
                  ))

rgs_file         <- file.path(data_dir, "EA3_rgs.10052019.rgs.csv")
relatedness_file <- file.path(data_dir, "relatedness_file.txt")
mf_pairs_file    <- file.path(data_dir, "spouse_pair_info", 
                                "UKB_out.mf_pairs_rebadged.csv")
alto_file <- file.path(data_dir, "spouse_pair_info", 
                         "UKB_out.all_living_together_rebadged.csv")
ashe_income_file <- file.path(data_dir, 
                      "SOC-income", "Occupation (4) Table 14.7a   Annual pay - Gross 2007.xls") 


import_famhist <- function (famhist_files, pgs_dir) {

  col_types <- cols(
            .default               = col_double(), 
            geno_measurement_plate = col_skip(),
            geno_measurement_well  = col_skip(),
            `53-0.0`               = col_date(),
            `53-1.0`               = col_date(),
            `53-2.0`               = col_date()
          )
  
  # this generates a "character" spec for all 22662-x.y columns
  fake_csv_string <- paste0("22662-0.", 0:18, collapse = ",")
  fake_csv_string <- paste0(fake_csv_string, "\n")
  col_types_22662 <- spec_csv(fake_csv_string, 
                                col_types = cols(
                                  .default = col_character()
                                )
                              )
  col_types$cols <- c(col_types$cols, col_types_22662$cols)
  
  fhl <-  lapply(famhist_files, read_csv, col_types = col_types)
  
  fhl[[1]] %<>% rename(f.eid = eid)
  fhl[-1] <- purrr::map(fhl[-1], ~ {
    names(.x) <- paste0("f.", names(.x))
    names(.x) <- gsub("\\-", ".", names(.x)) 
    .x
  })
  
  famhist <- purrr::reduce(fhl, left_join, by = "f.eid")

  # remove and rename some duplicated columns
  dupes <- grep("\\.x$", names(famhist), value = TRUE)
  for (dupe in dupes) {
    dupe_y <- sub("\\.x$", ".y", dupe)
    famhist[[dupe_y]] <- NULL
  }
  famhist %<>% rename_with(
                 .fn = ~ sub("\\.x$", "", .x), 
                 .cols = all_of(dupes)
               )
  
  # only self-identified, and genetically identified, white people
  famhist %<>% filter(f.21000.0.0 == 1001, ! is.na(genetic_ethnic_grouping))
  
  for (pgs_file in list.files(pgs_dir, pattern = "csv$", full.names = TRUE)) {
    score_name <- sub(".*UKB\\.AMC\\.(.*?)\\..*", "\\1", pgs_file, perl = TRUE)
    pgs <- read_delim(pgs_file, delim = " ", col_types = "dd")
    pgs %<>% filter(FID > 0) 
    names(pgs)[2] <- score_name # instead of "SCORE"
    famhist %<>% left_join(pgs, by = c("f.eid" = "FID"))
  }
  
  return(famhist)
}


import_score_names <- function (pgs_dir) {
  pgs_files <- list.files(pgs_dir, pattern = "csv$")
  score_names <-  sub(
                    ".*UKB\\.AMC\\.(.*?)\\..*", 
                    "\\1", 
                    pgs_files, 
                    perl = TRUE
                  )
  score_names <- setNames(score_names, score_names)

  score_names
}


clean_famhist <- function (famhist, score_names, sib_groups, ashe_income) {
  # we get very few extra cases from adding f.2946.1.0 etc, and it makes calculating
  # father's year of birth more complex
  
  # -10 means less than a year, we call that 0 full years
  famhist$f.699.0.0[famhist$f.699.0.0 == -10] <- 0
  famhist$f.699.1.0[famhist$f.699.1.0 == -10] <- 0
  famhist$f.699.2.0[famhist$f.699.2.0 == -10] <- 0
  
  # remove negatives
  famhist %<>% mutate(across(
      c(age_fulltime_edu, starts_with(c(
        "f.2946", "f.1845", "f.2754", "f.738",  "f.2764", "f.2405", "f.2734",
        "f.2149", "f.1873", "f.1883", "f.2784", "f.2794", "f.709", 
        "f.699", "f.3872", "f.728", "f.670", "f.680",
        "f.5057", "birth_lon", "birth_lat"
      ))), 
      negative_to_na
    )
  )
  # -7 means "never went to school"
  famhist$f.6138.0.0[famhist$f.6138.0.0 == -3] <- NA
  
  
  famhist$age_at_recruitment <- famhist$f.21022.0.0
  # questionnaire sex:
  famhist$female <- famhist$f.31.0.0 == 0 
  famhist$birth_year <- famhist$f.34.0.0
  famhist$birth_mon <- famhist$f.52.0.0
  
  # "Field 845 was collected from all participants except those who indicated 
  # they have a College or University degree, as defined by their answers to 
  # Field 6138". So, we impute this to be 21.
  famhist$age_fulltime_edu[is.na(famhist$age_fulltime_edu) & famhist$edu_qual.0.0 == 1] <- 21
  
  # TODO: ask Abdel how edu_qual was made up. Seems to be a "max" of all the
  # individual edu_qual.x.y's; these are just f.6138.x.y.
  # TODO: that doesn't work!!! Warn Abdel.
  # Currently using edu_qual.0.0
  famhist$university <- famhist$edu_qual.0.0 == 1
  famhist$income_cat <- famhist$f.738.0.0
  
  famhist$n_children <- pmax(famhist$f.2405.0.0, famhist$f.2405.1.0,
                             famhist$f.2405.2.0, famhist$f.2734.0.0, famhist$f.2734.1.0, 
                             famhist$f.2734.2.0, 
                             na.rm = TRUE
  )
  
  famhist$n_in_household <- famhist$f.709.0.0
  
  famhist$with_partner   <- famhist$f.6141.0.0 == 1
  # Many NAs, almost all from people living alone i.e. f.709 == 1
  famhist$with_partner[famhist$n_in_household == 1] <- FALSE
  famhist$with_partner[famhist$f.6141.0.0 == -3] <- NA
  

  famhist$age_fte_cat <- santoku::chop(famhist$age_fulltime_edu,
                                       c(16, 18),
                                       c("< 16", "16-18", "> 18"))

  # -7 means never went to school. We recode to 0 for simpliciy
  famhist$edu_qual[famhist$edu_qual == -7] <- 0
  famhist$edu_qual[famhist$edu_qual == -3] <- NA
  
  # we use pmax, assuming that people *can* have given birth for the first
  # time in between surveys.
  famhist$age_flb <-  pmax(
                        famhist$f.3872.0.0, famhist$f.3872.1.0, famhist$f.3872.2.0,
                        famhist$f.2754.0.0, famhist$f.2754.1.0, famhist$f.2754.2.0,
                        na.rm = TRUE
                      )
  famhist$age_flb_cat <- santoku::chop_equally(famhist$age_flb, 3, 
                                               labels = lbl_discrete("-"))
  
  # full brothers and sisters
  famhist$nbro <- pmax(famhist$f.1873.0.0, famhist$f.1873.1.0, 
                       famhist$f.1873.2.0, na.rm = TRUE)
  famhist$nsis <- pmax(famhist$f.1883.0.0, famhist$f.1883.1.0, 
                       famhist$f.1883.2.0, na.rm = TRUE)
  famhist$n_sibs <- famhist$nbro + famhist$nsis + 1
  # a few people give varying answers, we assume median is fine.
  # including later answers picks up c. 10K extra people.
  # some people have a non-integer median; we round down.
  famhist$birth_order <- floor(matrixStats::rowMedians(
    as.matrix(famhist[c("f.5057.0.0", "f.5057.1.0", "f.5057.2.0")]),
    na.rm = TRUE
  ))
  famhist$birth_order[famhist$n_sibs == 1] <- 0
  # number of older siblings, plus one:
  famhist$birth_order <- famhist$birth_order + 1
  # TODO: how does this come about??? Stupid answers?
  famhist$birth_order[famhist$birth_order > famhist$n_sibs] <- NA
  # TODO: why is birth_order often NaN?
  # TODO: why is f.5057 often NA when n_sibs == 1? And why often NA in general
  
  # TODO: minimum fath_age_birth is 3, moth_age_birth is 0... Why?
  famhist$fath_age <- famhist$f.2946.0.0
  famhist$moth_age <- famhist$f.1845.0.0
  famhist$fath_age_birth <- famhist$fath_age - famhist$age_at_recruitment
  famhist$moth_age_birth <- famhist$moth_age - famhist$age_at_recruitment

  # TODO: get f.20191, it is the touchscreen equivalent and will have more data
  famhist$fluid_iq <- rowMeans(famhist %>% select(starts_with("f.20016")), na.rm = TRUE)
  
  # first measurement has almost everyone
  # 1 = excellent, 2 = good, 3 = fair, 4 = poor
  # I reverse-code
  famhist$sr_health <- -1 * negative_to_na(famhist$f.2178.0.0)
  # "longstanding illness, disability or infirmity". 1 = TRUE
  famhist$illness   <- negative_to_na(famhist$f.2188.0.0)
  famhist$height <- famhist$f.50.0.0
  famhist$bmi <- famhist$f.21001.0.0
  
  famhist$num_jobs <- famhist$f.22599.0.0
  # job codes are f.22601.0.x
  # SOC2000 job codes are f.22617.0.x
  # start years   f.22602.0.x
  # end years     f.22603.0.x, which we don't have
  

  famhist[score_names] <- scale(famhist[score_names])

  # TODO: ask Abdel for job code f.20277
  # famhist %<>% left_join(ashe_income, by = c("f.20277" = "Code"))
  
  famhist %<>% 
        mutate(f.22617.0.0 = as.character(f.22617.0.0)) %>% 
        left_join(ashe_income, by = c("f.22617.0.0" = "Code")) %>% 
        select(-Description, -mean_pay) %>% 
        mutate(first_job_pay = median_pay/1000) %>% 
        select(-median_pay)
  
  famhist %<>% left_join(sib_groups, by = c("f.eid" = "id"))
  
  return(famhist)
}

