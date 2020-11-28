

data_dir       <- "../negative-selection-data"

pcs_file       <- file.path(data_dir, "UKB.HM3.100PCs.40310.txt")

pcs_40_file    <- file.path(data_dir, "ukb30545.40PCs.csv")

famhist_files <- file.path(data_dir, c(
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
