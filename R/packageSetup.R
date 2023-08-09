packageSetup <-
  function(){
    library(vroom)
    library(tidyverse)
    library(parallel)
    options(ggrepel.max.overlaps = Inf)
    
    # ideally would make these global variables?
    vocabulary.map <- vroom::vroom("../Data files/mg_phecode_map_5.5.2023.csv",
                                   .name = janitor::make_clean_names, delim = ",",
                                   col_types = c(vocabulary_id = "c", code = "c", phecode = "c"))
    
    rollup.map <- vroom::vroom("../Data files/mg_phecode_rollup_5.5.2023.csv",
                               .name = janitor::make_clean_names, delim = ",",
                               col_types = c(code = "c", phecode_unrolled = "c"))
    
    exclusion.map <- vroom::vroom("../Data files/mg_phecode_exclude_5.5.2023.csv",
                                  .name = janitor::make_clean_names, delim = ",",
                                  col_types = c(code = "c", exclusion_criteria = "c"))
    
    gender_restriction <- vroom::vroom("../Data files/mg_gender_restriction_5.5.2023.csv",
                                       .name = janitor::make_clean_names, delim = ",",
                                       col_types = c(phecode = "c", male_only = "logical", female_only = "logical"))
    
    annotate.phenotype.description <- vroom::vroom("../Data files/mg_pheinfo_5.5.2023.csv",
                                                   .name = janitor::make_clean_names, delim = ",",
                                                   col_types = c(phecode = "c", description = "c", groupnum = "double", group = "c", color = "c"))
    
    example_data_icd <- vroom::vroom("../Data files/person_icd.csv",
                                     .name = janitor::make_clean_names, delim = ",",
                                     col_types = c(id = "i", vocabulary_id = "c", code = "c", count = "i"))
    
    example_data_gender <- vroom::vroom("../Data files/test_gender.csv",
                                        .name = janitor::make_clean_names, delim = ",",
                                        col_types = c(id = "i", sex = "c"))
    
    example_data_genotypes <- vroom::vroom("../Data files/test_genotypes.csv",
                                           .name = janitor::make_clean_names, delim = ",",
                                           col_types = c(id = "i", rs_example = "i"))
    
  }