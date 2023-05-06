### Data Preprocessing
The flowchart below illustrates in detail how the data is consumed and processed by various R scripts at each stage. The coloured silos, rounded boxes, and grey boxes, represent data, R scripts, and directories respectively. The data is coloured according to the [medallion architecture](https://www.databricks.com/glossary/medallion-architecture), where bronze represents raw data, silver represents cleansed and processed data, and gold represents data which are consumption ready (i.e. directly used within statistical models).

#### Preprocessing Stage 1
```mermaid
  flowchart LR
    subgraph data
      d1[(aggregated_participants_count_art_coverage.csv)]
      d2[(aggregated_newlyregistered_count_art_coverage.csv)]
      d3[(aggregated_participants_count_unsuppressed.csv)]
      d4[(aggregated_newlyregistered_count_unsuppressed.csv)]
      d5[(sensitivity_specificity_art.csv)]
    end
    style d1 fill:#CD7F32
    style d2 fill:#CD7F32
    style d3 fill:#CD7F32
    style d4 fill:#CD7F32
    style d5 fill:#CD7F32
    style data fill:#F5F5F5

    subgraph misc
      r1(get_estimate_art_coverage_participants.R)
      r2(get_estimate_art_coverage_non_participants.R)
      r3(get_estimates_unsuppressed_proportion_participants.R)
      r4(get_estimates_unsuppressed_proportion_non_participants.R)
      r5(get_treatment_cascade_participants.R)
      r6(get_treatment_cascade_non_participants.R)
    end
    style misc fill:#F5F5F5

    subgraph fit
      f1[(RCCS_art_posterior_samples_221208.rds)]
      f2[(RCCS_art_posterior_samples_newlyregistered_221208.rds)] 
      f3[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds)]
      f4[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds)]

      f5[(RCCS_treatment_cascade_participants_posterior_samples_221208.rds)]
      f6[(RCCS_treatment_cascade_participants_estimates_221208.csv)]

      f7[(RCCS_treatment_cascade_nonparticipants_posterior_samples_221208.rds)]
      f8[(RCCS_treatment_cascade_nonparticipants_estimates_221208.csv)]      
    end
    style f1 fill:#C0C0C0
    style f2 fill:#C0C0C0
    style f3 fill:#C0C0C0
    style f4 fill:#C0C0C0

    style f5 fill:#FFD700
    style f6 fill:#FFD700
    style f7 fill:#FFD700
    style f8 fill:#FFD700
    style fit fill:#F5F5F5

    d1 -- load --> r1
    r1 -- create --> f1

    d2 -- load --> r2
    r2 -- create --> f2

    d3 -- load --> r3
    r3 -- create --> f3

    d4 -- load --> r4
    r4 -- create --> f4

    d5 -- load --> r5
    f1 -- load --> r5
    f3 -- load --> r5
    r5 -- create --> f5
    r5 -- create --> f6

    d5 -- load --> r6
    f2 -- load --> r6
    f4 -- load --> r6
    r6 -- create --> f7
    r6 -- create --> f8
```

#### Preprocessing Stage 2
```mermaid
  flowchart LR
    subgraph data
      d1[(RCCS_participation_221208.csv)]
      d2[(sensitivity_specificity_art.csv)]
      d3[(aggregated_count_hiv_positive.csv)]
      d4[(RCCS_census_eligible_individuals_221116.csv)]
    end
    style d1 fill:#CD7F32
    style d2 fill:#CD7F32
    style d3 fill:#CD7F32
    style d4 fill:#CD7F32
    style data fill:#F5F5F5

    subgraph misc
      r1(get_treatment_cascade_population.R)
      r2(get_estimates_prevalence.R)
      
      r3(get_unsuppressed_median_age.R)
      r4(get_unsuppressed_share_sex.R)
      r5(get_unsuppressed_ratio_sex.R)
    end
    style misc fill:#F5F5F5

    subgraph fit
      f1[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds)]
      f2[(RCCS_art_posterior_samples_221208.rds)]
      f3[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds)]
      f4[(RCCS_art_posterior_samples_newlyregistered_221208.rds)]

      f5[(RCCS_treatment_cascade_population_posterior_samples_221208.rds)]
      f6[(RCCS_treatment_cascade_population_estimates_221208.rds)]

      f7[(RCCS_prevalence_posterior_sample_221116.rds)]
      f8[(RCCS_prevalence_estimates_221116.csv)]
      
      f9[(RCCS_unsuppressed_median_age_221208.csv)]
      f10[(RCCS_unsuppressed_share_sex_221208.csv)]
      f11[(RCCS_unsuppressed_ratio_sex_221208.csv)]
    end
    style f1 fill:#C0C0C0
    style f2 fill:#C0C0C0
    style f3 fill:#C0C0C0
    style f4 fill:#C0C0C0

    style f5 fill:#C0C0C0
    style f6 fill:#FFD700

    style f7 fill:#C0C0C0
    style f8 fill:#FFD700

    style f9 fill:#FFD700
    style f10 fill:#FFD700
    style f11 fill:#FFD700
    
    style fit fill:#F5F5F5

    d1 -- load --> r1
    d2 -- load --> r1
    f1 -- load --> r1
    f2 -- load --> r1
    f3 -- load --> r1
    f4 -- load --> r1
    r1 -- create --> f5
    r1 -- create ----> f6

    d3 -- load --> r2
    r2 -- create --> f7
    r2 -- create ----> f8

    d4 -- load --> r3
    f5 -- load --> r3
    f7 -- load --> r3
    r3 -- create --> f9

    d4 -- load --> r4
    f5 -- load --> r4
    f7 -- load --> r4
    r4 -- create --> f10

    d4 -- load --> r5
    f5 -- load --> r5
    f7 -- load --> r5
    r5 -- create --> f11
```

#### Preprocessing (Sensitivity Analysis)
```mermaid
  flowchart LR
    subgraph data
      d1[(aggregated_participants_count_art_coverage_vl200.csv)]
      d2[(aggregated_newlyregistered_count_unsuppressed_vl200.csv)]
      d3[(aggregated_newlyregistered_count_art_coverage_vl200.csv)]
      d4[(aggregated_participants_count_unsuppressed_vl200.csv)]
      d5[(sensitivity_specificity_art_vl200.csv)]
    end
    style d1 fill:#CD7F32
    style d2 fill:#CD7F32
    style d3 fill:#CD7F32
    style d4 fill:#CD7F32
    style d5 fill:#CD7F32
    style data fill:#F5F5F5

    subgraph misc
      r1(get_estimates_art_coverage_participants_vl200.R)
      r2(get_estimates_art_coverage_non_participants_vl200.R)
      r3(get_estimates_unsuppressed_proportion_participants_vl200.R)
      r4(get_estimates_unsuppressed_proportion_non_participants_vl200.R)
      r5(get_treatment_cascade_participants_vl200.R)
      r6(get_treatment_cascade_non_participants_vl200.R)
    end
    style misc fill:#F5F5F5

    subgraph fit
      f1[(RCCS_art_posterior_samples_vl200_221208.rds)]
      f2[(RCCS_art_posterior_samples_newlyregistered_vl200_221208.rds)]
      f3[(RCCS_nonsuppressed_proportion_posterior_samples_vl_200_221121.rds)]
      f4[(RCCS_nonsuppressed_proportion_posterior_samples_vl_200_newlyregistered_221121.rds)]

      f5[(RCCS_treatment_cascade_participants_posterior_samples_vl200_221208.rds)]
      f6[(RCCS_treatment_cascade_participants_estimates_vl200_221208.csv)]

      f7[(RCCS_treatment_cascade_nonparticipants_posterior_samples_vl200_221208.rds)]
      f8[(RCCS_treatment_cascade_nonparticipants_estimates_vl200_221208.csv)]
    end
    style f1 fill:#C0C0C0
    style f2 fill:#C0C0C0
    style f3 fill:#C0C0C0
    style f4 fill:#C0C0C0
    
    style f5 fill:#FFD700
    style f6 fill:#FFD700
    style f7 fill:#FFD700
    style f8 fill:#FFD700
    style fit fill:#F5F5F5
    
    d1 -- load --> r1
    r1 -- create --> f1

    d2 -- load --> r2
    r2 -- create --> f2

    d3 -- load --> r3
    r3 -- create --> f3

    d4 -- load --> r4
    r4 -- create --> f4

    d5 -- load --> r5
    f1 -- load --> r5
    f3 -- load --> r5
    r5 -- create --> f5
    r5 -- create --> f6

    d5 -- load --> r6
    f2 -- load --> r6
    f4 -- load --> r6
    r6 -- create --> f7
    r6 -- create --> f8
```