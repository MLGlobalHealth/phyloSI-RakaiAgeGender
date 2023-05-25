### Data Preprocessing
The flowchart below illustrates in detail how the data is consumed and processed by various R scripts at each stage. The coloured silos, rounded boxes, and grey boxes, represent data, R scripts, and directories respectively. The data is coloured according to the [medallion architecture](https://www.databricks.com/glossary/medallion-architecture), where bronze represents raw data, silver represents cleansed and processed data, and gold represents data which are consumption ready (i.e. directly used within statistical models).

#### Preprocessing Stage 1
```mermaid
  flowchart LR
    subgraph data
      d5[(aggregated_count_hiv_positive.csv)]
      d1[(aggregated_participants_count_art_coverage.csv)]
      d2[(aggregated_newlyregistered_count_art_coverage.csv)]
      d3[(aggregated_participants_count_unsuppressed.csv)]
      d4[(aggregated_newlyregistered_count_unsuppressed.csv)]
    end
    style d1 fill:#CD7F32
    style d2 fill:#CD7F32
    style d3 fill:#CD7F32
    style d4 fill:#CD7F32
    style d5 fill:#CD7F32
    style data fill:#F5F5F5

    subgraph surveillance_pipeline_src
      r5(get_estimates_prevalence.R)
      r1(get_estimate_art_coverage_participants.R)
      r2(get_estimate_art_coverage_non_participants.R)
      r3(get_estimates_unsuppressed_proportion_participants.R)
      r4(get_estimates_unsuppressed_proportion_non_participants.R)
    end
    style surveillance_pipeline_src fill:#F5F5F5

    subgraph preprocessed_outputs_RCCS_surveillance
      f5[(RCCS_prevalence_posterior_sample_221116.rds)]
      
      f1[(RCCS_art_posterior_samples_221208.rds)]
      f2[(RCCS_art_posterior_samples_newlyregistered_221208.rds)] 
      f3[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds)]
      f4[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds)]      
    end

    subgraph final_RCCS_surveillance
      f6[(RCCS_prevalence_estimates_221116.csv)]
    end

    style f1 color:#000000
    style f2 color:#000000
    style f3 color:#000000
    style f4 color:#000000
    style f5 color:#000000
    style f6 color:#000000
    style f1 fill:#C0C0C0
    style f2 fill:#C0C0C0
    style f3 fill:#C0C0C0
    style f4 fill:#C0C0C0
    style f5 fill:#C0C0C0
    style f6 fill:#FFD700
    style preprocessed_outputs_RCCS_surveillance fill:#F5F5F5

    d1 -- load --> r1
    r1 -- create --> f1

    d2 -- load --> r2
    r2 -- create --> f2

    d3 -- load --> r3
    r3 -- create --> f3

    d4 -- load --> r4
    r4 -- create --> f4

    d5 -- load --> r5
    r5 -- create --> f5
    r5 -- create --> f6
```

#### Preprocessing Stage 2
```mermaid
  flowchart LR
    subgraph data
      d3[(sensitivity_specificity_art.csv)]
      d1[(RCCS_participation_221208.csv)]
    end
    style d1 fill:#CD7F32
    style d3 fill:#CD7F32
    style data fill:#F5F5F5

    subgraph surveillance_pipeline_src
      r_part(get_treatment_cascade_participants.R)
      r_npart(get_treatment_cascade_non_participants.R)
      r_pop(get_treatment_cascade_population.R)
    end
    style surveillance_pipeline_src fill:#F5F5F5

    subgraph preprocessed_outputs_RCCS_surveillance
      f_art_post[(RCCS_art_posterior_samples_221208.rds)]
      f_art_post_new[(RCCS_art_posterior_samples_newlyregistered_221208.rds)]
      f_nsup_prop_post[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds)]
      f_nsup_prop_post_new[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds)]
    end
    
    subgraph final_RCCS_surveillance
      f_part_est[(RCCS_treatment_cascade_participants_estimates_221208.csv)]
      f_npart_est[(RCCS_treatment_cascade_nonparticipants_estimates_221208.csv)]
      f_pop_est[(RCCS_treatment_cascade_population_estimates_221208.csv)]
    end

    style f_nsup_prop_post color:#000000
    style f_art_post color:#000000
    style f_nsup_prop_post_new color:#000000
    style f_art_post_new color:#000000
    style f_nsup_prop_post fill:#C0C0C0
    style f_art_post fill:#C0C0C0
    style f_nsup_prop_post_new fill:#C0C0C0
    style f_art_post_new fill:#C0C0C0
  


    style f_pop_est color:#000000
    style f_part_est color:#000000
    style f_npart_est color:#000000
    style f_pop_est fill:#FFD700
    style f_part_est fill:#FFD700
    style f_npart_est fill:#FFD700
    
    style preprocessed_outputs_RCCS_surveillance fill:#F5F5F5

    f_art_post -- load --> r_part
    f_nsup_prop_post -- load --> r_part
    d3 -- load --> r_part
    r_part -- create --> f_part_est

    f_art_post_new -- load --> r_npart
    f_nsup_prop_post_new -- load --> r_npart
    d3 -- load --> r_npart
    r_npart -- create --> f_npart_est

    d1 -- load --> r_pop
    f_art_post -- load --> r_pop
    f_nsup_prop_post -- load --> r_pop
    f_art_post_new -- load --> r_pop
    f_nsup_prop_post_new -- load --> r_pop
    d3 -- load --> r_pop
    r_pop -- create --> f_pop_est
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
  
    style f1 color:#000000
    style f2 color:#000000
    style f3 color:#000000
    style f4 color:#000000

    style f1 fill:#C0C0C0
    style f2 fill:#C0C0C0
    style f3 fill:#C0C0C0
    style f4 fill:#C0C0C0
  
    style f5 color:#000000
    style f6 color:#000000
    style f7 color:#000000
    style f8 color:#000000
    
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
