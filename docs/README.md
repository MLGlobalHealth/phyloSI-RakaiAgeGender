## Table of Contentes
- [Data Repository Structure](#data-repository-structure)
- [Code Repository Structure](#code-repository-structure)
- [RCCS Surveillance Analyses](#rccs-surveillance-analyses)
  - [Estimate HIV status and prevalence](#estimate-hiv-status-and-prevalence)
  - [Estimate ART use](#estimate-art-use)
  - [Estimate viral suppresion](#estimate-viral-suppresion)
  - [Estimate treatment cascade](#estimate-treatment-cascade)


### Data Repository Structure
We provide all the data required to produce the findings in our paper via Zenodo (TODO: Add link). The directory is structured as follows:
<pre>
.
├── data_RCCS_surveillence
├── data_incidence_rate
├── data_deep_sequence_phylogenies
├── preprocessed_outputs_RCCS_surveillence
├── preprocessed_outputs_deep_sequence_phylogenies
├── final_RCCS_surveillence
├── final_incidence_rate
├── final_sexual_partnership
├── final_transmission_flow
└── final_central_analysis
</pre>
`data_*` directories contain de-identified data which serves as the starting point of all analysis. `preprocessed_*` directories contain intermediate data files produced by the scripts within the `src` directory. `final_*` contains cleaned data or model outputs which are directly used to produce figures and tables. A list of discriptions for the data can be found at `docs/README.md`.

### GitHub Repository Structure
The contents of the repository and short discriptions are as follows:
<pre>
.
├── R                                                # R helper functions 
├── docs                                             # Documentation
├── figure_src                                       # Scripts to generate paper figures 
├── scripts                                          # R scripts for final analysis
├── src                                              # R scripts for sub-analysis 
├── stan_models                                      # Stan code
├── ...
├── config.R                                         # Path configurations
├── phyloSI-RakaiAgeGender-install.sh                # Installation
├── phyloSI-RakaiAgeGender-install-macosM1.sh        # Installation for Apple M1
├── phyloSI-RakaiAgeGender-preprocessing.sh          # Preprocess data and run sub-analyses
├── phyloSI-RakaiAgeGender-run_phyloflows-hpc.sh     # Run final analysis on Imperial HPC
├── phyloSI-RakaiAgeGender-run_phyloflows-laptop.sh  # Run final analysis on laptop
└── ...
</pre>

### RCCS Surveillance Analysis
The flowcharts below illustrates in detail how the data is consumed and processed by various R scripts at each stage. The silos, rounded boxes, and yellow boxes, represent data, R scripts, and directories respectively.

#### Estimate HIV status and prevalence
```mermaid
flowchart LR
  subgraph data_RCCS_surveillance
    d[(aggregated_count_hiv_positive.csv)]
  end
  
  subgraph src/surveillance_pipeline
    r(get_estimates_prevalence.R)
  end

  subgraph preprocessed_outputs_RCCS_surveillance
    f1[(RCCS_prevalence_posterior_sample_221116.rds)]
  end

  subgraph final_RCCS_surveillance
    f2[(RCCS_prevalence_estimates_221116.csv)]
  end

  d -- load --> r
  r -- create --> f1
  r -- create --> f2
```

#### Estimate ART use
```mermaid
  flowchart LR
    subgraph data_RCCS_surveillance
      d1[(aggregated_participants_count_art_coverage.csv)]
      d2[(aggregated_newlyregistered_count_art_coverage.csv)]
    end

    subgraph src/surveillance_pipeline_src
      r1(get_estimate_art_coverage_participants.R)
      r2(get_estimate_art_coverage_non_participants.R)
    end

    subgraph preprocessed_outputs_RCCS_surveillance
      f1[(RCCS_art_posterior_samples_221208.rds)]
      f2[(RCCS_art_posterior_samples_newlyregistered_221208.rds)]   
    end

    d1 -- load --> r1
    r1 -- create --> f1

    d2 -- load --> r2
    r2 -- create --> f2
```

#### Estimate viral suppresion
```mermaid
  flowchart LR
    subgraph data
      d3[(aggregated_participants_count_unsuppressed.csv)]
      d4[(aggregated_newlyregistered_count_unsuppressed.csv)]
    end

    subgraph src/surveillance_pipeline_src
      r3(get_estimates_unsuppressed_proportion_participants.R)
      r4(get_estimates_unsuppressed_proportion_non_participants.R)
    end

    subgraph preprocessed_outputs_RCCS_surveillance
      f3[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds)]
      f4[(RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds)]      
    end

    d3 -- load --> r3
    r3 -- create --> f3

    d4 -- load --> r4
    r4 -- create --> f4
```

#### Estimate treatment cascade
```mermaid
  flowchart LR
    subgraph data_RCCS_surveillance
      d3[(sensitivity_specificity_art.csv)]
      d1[(RCCS_participation_221208.csv)]
    end

    subgraph surveillance_pipeline_src
      r_part(get_treatment_cascade_participants.R)
      r_npart(get_treatment_cascade_non_participants.R)
      r_pop(get_treatment_cascade_population.R)
    end

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

## Data and script reference

### Sample Data
The table below lists the data files within `/data` and a brief description of its contents

<details>
<summary><b>Click to show table</b></summary>
  <table>
    <thead>
      <tr>
        <th>File name</th>
        <th>Description</th>
      </tr>
    </thead>
  <tbody>
    <tr>
      <td><code>pairsdata_toshare_d1_w11_netfrompairs_postponessrm.rds</code></td>
      <td>HIV source-recipient pairs</td>
    </tr>
    <tr>
      <td><code>RCCS_census_eligible_individuals_221116.csv</code></td>
      <td>Count of census eligible population by age, gender and round</td>
    </tr>
      <tr>
      <td><code>RCCS_participation_221208.csv</code></td>
      <td>Participation rates to the RCCS survey by age, gender and round</td>
      </tr>
    <tr>
      <td><code>aggregated_count_hiv_positive.csv</code></td>
      <td>Count of participants by hiv status, age, gender and round used in <code>surveillance_pipeline_src/get_estimates_prevalence.R</code> to estimate smooth proportion of hiv prevalence among population.</td>
    </tr>
    <tr>
      <td><code>aggregated_participants_count_art_coverage.csv</code></td>
      <td>count of participants by self-reported art use, age, gender and round, used in <code>surveillance_pipeline_src/get_estimates_art_coverage_participants.R</code> to estimate smooth proportion of art coverage among participants.</td>
    </tr>
    <tr>
      <td><code>aggregated_participants_count_unsuppressed.csv</code></td>
      <td>count of participants by viremic viral loads, age, gender and round, used in <code>surveillance_pipeline_src/get_estimates_unsuppressed_proportion_participants.R</code> to estimate smooth proportion of viral suppression among participants.</td>
    </tr>
    <tr>
      <td><code>aggregated_newlyregistered_count_art_coverage.csv</code></td>
      <td>count of first-time participants by self-reported art use, age, gender and round, used in <code>surveillance_pipeline_src/get_estimates_art_coverage_non_participants.R</code> to estimate smooth proportion of art coverage among first-time participants.</td>
    </tr>
    <tr>
      <td><code>aggregated_newlyregistered_count_unsuppressed.csv</code></td>
      <td>count of first-time participants by viremic viral loads, age, gender and round, used in <code>surveillance_pipeline_src/get_estimates_unsuppressed_proportion_non_participants.R</code> to estimate smooth proportion of viral suppression among first-time participants.</td>
    </tr>
    <tr>
      <td><code>seroconverter_cohort_R6R19.rds</code></td>
      <td>Individual-level seroconvert cohort data with information on rounds of enrollment, age, sex, hiv status</td>
      </tr>
      <tr>
      <td><code>Rakai_incpredictions_inland_221107.csv</code></td>
      <td>HIV incidence rates estimates by age, gender and round</td>
    </tr>
    <tr>
      <td><code>inland_R015_cntcts_rate_1130b.rds</code></td>
      <td>Sexual contact rate estimates by age and gender for round 15</td>
    </tr>
       <tr>
      <td><code>Rakai_Pangea2_RCCS_Metadata_randomized.RData  data/sequences_collection_dates_randomized.rds</code></td>
      <td>Individual-level meta data containing Anonymised IDs, round sex, community type, and randomized visit dates,  birthdays, and test dates  </td>
    </tr>
         <tr>
      <td><code>sequences_collection_dates_randomized.rds</code></td>
      <td> Blood samples collection dates randomized within 6 months. </td>
    </tr>
    </tbody>
  </table>
</details>

### Generated Data
<details>
<summary><b>Click to show table</b></summary>
<table>
  <thead>
    <tr>
      <th>File name</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>RCCS_prevalence_estimates_DATE.csv</code></td>
      <td>Smooth proportion of HIV prevalence in population by age, gender and round</td>
    </tr>
    <tr>
      <td><code>RCCS_treatment_cascade_population_estimates_DATE.csv</code></td>
      <td>Smooth estimates of ART coverage and viral suppression in population by age, gender and round</td>
    </tr>
    <tr>
      <td><code>RCCS_treatment_cascade_participants_estimates_DATE.csv</code></td>
      <td>Smooth estimates of ART coverage and viral suppression in participants by age, gender and round</td>
    </tr>
    <tr>
      <td><code>RCCS_treatment_cascade_nonparticipants_estimates_DATE.csv</code></td>
      <td>Smooth estimates of ART coverage and viral suppression in first-time participants by age, gender and round</td>
    </tr>
  </tbody>
</table>
</details>

### Scripts

The following table list the actions performed by each of the scripts in `surveillance_pipeline_src/`:
<details>
  <summary><b>Click to show table</b></summary>
  <table>
    <thead>
      <tr>
        <th>Action</th>
        <th>Script</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>Obtain smooth HIV prevalence by age, gender, and round</td>
        <td><code>get_estimates_prevalence.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth ART coverage in participants by age, gender, and round</td>
        <td><code>get_estimates_art_coverage_participants.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth ART coverage in first-time participants by age, gender, and round</td>
        <td><code>get_estimates_art_coverage_non_participants.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth viral suppression in participants by age, gender, and round</td>
        <td><code>get_estimates_unsuppressed_proportion_participants.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth viral suppression in first-time participants by age, gender, and round</td>
        <td><code>get_estimates_unsuppressed_proportion_non_participants.R</code></td>
      </tr>
      <tr>
        <td>Combine ART coverage and viral suppression estimates to obtain treatment cascade in participants</td>
        <td><code>get_treatment_cascade_participants.R</code></td>
      </tr>
      <tr>
        <td>Combine ART coverage and viral suppression estimates to obtain treatment cascade in first-time participants</td>
        <td><code>get_treatment_cascade_non_participants.R</code></td>
      </tr>
      <tr>
        <td>Combine ART coverage and viral suppression estimates to obtain treatment cascade in population</td>
        <td><code>get_treatment_cascade_population.R</code></td>
      </tr>
    </tbody>
  </table>
</details>
