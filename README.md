# mechanistic-trimming 
The goal of this project is to use model-based statistical inference to identify the sequence-level determinants of nucleotide trimming during V(D)J recomination of adaptive immune receptor loci.

# Install
Everything R 4.1.3 based. R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
conda env create -f environment.yml
conda activate mechanistic-trimming 
```

You will also need to install IGoR if you wish to annotate sequences using IGoR (Marcou et.al Nature Communications 2018)

# Requirements: 
Most of these analyses can be run on any machine.
However, some of the data preparation steps, such as sequence annotation using IGoR (Marcou et.al Nature Communications 2018), are computationally intensive and require a cluster to run efficiently.
This sequence annotation script is written specifically for a cluster set up to use the Slurm job scheduler. 
(Some minor modifications to the [sequence annotation script](annotate_with_igor.sh) could allow this step to be run locally or using a different cluster workload manager. 

# Analysis outline: 

__Table of Contents:__

* [Model training](#model-training)
    * [Summary of model types](#model-types)
* [Quantifying coefficient significance](#quantifying-coefficient-significance)
* [Model evaluation](#model-evaluation)
    * [Summary of model evaluation options](#model-evaluation-options)
* [Model validation](#model-validation)
    * [Summary of model validation options](#model-validation-options)
* [Plot results](#plot-results)
* [Supplementary analyses](#supplementary-analyses)
    * [Meaures the relative weights of parameters across different testing data sets](#parameter-relative-weight-analysis)
    * [Artemis SNP analysis](#artemis-snp-analysis)

## Model training

0. Download the training cohort data set using the [link](https://doi.org/10.21417/B7001Z) provided in the original publication
1. Annotate these sequences using IGoR. You can run the [annotation script](annotate_with_igor.sh) to do this. As described above, this script is written specifically for a cluster set up to use the Slurm job scheduler. This script takes five arguments:
    1. the raw file path (wherever you downloaded the files in step 0 above)
    2. a temporary directory location (wherever you want the intermediate annotation files to be stored)
    3. an output directory location (wherever you want the final annotation files to be stored--you will eventually save this location as `TCR_REPERTOIRE_DATA_igor` within the [config](config/config.R) file in step 3)
    4. the number of possible rearrangment scenarios you want to sample from for each sequence (we used 10 scenarios)
    5. the number of CPUs you want to use
2. Download TRB gene name and germline sequences from IMGT (we downloaded these data in December 2021); save these data to a file--you will eventually save the location of this file as `WHOLE_NUCSEQS_igor` within the [config](config/config.R) file in step 3 
3. Edit [config](config/config.R) file to be project and/or computer specific. See the [README](config/README.md) for more details.
4. Train model using the [model fitting script](fit_model.sh). This script takes 12 arguments, and can be run locally or on a cluster: 
    1. the annotation type--for the general model, this should be `igor`
    2. the trimming type for model fitting (either `v_trim` or `j_trim`)
    3. productivity of sequences (either `productive`, `nonproductive`, or `both`)
    4. motif type--for the general model, this should be `unbounded`
    5. the number of CPUs you want to use
    6. model group (either `all_subjects` or `individual_subjects` depending on whether you want to train a model across all subjects, or a unique model for each individual)
    7. gene weight type (either `p_gene_marginal` or `p_gene_given_subject`); for the general model use `p_gene_marginal`--this will ensure that all genes are weighted the same across all subjects in the construction of the likelihood
    8. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    9. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    10. upper trim bound--the upper bound on trimming lengths; for the general model, we use a bound of 14
    11. model type--see the table below for [model options](#model-types)
    12. (optional) for models with base-count terms, specify the number of nucleotides to be included in the base-count 5' of the trimming site; for the general model, we use 10 nucleotides

This trained model will be stored in the [models](models/) directory. 

### Model types

Here is a summary of some of the model options. You can find a description of additional model type options [here](scripts/model_formula_functions/README.md)

| Long model name (same as manuscript) | Code argument                      | Features parameterized                                                                                                                        | Other notes                                                                                                                                                                                                                    |
|--------------------------------------|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| _null_                               | `null`                             | none                                                                                                                                          |                                                                                                                                                                                                                                |
| _motif_                              | `motif`                            | a sequence motif                                                                                                                              | the motif size is specified by the `left_motif_count` and `right_motif_count` arguments                                                                                                                                        |
| _DNA-shape_                          | `dna_shape-std`                    | DNA-shape parameters for nucleotides contained in a specified sequence window                                                                 | the sequence window size is specified by the `left_motif_count` and `right_motif_count` arguments                                                                                                                              |
| _distance_                           | `linear-distance`                  | sequence-independent distance from the end of the gene (integer-valued)                                                                       |                                                                                                                                                                                                                                |
| _two-side base-count_                | `two-side-base-count`              | count of AT and GC nucleotides on either side of the trimming site                                                                            | the total number of nucleotides included in the 5' AT and GC counts is specified by the 12th argument of the model training script                                                                                             |
| _motif + distance_                   | `motif_linear-distance`            | a sequence motif and sequence-independent distance from the end of the gene (integer-valued)                                                  | the motif size is specified by the `left_motif_count` and `right_motif_count` arguments                                                                                                                                        |
| _motif + two-side base-count beyond_ | `motif_two-side-base-count-beyond` | a sequence motif and the count of AT and GC nucleotides on either side of the trimming site (including only nucleotides outside of the motif) | the motif size is specified by the `left_motif_count` and `right_motif_count` arguments and the total number of nucleotides included in the 5' AT and GC counts is specified by the 12th argument of the model training script |

## Quantifying coefficient significance

If you would like to measure the significance of each model coefficient, you can use [this script](evaluate_coef_significance.sh).
This script takes 11 arguments:

1. the annotation type--for the general model, this should be `igor`
2. the trimming type for model fitting (either `v_trim` or `j_trim`)
3. productivity of sequences (either `productive`, `nonproductive`, or `both`)
4. motif type--for the general model, this should be `unbounded`
5. the number of CPUs you want to use
7. gene weight type (either `p_gene_marginal` or `p_gene_given_subject`); for the general model use `p_gene_marginal`--this will ensure that all genes are weighted the same across all subjects in the construction of the likelihood
8. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
9. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
10. upper trim bound--the upper bound on trimming lengths; for the general model, we use a bound of 14
11. model type--see the table below for [model options](#model-types)
12. (optional) for models with base-count terms, specify the number of nucleotides to be included in the base-count 5' of the trimming site; for the general model, we use 10 nucleotides

This analysis will save a file containing the results in the indicated `OUTPUT_PATH` as specified in the [config](config) files

## Model evaluation

If you would like to evaluate the performance of a model using various subsets of the training data set, you can use the [model evaluation script](evaluate_model.sh). This script also takes 12 arguments which are similar to the model training script:

1. the annotation type--for the general model, this should be `igor`
2. the trimming type for model fitting (either `v_trim` or `j_trim`)
3. productivity of sequences (either `productive`, `nonproductive`, or `both`)
4. motif type--for the general model, this should be `unbounded`
5. the number of CPUs you want to use
7. gene weight type (either `p_gene_marginal` or `p_gene_given_subject`); for the general model use `p_gene_marginal`--this will ensure that all genes are weighted the same across all subjects in the construction of the likelihood
8. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
9. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
10. upper trim bound--the upper bound on trimming lengths; for the general model, we use a bound of 14
11. model type--see the table above for [model options](#model-types)
11. model evaluation type--see the table below for [evaluation options](#model-evaluation-options)
12. (optional) for models with base-count terms, specify the number of nucleotides to be included in the base-count 5' of the trimming site; for the general model, we use 10 nucleotides

**Note: all output files will be located at the indicated `OUTPUT_PATH` as specified in the [config](config) files**

### Model evaluation options

Here is a summary of the model evaluation arguments. More details can be found [here](scripts/model_evaluation_type_functions/README.md)

| Long evaluation name (same as manuscript)                  | Code argument             | Summary of held-out group                                                                                                                         |
|------------------------------------------------------------|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| _full V-gene training data set_                            | `log_loss`                | Loss computed across the full training data set                                                                                                   |
| _many held-out subsets of the V-gene training data set_    | `expected_log_loss`       | Expected loss computed across 20 random different held-out subsets of the full training data set                                                  |
| _"most different" cluster of V-genes (terminal sequences)_ | `v_gene_family_loss`      | Loss computed across the group of V-genes which are most-different according to the hamming distances of the last 25 nucleotides of each sequence |
| _"most different" cluster of V-genes (full sequences)_     | `full_v_gene_family_loss` | Loss computed across the group of V-genes which are most-different according to the hamming distances of each full sequence                       |
| _full J-gene data set_                                     | `log_loss_j_gene`         | Loss computed across the full J-gene data set from the training data set                                                                          |

## Model validation

If you would like to use the model to make predictions on a new data set and/or validate the model on a testing data set, you can follow these steps:
    
1. Download the processed testing data set, and make sure it contains the same column names as the training data set 
2. Depending on the locus of the testing data, download the appropriate gene name and germline sequences from IMGT (we downloaded these data in August 2022); save these data to a file--you will save the location of this file in the next step
3. Edit [config](config/config.R) file to be project and/or computer specific (be sure to add the path to the file from the last step)
2. Run the [model validation script](validate_model.sh). This script takes 17 arguments which force you to specify model parameters and validation data set specifics:

    1. the annotation type--for the general model, this should be `igor`
    2. the trimming type used in model fitting (either `v_trim` or `j_trim`)
    3. productivity of sequences used in model fitting (either `productive`, `nonproductive`, or `both`)
    4. motif type--for the general model, this should be `unbounded`
    5. the number of CPUs you want to use
    6. model group (either `all_subjects` or `individual_subjects` depending on whether the model was trained across all subjects, or uniquely for each individual)
    7. gene weight type used for model training (either `p_gene_marginal` or `p_gene_given_subject`); for the general model use `p_gene_marginal`--this will ensure that all genes are weighted the same across all subjects in the construction of the likelihood
    8. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    9. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
    10. upper trim bound--the upper bound on trimming lengths; for the general model, we use a bound of 14
    11. model type--see the table below for model options
    12. for models with base-count terms, specify the number of nucleotides to be included in the base-count 5' of the trimming site; for the general model, we use 10 nucleotides (you can specify `NA` otherwise)
    13. the directory storing the validation data set
    14. the validation type--see the [summary](#model-validation-options) below for details
    15. validation trimming type (either `v_trim` or `j_trim`)--this does not need to be the same as what the model was trained with
    16. productivity of sequences for model validation (either `productive`, `nonproductive`, or `both`)
    17. gene weight type used for loss calculation (either `p_gene_marginal` or `p_gene_given_subject`); for the general use `p_gene_given_subject`

**Note: all output files will be located at the indicated `OUTPUT_PATH` as specified in the [config](config) files**

### Model validation options

Summary of the model validation data sets used in our analysis.

| Locus name  | Code argument           | Data used in the manuscript                                                                                                                                                                                                                                 |
|-------------|-------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| _TRB locus_ | `validation_data_beta`  | TRB testing data set; processed data can be downloaded [here](https://doi.org/10.5281/zenodo.5719516)                                                                                                                                                       |
| _TRA locus_ | `validation_data_alpha` | TRA testing data set; processed data can be downloaded [here](https://doi.org/10.5281/zenodo.5719516)                                                                                                                                                       |
| _TRG locus_ | `validation_data_gamma` | TRG testing data set; processed data can be downloaded [here](https://doi.org/10.21417/B7SG6T)                                                                                                                                                              |
| _TRD locus_ | `validation_data_delta` | TRD not analyzed in the manuscript                                                                                                                                                                                                                          |
| _IGH locus_ | `validation_data_igh`   | productive IGH testing data set can be downloaded [here](https://figshare.com/articles/preprint/Functional_antibodies_exhibit_light_chain_coherence/19617633); nonproductive IGH testing data set can be downloaded [here](www.github.com/briney/grp_paper) |

## Plot results

Plot [figures](plotting_scripts/manuscript_plots) from the manuscript.

Also, have a look at the plotting [README](plotting_scripts/README.md) for more details.

## Supplementary analyses

### Parameter relative weight analysis

If you would like to measure the relative weights of the motif and base-count-beyond parameters within the _motif + two-side base-count-beyond_ model for a specific validation data set, you can run [this script](evaluate_rel_importance.sh) locally or on a cluster. This script takes 15 arguments:

1. the annotation type--for the general model, this should be `igor`
2. the trimming type used in model fitting (either `v_trim` or `j_trim`)
3. productivity of sequences used in model fitting (either `productive`, `nonproductive`, or `both`)
4. motif type--for the general model, this should be `unbounded`
5. the number of CPUs you want to use
6. model group (either `all_subjects` or `individual_subjects` depending on whether the model was trained across all subjects, or uniquely for each individual)
7. gene weight type used for model training (either `p_gene_marginal` or `p_gene_given_subject`); for the general model use `p_gene_marginal`--this will ensure that all genes are weighted the same across all subjects in the construction of the likelihood
8. left motif count (an integer between 0 and 6)--specifies the number of nucleotides 5' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
9. right motif count (an integer between 0 and 6)--specifies the number of nucleotides 3' of the trimming site to be included in the motif parameter (for models without motif parameters, this should be 0)
10. upper trim bound--the upper bound on trimming lengths; for the general model, we use a bound of 14
12. the number of nucleotides to be included in the base-count 5' of the trimming site; for the general model, we use 10 nucleotides
13. the directory storing the validation data set
14. the validation type--see the [summary](#model-validation-options) for details
15. validation trimming type (either `v_trim` or `j_trim`)--this does not need to be the same as what the model was trained with
16. productivity of sequences for model validation (either `productive`, `nonproductive`, or `both`)

This analysis will save a file containing the results in the indicated `OUTPUT_PATH` as specified in the [config](config) files

### Artemis SNP analysis

If you would like to evaluate whether model coefficients vary in the context of an Artemis locus SNP (rs41298872), you can follow the instructions in the [quantifying coefficient significance](#quantifying-coefficient-significance) section using `motif_two-side-base-count-beyond_snp-interaction-20717849` as the model type. 

# About the analysis

With this analysis, we want to quantify the sequence-level determinants of nucleotide trimming during V(D)J recombination.
See the manuscript for specific model and methods details: 

Russell, Magdalena L., Noah Simon, Philip Bradley, and Frederick A. Matsen IV. 2022. "Statistical inference reveals the role of length, breathing, and nucleotide identity in V(D)J nucleotide trimming." In preparation.
