import sys
sys.path.append('/home/mrussel2/microhomology/jax_scripts/')
import os
import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import importlib
import pickle
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_model_classes import ConditionalLogisticRegressor, ConditionalLogisticRegressionPredictor
from config import MOD_OUTPUT_PATH, MOD_PROJECT_PATH
import variable_configuration

ANNOTATION_TYPE = sys.argv[1]
PARAM_GROUP = sys.argv[2]
LEFT_NUC_MOTIF_COUNT = int(sys.argv[3])
RIGHT_NUC_MOTIF_COUNT = int(sys.argv[4])
MODEL_TYPE = sys.argv[5]
L2 = sys.argv[6]
L2 = (L2.lower() == 'true')
NCPU = int(sys.argv[7])
ANNOTATION_TYPE_VALIDATION = sys.argv[8]

# initialize parallelized pandas
pandarallel.initialize(nb_workers=NCPU, progress_bar=True)

# set global variables
param_config = importlib.import_module(f"param_group_configs.{PARAM_GROUP}")
params = variable_configuration.global_paramaters(param_config,
                                                  MOD_OUTPUT_PATH,
                                                  MOD_PROJECT_PATH,
                                                  ANNOTATION_TYPE,
                                                  PARAM_GROUP,
                                                  LEFT_NUC_MOTIF_COUNT,
                                                  RIGHT_NUC_MOTIF_COUNT,
                                                  MODEL_TYPE)

# set model type specific parameters
model_config = importlib.import_module(f"model_type_configs.{MODEL_TYPE}")
model_params = variable_configuration.model_specific_parameters(param_config,
                                                                model_config,
                                                                LEFT_NUC_MOTIF_COUNT,
                                                                RIGHT_NUC_MOTIF_COUNT)
model_params = model_params.process_model_parameters()
print('loaded parameters')

# read in data
processed_data_filename = params.R_processed_data_path(annotation = ANNOTATION_TYPE_VALIDATION)
processed_data = pd.read_csv(processed_data_filename, sep = '\t')
print('read in data')

# validate model
model_filename = params.model_output_path(L2)
evaluator = ConditionalLogisticRegressionEvaluator(model_filename,
                                                   processed_data)

result = evaluator.compile_results_df(params.left_nuc_motif_count,
                                      params.right_nuc_motif_count,
                                      params.motif_type,
                                      params.gene_weight_type,
                                      params.upper_trim_bound,
                                      params.lower_trim_bound,
                                      params.insertions,
                                      params.model_type,
                                      10,
                                      True)

result['validation_data_type'] = ANNOTATION_TYPE_VALIDATION

# write predictions and coefficients
path = params.model_eval_results_path(L2)

if os.path.isfile(path):
    # Read the file as a Pandas DataFrame
    df = pd.read_csv(path)
    result = pd.concat([df, result], axis = 0, fill_value = 'NA')

result.to_csv(path, sep='\t', index=False)
print('finished validating model')
