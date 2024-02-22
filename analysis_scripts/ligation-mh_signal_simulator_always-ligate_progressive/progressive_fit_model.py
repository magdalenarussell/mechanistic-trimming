import sys
sys.path.append('/home/mrussel2/microhomology/mechanistic-trimming/jax_scripts/')
import glob
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
processed_data_filename = params.R_processed_data_path()
basename = os.path.dirname(processed_data_filename)

# Pattern to match files
pattern = 'progressively*.tsv'
search_pattern = os.path.join(basename, pattern)
files = glob.glob(search_pattern)

# sort files
root_files = [os.path.basename(i)[25:-4] for i in files]
g, num= zip(*(item.split('_') for item in root_files))
g = list(g)
num = [int(i) for i in num]

file_info = {'file_path' : files, 'gene_pair_count' : num, 'Vgene_addition': g}
file_info_df = pd.DataFrame(file_info)
file_info_df = file_info_df.sort_values('gene_pair_count').reset_index(drop = True)

prog_data = pd.DataFrame()
if MODEL_TYPE == 'ligation-mh':
    old_coefs = jnp.array([[0.01]])
elif MODEL_TYPE == 'motif_two-side-base-count-beyond_ligation-mh':
    old_coefs = jnp.full((25,1), 0.01)

for file in file_info_df.file_path:
    processed_data = pd.read_csv(file, sep = '\t')
    print('read in data: ', file)
    # initialize model 
    model = ConditionalLogisticRegressor(training_df = processed_data,
                                         variable_colnames = model_params.variable_colnames,
                                         count_colname = model_params.count_colname,
                                         group_colname = model_params.group_colname,
                                         repeat_obs_colname = model_params.repeat_obs_colname,
                                         choice_colname = model_params.choice_colname,
                                         params = params)
    model.initial_coefs = old_coefs
    print('initialized model')
    # train model
    model = model.train_model(l2=L2, maxiter=1000, tolerance=1e-8)
    old_coefs = model.coefs
    print('trained model')
    coefs = model.get_coefficients_df()
    coefs['progressive_gene_pair_count'] = int(file_info_df[file_info_df.file_path == file].gene_pair_count.iloc[0])
    coefs['training_error'] = float(model.training_info.state.error)
    prog_data = pd.concat([prog_data, coefs], axis=0).reset_index(drop=True)


coefs_filename = params.trained_coefs_path(L2)
coefs_filename2 = os.path.dirname(coefs_filename) + '/progressively_trained_coefs_L2False.tsv'

prog_data.to_csv(coefs_filename2, sep='\t', index=False)
print('finished processing model predictions')
