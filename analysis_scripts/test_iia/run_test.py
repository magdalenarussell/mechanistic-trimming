import sys
sys.path.append('/home/mrussel2/microhomology/mechanistic-trimming/jax_scripts/')
import os
import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import importlib
import pickle
from scipy.stats import chi2
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_model_classes import ConditionalLogisticRegressor, ConditionalLogisticRegressionPredictor
from config import MOD_OUTPUT_PATH, MOD_PROJECT_PATH
import variable_configuration

ANNOTATION_TYPE = 'igor_mh_sim_alpha_from_uniform_MHprob0'
PARAM_GROUP = 'nonproductive_v-j_trim_ligation-mh'
LEFT_NUC_MOTIF_COUNT = 1
RIGHT_NUC_MOTIF_COUNT = 2
MODEL_TYPE = 'motif_two-side-base-count-beyond_ligation-mh'
L2 = False
NCPU = int(sys.argv[1])

model_filename = params.model_output_path(L2)
subset_filename = os.path.dirname(model_filename) + '/trained_model_L2False_missing_1_ligation_mh_cases.pkl'

with open(model_filename, 'rb') as file:
    full_model = pickle.load(file)

with open(subset_filename, 'rb') as file:
    subset_model = pickle.load(file)

# Calculate Hausman test statistic
## get parameter estimates
beta_full = full_model.coefs.flatten()
beta_reduced = subset_model.coefs.flatten()

## Covariance matrices
cov_full = full_model.cov_matrix
cov_reduced = subset_model.cov_matrix

## Hausman test statistic calculation
beta_diff = beta_full - beta_reduced
cov_diff = cov_full - cov_reduced
inv_cov_diff = np.linalg.inv(cov_diff)
H_statistic = np.dot(np.dot(beta_diff.T, inv_cov_diff), beta_diff)

print("Hausman Test Statistic:", H_statistic)

# calculate P=value
df = len(beta_full)
p_value = chi2.sf(H_statistic, df)

print("P-Value:", p_value)
