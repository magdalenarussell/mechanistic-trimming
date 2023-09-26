import sys
sys.path.append('/home/mrussel2/microhomology/jax_scripts/')
import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import importlib
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_model_classes import ConditionalLogisticRegressor, ConditionalLogisticRegressionPredictor
from config import MOD_OUTPUT_PATH
import variable_configuration

ANNOTATION_TYPE = sys.argv[1]
PARAM_GROUP = sys.argv[2]
LEFT_NUC_MOTIF_COUNT = int(sys.argv[3])
RIGHT_NUC_MOTIF_COUNT = int(sys.argv[4])
MODEL_TYPE = sys.argv[5]
L2 = bool(sys.argv[6])

# set global variables
param_config = importlib.import_module(f"param_groups.{PARAM_GROUP}")

params = variable_configuration.global_paramaters(param_config, MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, LEFT_NUC_MOTIF_COUNT, RIGHT_NUC_MOTIF_COUNT, MODEL_TYPE).set()

# read in data
processed_data = pd.read_csv(params.R_processed_data_path(), sep = '\t')

#TODO set vars 
## TODO write model configs maybe a class file for each model type that sets variables?
# initialize model 
model = ConditionalLogisticRegressor(training_df = processed_data,
                                     variable_colnames=['height','width', 'rating'],
                                     count_colname='count',
                                     group_colname='color',
                                     repeat_obs_colname='garden',
                                     choice_colname='species')

# TODO write trained model (pickle?)
# train model
model = model.train_model(l2=L2)

# example of getting results
predictor = ConditionalLogisticRegressionPredictor(model=model,
                                                   new_df=processed_data,
                                                   variable_colnames=['height','width', 'rating'],
                                                   count_colname='count',
                                                   group_colname='color',
                                                   repeat_obs_colname='garden',
                                                   choice_colname='species')

#TODO write predictions, loss, coefficients
y_hat = predictor.predict()
print(y_hat)
loss = predictor.compute_loss()
print(loss)
c = predictor.get_coefficients()
print(c)

