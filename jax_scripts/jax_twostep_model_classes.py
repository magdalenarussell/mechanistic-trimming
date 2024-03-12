import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
import dill
import pickle
from pandarallel import pandarallel
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold
from jax_model_classes import DataTransformer

class TwoStepDataTransformer(DataTransformer):
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.choice2_colname = choice2_colname
        self.choice1_variable_colnames = choice1_variable_colnames
        self.choice2_variable_colnames = choice2_variable_colnames
        self.original_choice1_variable_colnames = choice1_variable_colnames
        self.original_choice2_variable_colnames = choice2_variable_colnames
        self.original_choice2_colname = choice2_colname
        self.input_choice1_variable_colnames = choice1_variable_colnames
        self.input_choice2_variable_colnames = choice2_variable_colnames
        self.input_choice2_colname = choice2_colname
        self.coefs = None

    def preprocess_data(self, df, pretrain=True):
        """
        Preprocesses the training data, transforming columns as necessary.

        Returns:
            pd.DataFrame: The preprocessed training DataFrame.
        """
        # filter for possible sites
        if 'ligation_mh' in self.input_choice_colname:
            df = self.filter_input_domain_space(df)
        if 'ligation_mh' in self.input_choice2_colname:
            df = self.filter_input_domain_space(df)

        # Transform group column lists into strings
        df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # remove zero counts
        if pretrain:
            df = self.remove_zero_set_counts(df)

        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")

        # Transform choice column lists into strings
        df = self.expand_multivariable(df, "original_choice_colname", "choice_colname")
        df = self.expand_multivariable(df, "original_choice2_colname", "choice2_colname")

        nonrepeat_groups = self.group_colname

        if self.repeat_obs_colname != None:
            self.group_colname = [self.group_colname] + [self.repeat_obs_colname]
            self.original_group_colname = self.group_colname
            df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # Transform group and choice columns into integer type
        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")
        df = self.transform_categorical_response_vars(df, "original_choice_colname", "choice_colname")
        df = self.transform_categorical_response_vars(df, "original_choice2_colname", "choice2_colname")


        # Transform categorical variable columns into contrast columns
        for col in self.variable_colnames:
            df = self.transform_categorical_vars(df, col, pretrain)

        # check for within-choice-set variance
        self.check_within_set_variance(df, pretrain)

        return(df)


    def get_matrices(self, df, pretrain=True, replace_object=None, return_df=False):
        """
        Prepares data matrices for twostep modeling.

        Returns:
            tuple: A tuple containing four data matrices - variables, counts, non-repeat groups, and mask.
        """
        #TODO may need to create a second mask for nonproductive filtering!!!
        df = self.preprocess_data(df, pretrain)

        # Fill in missing counts with zero
        df[self.count_colname] = df[self.count_colname].fillna(0).astype(float)

        # Get matrix shapes
        groups = pd.unique(df[self.group_colname])
        choices = pd.unique(df[self.choice_colname])
        choice2s = pd.unique(df[self.choice2_colname])

        choice1_final_shape = (len(groups), len(choices), len(self.choice1_variable_colnames))
        choice1_int_shape = (len(groups), len(self.choice1_variable_colnames), len(choices))
        choice2_final_shape = (len(groups), len(choices), len(choice2s), len(self.choice2_variable_colnames))
        choice2_int_shape = (len(groups), len(self.choice2_variable_colnames), len(choices), len(choice2s))
        counts_shape = (len(groups), len(choices), len(choice2s), 1)

        # map variables to indices
        choice1_var_mapping = {key: value for key, value in enumerate(self.choice1_variable_colnames)}
        choice1_var_mapping_df = pd.DataFrame(list(choice1_var_mapping.items()), columns=['VarIndex', 'VarName'])
        choice2_var_mapping = {key: value for key, value in enumerate(self.choice2_variable_colnames)}
        choice2_var_mapping_df = pd.DataFrame(list(choice2_var_mapping.items()), columns=['VarIndex', 'VarName'])


        # fill in df with all combos
        df['indicator'] = 1.0
        all_combinations = pd.MultiIndex.from_product([choices, choice2s],
                                                      names=[self.choice_colname, self.choice2_colname]).to_frame(index=False)
        df_full = pd.merge(all_combinations, df, on=[self.choice_colname, self.choice2_colname], how='left').fillna(0)

        # get counts matrix
        pivot_counts = df_full.pivot_table(index=[self.group_colname],
                                      columns=[self.choice_colname, self.choice2_colname],
                                      values=self.count_colname,
                                      fill_value=0)
        counts_mat = jnp.array(pivot_counts).reshape(counts_shape)

        # get variable matrix
        choice1_pivot_vars = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname],
                                    values=choice1_var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        choice1_pivot_vars = choice1_pivot_vars.reindex(choice1_var_mapping_df.VarName.tolist(), axis=1, level = 0)

        choice2_pivot_vars = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname, self.choice2_colname],
                                    values=choice2_var_mapping_df.VarName.tolist(),
                                    fill_value=0)
        # reorder columns to reflect correct order
        choice2_pivot_vars = choice2_pivot_vars.reindex(choice2_var_mapping_df.VarName.tolist(), axis=1, level = 0)

        choice1_mat = jnp.array(choice1_pivot_vars).reshape(choice1_int_shape)
        choice1_mat = choice1_mat.transpose((0, 2, 1))
        assert choice1_mat.shape == choice1_final_shape, "choice1 variable matrix is the incorrect dimension"

        choice2_mat = jnp.array(choice2_pivot_vars).reshape(choice2_int_shape)
        choice2_mat = choice2_mat.transpose((0, 2, 3, 1))
        assert choice2_mat.shape == choice2_final_shape, "choice2 variable matrix is the incorrect dimension"

        # Assuming groups is a list of unique group identifiers
        group_indices = {group: idx for idx, group in enumerate(groups)}
        nonrepeat_groups_mat = jnp.array([group_indices[group] for group in df[self.group_colname].unique()])

        # get mask matrix (1 for valid entries or 0 for nonvalid)
        pivot_mask = df_full.pivot_table(index=[self.group_colname],
                                    columns=[self.choice_colname, self.choice2_colname],
                                    values=['indicator'],
                                    fill_value=0)
        mask_mat = jnp.array(pivot_mask).reshape(counts_shape)

        if replace_object is not None:
            setattr(self, replace_object, df)

        if return_df:
            return choice1_mat, choice2_mat, counts_mat, nonrepeat_groups_mat, mask_mat, df
        else:
            return choice1_mat, choice2_mat, counts_mat, nonrepeat_groups_mat, mask_mat

    def get_random_coefs(self):
        """
        Generates random coefficients for model initialization.

        Returns:
            ndarray: Randomly generated coefficients as a NumPy array.
        """
        # Set random coefficients for model initialization
        key = jax.random.PRNGKey(123)
        coefs = jax.random.normal(key, shape=(len(self.choice1_variable_colnames) + len(self.choice2_variable_colnames), 1))
        return coefs


class TwoStepConditionalLogisticRegressor(TwoStepDataTransformer):
    def __init__(self, training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params, l2kfold=10):
        super().__init__(training_df, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params)
        if training_df != None:
            self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix, self.mask_matrix = self.get_matrices(training_df, replace_object='training_df')
        self.initial_coefs = self.get_random_coefs()
        self.coefs = None
        self.training_info = None
        self.maxiter = 1000
        self.tolerance = 1e-6
        self.step = 0.1
        self.l2reg = 0
        self.l2kfold = None
        self.l2reg_grid = None

    # Get probability for input parameters given coefficients
    def get_indiv_prob(self, choice1_variables, choice2_variables, mask, coefs=None):
        if coefs is None:
            coefs = self.coefs

        assert coefs is not None, "Need to train the model before making predictions!"

        choice1_coefs = coefs[0:len(self.choice1_variable_colnames)]
        choice2_coefs = coefs[len(self.choice1_variable_colnames):]

        # Compute the logits for each choice
        choice1_cov = jnp.dot(choice1_variables, choice1_coefs)
        choice1_reshape = jnp.squeeze(choice1_cov)
        choice2_cov = jnp.dot(choice2_variables, choice2_coefs)
        choice2_reshape = jnp.squeeze(choice2_cov)
        reshaped_mask = jnp.squeeze(mask)

        # Calculate the probability of the observed choices
        # Dimensions of this matrix are groups x choices
        # replace missing choices with -INF so that they will not count towards probability
        choice1_probs = jax.nn.softmax(choice1_reshape)
        choice2_probs = jax.nn.softmax(jnp.where(reshaped_mask, choice2_reshape, jnp.NINF))
        # choice2_probs = jnp.where(jnp.isnan(choice2_probs), 0, choice2_probs)

        return choice1_probs, choice2_probs

    def get_joint_prob(self, choice1_variables, choice2_variables, mask, coefs=None):
        choice1_probs, choice2_probs = self.get_indiv_prob(choice1_variables, choice2_variables, mask, coefs)
        choice1_probs_reshape = choice1_probs.reshape(choice1_probs.shape[0], choice1_probs.shape[1], 1)
        prob = choice1_probs_reshape * choice2_probs
        return prob

    # Get cross-entropy loss
    def cross_entropy(self, choice1_probs, choice2_probs, counts):
        """
        Calculate the cross-entropy loss.

        Args:
            probs (ndarray): Choice probabilities.
            counts (ndarray): Counts of choices.

        Returns:
            float: Cross-entropy loss.
        """
        counts1 = counts.sum(axis = 2)
        counts2 = counts
        counts1_reshape = jnp.squeeze(counts1)
        counts2_reshape = jnp.squeeze(counts2)

        loss1 = -jnp.nansum(jnp.log(jnp.where(choice1_probs==0, 1, choice1_probs)) * counts1_reshape)
        loss2 = -jnp.nansum(jnp.log(jnp.where(choice2_probs==0, 1, choice2_probs)) * counts2_reshape)

        loss = loss1 +  loss2

        return loss

    def l2regularization(self, coefs, size, l2reg):
        """
        Compute L2 regularization term for only non-base-count and non-motif coefficients.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            size (int): Size of the training data
            l2reg (float): L2 regularization strength.

        Returns:
            float: L2 regularization term.
        """
        def calculate_coef_sum(choice_var_colnames, coefs):
            var_list = [var for var in choice_var_colnames if 'base_count' not in var or 'interaction' in var]
            var_list = [var for var in var_list if 'motif' not in var]

            mh_list = [var in var_list for var in choice_var_colnames]
            mh_jnp = jnp.array(mh_list).reshape(-1, 1)

            coef_subset = coefs * mh_jnp

            c = jnp.nansum(coef_subset**2)
            return c

        c = 0
        # add coefs from choice1 variables
        c += calculate_coef_sum(self.choice1_variable_colnames + self.choice2_variable_colnames, coefs)
        return(0.5*(1/size)*l2reg*c)

    # Compute the loss function
    def loss_fn(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
        """
        Compute the loss function for optimization.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            variables (ndarray): Data matrix of variables.
            counts (ndarray): Counts of choices.
            mask (ndarray): Data matrix of mask.
            l2reg (float): L2 regularization strength.

        Returns:
            float: Total loss, including cross-entropy and L2 regularization.
        """
        choice1_probs, choice2_probs = self.get_indiv_prob(choice1_variables, choice2_variables, mask, coefs)
        if choice1_probs.ndim == 1:
            choice1_probs = choice1_probs.reshape((choice1_probs.shape[0], 1))
        if choice2_probs.ndim == 1:
            choice2_probs = choice2_probs.reshape((choice2_probs.shape[0], 1))

        size = jnp.count_nonzero(mask).item()

        loss = self.cross_entropy(choice1_probs, choice2_probs, counts) + self.l2regularization(coefs, size, l2reg)
        return loss

    def fit(self, choice1_variable_matrix, choice2_variable_matrix, counts_matrix, mask_matrix, l2reg, maxiter, tol, step, initial_coefs):
        """
        Fit the conditional logistic regression model using gradient descent.

        Args:
            variable_matrix (ndarray): Data matrix of variables.
            counts_matrix (ndarray): Counts of choices.
            mask_matrix (ndarray): Data matrix of masks.
            l2reg (float): L2 regularization strength.
            maxiter (int): Maximum number of optimization iterations.
            tolerance (float): tolerance for optimization.
            initial_coefs (ndarray): Initial coefficients for model training.

        Returns:
            OptimizationResult: Result of the optimization process.
        """
        assert counts_matrix is not None, "counts column is missing"

        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.BFGS(fun=self.loss_fn, maxiter=maxiter, tol=tol, verbose=True)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         choice1_variables=choice1_variable_matrix,
                         choice2_variables=choice2_variable_matrix,
                         counts=counts_matrix,
                         mask=mask_matrix,
                         l2reg=l2reg)
        return(res)

    def cv_loss(self, fold_count, l2reg):
        assert self.counts_matrix is not None, "counts column is needed"

        kf = GroupKFold(n_splits=fold_count)
        scores = []

        for train_index, val_index in kf.split(X=self.choice1_variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
            train_c1_data, train_c2_data, val_c1_data, val_c2_data = self.choice1_variable_matrix[train_index], self.choice2_variable_matrix[train_index], self.choice1_variable_matrix[val_index], self.choice2_variable_matrix[val_index]
            train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]
            train_mask, val_mask = self.mask_matrix[train_index], self.mask_matrix[val_index]
            train_counts = self.reset_weighted_observations(train_counts)
            val_counts = self.reset_weighted_observations(val_counts)

            # Train the model on the training data
            model = self.fit(train_c1_data, train_c2_data, train_counts, train_mask, l2reg, self.maxiter, self.tolerance, self.step, self.initial_coefs)

            # Compute the loss on the validation data
            loss = self.loss_fn(model.params,
                                val_c1_data,
                                val_c2_data,
                                val_counts,
                                val_mask)

            # Store the loss as a score (lower score is better)
            scores.append(float(loss))
        return(scores)

    def grid_search_cv(self, l2kfold, l2reg_values=[10**i for i in range(-5, 5)] + [0]):
        """
        Perform grid search cross-validation for hyperparameter tuning.

        Args:
            l2kfold (int): Number of folds for cross-validation.
            l2reg_values (list): List of L2 regularization strengths to search.

        Returns:
            pd.DataFrame: DataFrame with hyperparameter values and corresponding mean cross-validation loss.
        """
        results = pd.DataFrame({'l2reg':[], 'mean_CV_loss':[]})

        for l2reg in l2reg_values:
            # Perform cross-validation
            scores = self.cv_loss(l2kfold, l2reg)

            # Calculate the mean cross-validation score
            mean_score = np.mean(scores)

            temp = pd.DataFrame({'l2reg':[l2reg], 'mean_CV_loss':[mean_score]})
            results = pd.concat([results, temp], ignore_index=True)

        return(results)

    def get_l2reg(self):
        """
        Determine the optimal L2 regularization strength.

        Returns:
            float: Optimal L2 regularization strength.
        """
        self.l2kfold = min(len(np.unique(self.nonrepeat_grp_matrix)), 10)
        grid = self.grid_search_cv(self.l2kfold)
        self.l2reg_grid = grid
        min_obs = grid[grid.mean_CV_loss == min(grid.mean_CV_loss)]
        return(float(min_obs.l2reg.iloc[0]))

    def train_model(self, l2=False, l2reg_value=None, maxiter=None, tolerance=None, step=None):
        """
        Train the conditional logistic regression model with optional L2 regularization.

        Args:
            l2 (bool): Whether to use L2 regularization.
            maxiter (int, optional): Maximum number of optimization iterations. If not provided, the default value is used.
            tolerance (float, optional): tolerance for optimization. If not provided, the default value is used.

        Returns:
            self: The trained ConditionalLogisticRegressor instance.
        """
        assert self.counts_matrix is not None, "counts column is needed"

        if l2 is False:
            self.l2reg = 0
        else:
            if l2reg_value is None:
                self.l2reg = self.get_l2reg()
            else:
                self.l2reg = l2reg_value

        if maxiter is None:
            maxiter = self.maxiter
        else:
            self.maxiter = maxiter

        if tolerance is None:
            tolerance = self.tolerance
        else:
            self.tolerance = tolerance

        if step is None:
            step = self.step
        else:
            self.step = step

        res = self.fit(self.choice1_variable_matrix,
                       self.choice2_variable_matrix,
                       self.counts_matrix,
                       self.mask_matrix,
                       self.l2reg,
                       maxiter,
                       tolerance,
                       step,
                       self.initial_coefs)

        self.coefs = res.params
        self.training_info = res
        self.maxiter = maxiter
        self.tolerance = tolerance
        self.step = step
        self.cov_matrix = self.get_cov_matrix(self.coefs, self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        self.standard_errors = self.get_errors(self.coefs, self.choice1_variable_matrix, self.choice2_variable_matrix, self.counts_matrix, self.mask_matrix, self.l2reg)
        return self

    def get_hessian(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
         # Wrapper function
        def wrapper_loss_fn(coefs):
            return self.loss_fn(coefs, choice1_variables, choice2_variables, counts, mask, l2reg)

        hessian_fn = jax.hessian(wrapper_loss_fn, argnums=0)
        hessian_matrix = hessian_fn(coefs.reshape(-1))
        return hessian_matrix

    def get_cov_matrix(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
        hess_mat = self.get_hessian(coefs, choice1_variables, choice2_variables, counts, mask, l2reg)
        cov_matrix = jnp.linalg.inv(hess_mat)
        return cov_matrix

    def get_errors(self, coefs, choice1_variables, choice2_variables, counts, mask, l2reg=0):
        cov = self.get_cov_matrix(coefs, choice1_variables, choice2_variables, counts, mask, l2reg)
        standard_errors = np.sqrt(np.diag(cov))
        return standard_errors

    def save_model(self, file_path):
        assert self.coefs is not None, "need to train model before saving"
        self.training_df = None
        self.choice1_variable_matrix = None
        self.choice2_variable_matrix = None
        self.counts_matrix = None
        self.nonrepeat_grp_matrix = None
        with open(file_path, 'wb') as file:
            # Serialize and save the object to the file
            dill.dump(self, file)


class TwoStepConditionalLogisticRegressionPredictor(TwoStepDataTransformer):
    """
    A class for making predictions using a trained conditional logistic regression model.

    Args:
        model (ConditionalLogisticRegressor): A trained conditional logistic regression model.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (list): List of column names representing group identifiers.
        choice_colname (list): List of column names representing choice identifiers.

    Attributes:
        model (ConditionalLogisticRegressor): The trained conditional logistic regression model.
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (list): List of original column names representing group identifiers.
        original_choice_colname (list): List of original column names representing choice identifiers.
        variable_matrix (numpy.ndarray): The matrix of predictor variables for prediction.
        counts_matrix (numpy.ndarray): The matrix of counts for prediction.
        nonrepeat_grp_matrix (numpy.ndarray): The matrix of non-repeated group identifiers.

    Methods:
        predict(): Make predictions using the model.
        compute_loss(): Compute the loss for the given data.
        get_coefficients(): Get the model coefficients as a dictionary.
    """
    def __init__(self, model, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params):
        super().__init__(None, variable_colnames, choice1_variable_colnames, choice2_variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, choice2_colname, params)
        self.model = model
        if not isinstance(model, TwoStepConditionalLogisticRegressor):
            raise TypeError("'model' must be a TwoStepConditionalLogisticRegressor object")

    # get probability for input parameters given coefficients
    def predict(self, new_df):
        """
        Make predictions using the trained model.

        Returns:
            pd.DataFrame: A DataFrame containing predicted probabilities merged with the original data.
        Raises:
            ValueError: If the number of variable columns in the input DataFrame doesn't match the model's coefficients.
        """
        original_new_df = new_df
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix, new_df = self.get_matrices(new_df, pretrain=False, return_df=True)

        if not choice1_variable_matrix.shape[-1] + choice2_variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_joint_prob(choice1_variable_matrix, choice2_variable_matrix, mask_matrix, self.model.coefs)
        # transform probs to a dataframe
        choice1_cols = self.get_mapping_dict(new_df, self.choice_colname)
        choice2_cols = self.get_mapping_dict(new_df, self.choice2_colname)
        group_cols = self.get_mapping_dict(new_df, self.group_colname)
        if len(group_cols) == 1:
            probs = probs.reshape((len(group_cols), len(choice1_cols), len(choice2_cols)))

        # Convert the 3D array to a 2D array where each "page" is flattened out
        # and then stack them vertically
        array_2d = probs.reshape(-1, probs.shape[2])

        # Create a MultiIndex representing the first two dimensions (depth and rows)
        # np.repeat and np.tile help to properly align the indices with the reshaped array
        index = pd.MultiIndex.from_tuples([(i, j) for i in range(probs.shape[0]) for j in range(probs.shape[1])], names=[self.group_colname, self.choice_colname])

        # Create the DataFrame
        prob_df = pd.DataFrame(array_2d, index=index, columns=[i for i in range(probs.shape[2])])
        prob_df = prob_df.reset_index()

        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname, self.choice_colname],
                            var_name=self.choice2_colname,
                            value_name='predicted_prob')

        # Invert the dictionary
        choice1_mapping = {v: k for k, v in choice1_cols.items()}
        choice2_mapping = {v: k for k, v in choice2_cols.items()}
        group_mapping = {v: k for k, v in group_cols.items()}

        # Use the map function to replace values in 'Column1'
        melted_df[self.group_colname] = melted_df[self.group_colname].map(group_mapping)
        melted_df[self.choice_colname] = melted_df[self.choice_colname].map(choice1_mapping)
        melted_df[self.choice2_colname] = melted_df[self.choice2_colname].map(choice2_mapping)


        # merge predicted probabilities with original df
        merged_df = pd.merge(new_df, melted_df,
                             on=[self.group_colname, self.choice_colname, self.choice2_colname],
                             how='inner')

        # fill in out of domain values
        common_cols = list(set(original_new_df) & set(merged_df))
        final_df = pd.merge(merged_df, original_new_df, how='outer', on=common_cols)

        # Fill missing values with a specified value
        final_df.fillna(0, inplace=True)
        return(final_df)

    def compute_loss(self, new_df):
        """
        Compute the loss for the given data.

        Returns:
            float: The computed loss.
        """
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(new_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(float(loss))


class TwoStepConditionalLogisticRegressionEvaluator(TwoStepDataTransformer):
    def __init__(self, model_path, params, training_df = None, validation_df = None):
        self.model = self.load_model(model_path)
        self.model.training_df = training_df
        self.validation_df = validation_df
        if not isinstance(self.model, TwoStepConditionalLogisticRegressor):
            raise TypeError("'model' must be a TwoStepConditionalLogisticRegressor object")
        super().__init__(self.model.training_df, self.model.input_variable_colnames, self.model.input_choice1_variable_colnames, self.model.input_choice2_variable_colnames, self.model.count_colname, self.model.input_group_colname, self.model.repeat_obs_colname, self.model.input_choice_colname, self.model.input_choice2_colname, params)
        self.log_loss = None
        self.expected_log_loss = None

    def load_model(self, file_path):
        with open(file_path, 'rb') as file:
            model = dill.load(file)
        assert model.coefs is not None, "model is not trained"
        return(model)

    def calculate_log_loss(self):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the loss on the training data
        loss = self.model.loss_fn(self.model.coefs,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def calculate_expected_log_loss(self, fold_count=20):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the expected loss on the training data
        self.model.choice1_variable_matrix = choice1_variable_matrix
        self.model.choice2_variable_matrix = choice2_variable_matrix
        self.model.counts_matrix = counts_matrix
        self.model.nonrepeat_grp_matrix = nonrepeat_grp_matrix

        expected = self.model.cv_loss(fold_count, self.model.l2reg)
        e_loss = sum((1/fold_count) * np.array(expected))
        return(e_loss)

    def calculate_validation_log_loss(self):
        assert self.validation_df is not None, 'No input validation dataframe provided'
        choice1_variable_matrix, choice2_variable_matrix, counts_matrix, nonrepeat_grp_matrix, mask_matrix = self.get_matrices(self.validation_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  choice1_variable_matrix,
                                  choice2_variable_matrix,
                                  counts_matrix,
                                  mask_matrix)
        return(loss)

    def compile_evaluation_results_df(self, calculate_validation_loss = False, calculate_expected_loss=False):
        result = {'training_annotation_type':[self.params.annotation_type],
                  'productivity':[self.params.productivity],
                  'motif_length_5_end':[self.params.left_nuc_motif_count],
                  'motif_length_3_end':[self.params.right_nuc_motif_count],
                  'motif_type':[self.params.motif_type],
                  'gene_weight_type':[self.params.gene_weight_type],
                  'upper_bound':[self.params.upper_trim_bound],
                  'lower_bound':[self.params.lower_trim_bound],
                  'insertion_bound':[self.params.insertions],
                  'model_type':[self.params.model_type],
                  'base_count_5end_length':[10],
                  'model_parameter_count':[len(self.model.coefs)]}

        results_df = pd.DataFrame(result)

        final = pd.DataFrame()

        if calculate_expected_loss is True:
            self.expected_log_loss = self.calculate_expected_log_loss()

            # add expected loss result
            e = results_df.copy()
            e['loss_type'] = 'Expected log loss across training data'
            e['log_loss'] = self.expected_log_loss

            final = pd.concat([final, e], axis = 0)

        if calculate_validation_loss is True:
            val_loss = self.calculate_validation_log_loss()
            val = results_df.copy()
            val['loss type'] = 'Log loss on validation data'
            val['log_loss'] = val_loss

            final = pd.concat([final, val], axis = 0)

        return(final)
