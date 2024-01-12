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

class DataPreprocessor():
    """
    A class for preprocessing and transforming data.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.

    Attributes:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.

    Methods:
        get_mapping_dict(training_df, col):
            Create a mapping dictionary for unique values in a column.

        transform_categorical_response_vars(training_df, col, new_col):
            Transform a categorical variable into integers.

        get_contrast_matrix(training_df, col):
            Get the contrast matrix for categorical variables.

        get_dropped_contrast_var(training_df, original_col):
            Get the names of dropped contrast variables after transformation.

        transform_categorical_vars(training_df, col):
            Transform categorical variable columns into contrast columns.

        expand_multivariable(training_df, col, new_col):
            Expand multi-variable columns into a single string column.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        self.training_df = training_df
        self.variable_colnames = variable_colnames
        self.count_colname = count_colname
        self.group_colname = group_colname
        self.repeat_obs_colname = repeat_obs_colname
        self.choice_colname = choice_colname
        self.params = params

    def get_mapping_dict(self, training_df, col):
        """
        Create a mapping dictionary for unique values in a column.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the column for which the mapping is generated.

        Returns:
            dict: A dictionary mapping unique values to their corresponding indices.
        """
        # Create a mapping dictionary for unique values in a column
        unique = sorted(list(pd.unique(training_df[col])))
        mapping = {key: index for index, key in enumerate(unique)}
        return mapping

    def transform_categorical_response_vars(self, training_df, col, new_col):
        """
        Transform a categorical variable into integers.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the column to be transformed.
            new_col (str): The name of the new integer-encoded column.

        Returns:
            pd.DataFrame: The DataFrame with the integer-encoded column added.
        """
        # Transform a categorical variable into integers
        col = getattr(self, col)
        if training_df[col].iloc[1] != int:
            new_col_name = col + '_int'
            setattr(self, new_col, new_col_name)
            var_mapping = self.get_mapping_dict(training_df, col)
            training_df[new_col_name] = training_df[col].map(var_mapping)
        return training_df

    def get_contrast_matrix(self, training_df, col, pretrain=True):
        """
        Get the contrast matrix for categorical variables.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the categorical column.

        Returns:
            pd.DataFrame: The contrast matrix for the categorical column.
        """
        first = training_df[col].iloc[0]
        if not isinstance(first, int) and not isinstance(first, float) and not isinstance(first, np.int64):
            # create sum contrasts matrix
            unique = sorted(list(pd.unique(training_df[col])))
            if not pretrain:
                unique = ['-', 'A', 'C', 'G', 'T']

            contrast = Sum().code_without_intercept(unique).matrix

            # Create an identity matrix with the number of categories
            identity_matrix = np.eye(len(unique))

            # Apply the contrast matrix to the identity matrix
            encoded_values = np.dot(identity_matrix, contrast)

            # Create a new DataFrame with the encoded values
            encoded_data = pd.DataFrame(encoded_values, columns=[f"{col}_{i}" for i in unique[:-1]])
            encoded_data[col] = unique
            nonbase_var = col + '_-'
            if nonbase_var in encoded_data.columns:
                encoded_data[nonbase_var] = 0.0
                encoded_data.loc[encoded_data[col] == '-', encoded_data.columns != col] = 0.0
        else:
            encoded_data = pd.DataFrame({col:unique})
        return(encoded_data)

    def get_dropped_contrast_var(self, training_df, original_col):
        """
        Get the names of dropped contrast variables after transformation.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            original_col (str): The name of the original categorical column.

        Returns:
            tuple: A tuple containing two lists - the names of contrast variables and the name of the dropped contrast variable.
        """
        unique = sorted(list(pd.unique(training_df[original_col])))
        contrast_vars=[f"{original_col}_{i}" for i in unique]
        return(contrast_vars[:-1], contrast_vars[-1])

    def transform_categorical_vars(self, training_df, col, pretrain=True):
        """
        Transform categorical variable columns into contrast columns.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the categorical column to be transformed.

        Returns:
            pd.DataFrame: The DataFrame with categorical variables transformed into contrast columns.
        """
        assert col in self.variable_colnames, "Input column name is not a variable name"
        first = training_df[col].iloc[0]
        if not isinstance(first, int) and not isinstance(first, float) and not isinstance(first, np.int64):
            contrast_df = self.get_contrast_matrix(training_df, col, pretrain)
            training_df = pd.merge(training_df, contrast_df, on=col, how='inner')
            new_cols = [x for x in list(contrast_df.columns) if x != col]
            self.variable_colnames = [x for x in self.variable_colnames if x != col] + new_cols
        return(training_df)

    def expand_multivariable(self, training_df, col, new_col):
        """
        Expand multi-variable columns into a single string column.

        Args:
            training_df (pd.DataFrame): The DataFrame containing the data.
            col (str): The name of the multi-variable column.
            new_col (str): The name of the new string column.

        Returns:
            pd.DataFrame: The DataFrame with the multi-variable column expanded into a string column.
        """
        # Expand multi-variable columns into a single string column
        if type(getattr(self, col)) == list:
            new_col_name = '_'.join(getattr(self, col))
            setattr(self, new_col, new_col_name)
            training_df[new_col_name] = training_df.parallel_apply(lambda row: '_'.join(row[getattr(self, col)].astype(str)), axis=1)
            setattr(self, col, new_col_name)
        return training_df

    def check_within_set_variance(self, training_df, pretrain=True):
        if pretrain:
            for col in self.variable_colnames:
                var_counts = training_df.groupby([self.group_colname, col]).size().reset_index(name = 'N')
                unique_var_counts = var_counts.groupby([self.group_colname]).size().reset_index(name = 'N')
                if 1 in unique_var_counts.N.unique():
                    self.variable_colnames.remove(col)
                    self.variable_colnames = self.variable_colnames
                    print('removing ' + col + ' from model due to insufficient within-choice set variance')
        else:
            # remove missing NT variables
            self.variable_colnames = [var for var in self.variable_colnames if '_-' not in var]


    def remove_zero_set_counts(self, training_df):
        count_sums = training_df.groupby([self.group_colname])[self.count_colname].sum()
        if 0 in count_sums.unique():
            zeros = count_sums[count_sums == 0]
            zero_groups = list(zeros.index)
            training_df = training_df[~training_df[self.group_colname].isin(zero_groups)]
            print('removing ' + str(zero_groups) + ' groups from model due to zero counts')
        return(training_df)

    def filter_input_domain_space(self, df):
        domain_file = self.params.R_input_domain_data_path()
        domain_data = pd.read_csv(domain_file, sep = '\t')

        # fill in zeros
        # Create a DataFrame of all unique pairs
        unique_genes = domain_data[['v_gene_group', 'j_gene_group']].drop_duplicates()
        unique_trims = domain_data[['v_trim', 'j_trim']].drop_duplicates()

        unique_pairs = unique_genes.merge(unique_trims, how='cross')
        unique_pairs['ligation_mh'] = 0

        # merge
        domain_data_subset = domain_data[unique_pairs.columns]
        tog = domain_data_subset.merge(unique_pairs, how='outer').drop_duplicates()

        # filter input df
        filtered_df = pd.merge(tog, df, how='inner', on=tog.columns.tolist())

        return filtered_df


class DataTransformer(DataPreprocessor):
    """
    A class for data preprocessing and transformation.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.

    Attributes:
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (str): The original column name representing group identifiers.
        original_choice_colname (str): The original column name representing choice identifiers.
        coefs (ndarray): Coefficients for the model.

    Methods:
        get_random_coefs():
            Generates random coefficients for model initialization.

        preprocess_data():
            Preprocesses the training data, transforming columns as necessary.

        get_matrices():
            Prepares data matrices for modeling.

        get_coefficients(coefs=None):
            Retrieves model coefficients as a dictionary.

    Inherits Attributes and Methods from DataPreprocessor class.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname
        self.coefs = None

    def get_random_coefs(self):
        """
        Generates random coefficients for model initialization.

        Returns:
            ndarray: Randomly generated coefficients as a NumPy array.
        """
        # Set random coefficients for model initialization
        key = jax.random.PRNGKey(123)
        coefs = jax.random.normal(key, shape=(len(self.variable_colnames), 1))
        return coefs

    def preprocess_data(self, df, pretrain=True):
        """
        Preprocesses the training data, transforming columns as necessary.

        Returns:
            pd.DataFrame: The preprocessed training DataFrame.
        """
        # filter for possible sites
        if 'ligation_mh' in self.input_choice_colname:
            df = self.filter_input_domain_space(df)

        # Transform group column lists into strings
        df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # remove zero counts
        if pretrain:
            df = self.remove_zero_set_counts(df)

        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")

        # Transform choice column lists into strings
        df = self.expand_multivariable(df, "original_choice_colname", "choice_colname")
        nonrepeat_groups = self.group_colname

        if self.repeat_obs_colname != None:
            self.group_colname = [self.group_colname] + [self.repeat_obs_colname]
            self.original_group_colname = self.group_colname
            df = self.expand_multivariable(df, "original_group_colname", "group_colname")

        # Transform group and choice columns into integer type
        df = self.transform_categorical_response_vars(df, "original_group_colname", "group_colname")
        df = self.transform_categorical_response_vars(df, "original_choice_colname", "choice_colname")

        # Transform categorical variable columns into contrast columns
        for col in self.variable_colnames:
            df = self.transform_categorical_vars(df, col, pretrain)

        # check for within-choice-set variance
        self.check_within_set_variance(df, pretrain)

        return(df)

    def reset_weighted_observations(self, counts_mat):
        mat_sum = jnp.nansum(counts_mat)
        reset_mat = counts_mat/mat_sum
        return(reset_mat)

    def get_matrices(self, df, pretrain=True, replace_object=None, return_df=False):
        """
        Prepares data matrices for modeling.

        Returns:
            tuple: A tuple containing three data matrices - variables, counts, and non-repeat groups.
        """
        df = self.preprocess_data(df, pretrain)

        # Create three-dimensional matrix
        groups = pd.unique(df[self.group_colname])
        choices = pd.unique(df[self.choice_colname])

        shape = (len(groups), len(choices), len(self.variable_colnames))
        counts_shape = (len(groups), len(choices), 1)
        nonrep_groups_shape = (len(groups), 1)

        mat = np.empty(shape)
        counts_mat = np.empty(counts_shape)
        nonrepeat_groups_mat = np.empty(nonrep_groups_shape)
        var_mapping = {key: value for key, value in enumerate(self.variable_colnames)}

        nonrepeat_groups = self.group_colname
        # Fill the 3D arrays with the training_df data
        for i in groups:
            temp1 = df[df[self.group_colname] == i]
            nonrepeat_groups_mat[i,] = pd.unique(temp1[nonrepeat_groups])
            for j in choices:
                temp2 = temp1[temp1[self.choice_colname] == j]
                if self.count_colname in temp2.columns:
                    if temp2.shape[0] == 0:
                        counts_mat[i, j, 0] = 0
                    else:
                        counts_mat[i, j, 0] = float(temp2[self.count_colname].iloc[0])
                else:
                    counts_mat = None
                for k in var_mapping.keys():
                    if temp2.shape[0] == 0:
                        mat[i, j, k] = 0
                    else:
                        mat[i, j, k] = float(temp2[var_mapping[k]].iloc[0])

        if replace_object is not None:
            setattr(self, replace_object, df)

        if return_df:
            return mat, counts_mat, nonrepeat_groups_mat, df
        else:
            return mat, counts_mat, nonrepeat_groups_mat

    def get_coefficients(self, coefs=None):
        if coefs is None:
            coefs = self.coefs
        coefs = coefs.flatten().tolist()
        d = {n: c for n, c in zip(self.variable_colnames, coefs)}
        for col in list(set(self.original_variable_colnames) - set(self.variable_colnames)):
            contrast_vars, missing_var = self.get_dropped_contrast_var(self.training_df, col)
            for element in contrast_vars:
                if element not in self.variable_colnames:
                    contrast_vars.remove(element)
            d[missing_var] = -1*sum(d[var] for var in contrast_vars)
        return d

    def get_coefficients_df(self):
        """
        Get the model coefficients as a DataFrame.

        Returns:
            dict: A dictionary containing variable names as keys and their corresponding coefficients as values.
        """
        df = pd.DataFrame.from_dict(self.get_coefficients(), columns = ['value'], orient = 'index')
        df['coefficient'] = df.index
        df.reset_index(drop=True, inplace=True)

        # get bases
        df['base'] = None
        df.loc[df.coefficient.str.contains('motif'), 'base'] = df.coefficient.str.split('_').str[-1]
        df.loc[df.coefficient.str.contains('base_count'), 'base'] = df.coefficient.str.split('_').str[-1]
        df.loc[df.coefficient.str.contains('base_count') & df.coefficient.str.contains('prop'), 'base'] = df.coefficient.str.split('_').str[-2]


        # get positions
        df['position'] = None
        df.loc[df.coefficient.str.contains('motif'), 'position'] = df.coefficient.str.split('_').str[-2]
        df.loc[df.coefficient.str.contains('mh_prop'), 'position'] = df.coefficient.str.split('_').str[3] + df.coefficient.str.split('_').str[4]
        df.loc[df.coefficient.str.contains('mh_count'), 'position'] = df.coefficient.str.split('_').str[3] + df.coefficient.str.split('_').str[4]


        # get side
        df['side'] = None
        df.loc[df.coefficient.str.contains('motif'), 'side'] = df.coefficient.str.split('_').str[3]
        df.loc[df.coefficient.str.contains('base_count'), 'side'] = df.coefficient.str.split('_').str[2]
        df.loc[df.coefficient.str.contains('mh_prop'), 'side'] = df.coefficient.str.split('_').str[2]
        df.loc[df.coefficient.str.contains('mh_count'), 'side'] = df.coefficient.str.split('_').str[2]

        # fix trimming specific
        df['trim_type'] = None
        df.loc[df.coefficient.str.contains('_trim'), 'trim_type'] = df.coefficient.str.split('_trim').str[0] + '_trim'
        df.loc[df.coefficient.str.contains('_length'), 'trim_type'] = df.coefficient.str.split('_length').str[0].str.split('_').str[-1] + '_trim'

        # simplify coefficients
        df.loc[df.coefficient.str.contains('interaction') & df.coefficient.str.contains('mh_prop'), 'coefficient'] = 'mh_prop_length_interaction'
        df.loc[df.coefficient.str.contains('interaction') & df.coefficient.str.contains('mh_count'), 'coefficient'] = 'mh_count_length_interaction'
        df.loc[df.coefficient.str.contains('motif'), 'coefficient'] = 'motif'
        df.loc[df.coefficient.str.contains('base_count'), 'coefficient'] = 'base_count'
        df.loc[df.coefficient.str.contains('mh_prop') & ~df.coefficient.str.contains('interaction'), 'coefficient'] = 'mh_prop'
        df.loc[df.coefficient.str.contains('mh_count') & ~df.coefficient.str.contains('interaction'), 'coefficient'] = 'mh_count'

        return(df)



class ConditionalLogisticRegressor(DataTransformer):
    """
    A class for training a conditional logistic regression model and performing cross-validation for hyperparameter tuning.

    Args:
        training_df (pd.DataFrame): The training DataFrame.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (str): The column name representing group identifiers.
        repeat_obs_colname (str): The column name representing repeated observation identifiers.
        choice_colname (str): The column name representing choice identifiers.
        l2kfold (int): Number of folds for L2 regularization hyperparameter tuning.

    Attributes:
        original_variable_colnames (list): List of original variable column names.
        original_group_colname (str): The original column name representing group identifiers.
        original_choice_colname (str): The original column name representing choice identifiers.
        variable_matrix (ndarray): Data matrix of variables.
        counts_matrix (ndarray): Data matrix of counts.
        nonrepeat_grp_matrix (ndarray): Data matrix of non-repeating groups.
        initial_coefs (ndarray): Initial coefficients for model training.
        coefs (ndarray): Trained model coefficients.
        training_info: Information about the training process.
        maxiter (int): Maximum number of iterations for optimization.
        tolerance (float): tolerance for optimization.
        l2reg (float): L2 regularization strength.
        l2kfold (int): Number of folds for L2 regularization hyperparameter tuning.

    Methods:
        get_prob(variables, coefs=None):
            Computes choice probabilities for given variables and coefficients.

        cross_entropy(probs, counts):
            Calculates the cross-entropy loss.

        l2regularization(coefs, size, l2reg):
            Computes L2 regularization term.

        loss_fn(coefs, variables, counts, l2reg=0):
            Computes the loss function for optimization.

        fit(variable_matrix, counts_matrix, l2reg, maxiter, tolerance, initial_coefs):
            Fits the conditional logistic regression model using gradient descent.

        grid_search_cv(l2kfold, l2reg_values=[10**i for i in range(-2, 8)] + [0]):
            Performs grid search cross-validation for hyperparameter tuning.

        get_l2reg():
            Determines the optimal L2 regularization strength.

        train_model(l2=False, maxiter=None, tolerance=None):
            Trains the conditional logistic regression model with optional L2 regularization.

    Inherits Attributes and Methods from DataTransformer class.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params, l2kfold=10):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname
        self.variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix = self.get_matrices(training_df, replace_object='training_df')
        self.initial_coefs = self.get_random_coefs()
        self.coefs = None
        self.training_info = None
        self.maxiter = 1000
        self.tolerance = 1e-8
        self.step = 0.1
        self.l2reg = 0
        self.l2kfold = None
        self.l2reg_grid = None

    # Get probability for input parameters given coefficients
    def get_prob(self, variables, coefs=None):
        """
        Compute choice probabilities for given variables and coefficients.

        Args:
            variables (ndarray): Data matrix of variables.
            coefs (ndarray, optional): Coefficients for the logistic regression model. If not provided, the trained coefficients will be used.

        Returns:
            ndarray: Choice probabilities for each choice given the variables.
        """
        if coefs is None:
            coefs = self.coefs
        assert coefs is not None, "Need to train the model before making predictions!"
        # Compute the logits for each choice
        cov = jnp.dot(variables, coefs)
        reshape = jnp.squeeze(cov)

        # Calculate the probability of the observed choices
        # Dimensions of this matrix are groups x choices
        # replace missing choices with -INF so that they will not count towards probability
        probs = jnp.where(reshape == 0, 0, jax.nn.softmax(jnp.where(reshape==0, jnp.NINF, reshape)))
        return probs

    # Get cross-entropy loss
    def cross_entropy(self, probs, counts):
        """
        Calculate the cross-entropy loss.

        Args:
            probs (ndarray): Choice probabilities.
            counts (ndarray): Counts of choices.

        Returns:
            float: Cross-entropy loss.
        """
        counts_reshape = jnp.squeeze(counts)
        loss = -jnp.sum(jnp.log(jnp.where(probs==0, 1, probs)) * counts_reshape)
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
        var_list = [var for var in self.variable_colnames if 'base_count' not in var or 'interaction' in var]
        var_list = [var for var in var_list if 'motif' not in var]

        mh_list = [var in var_list for var in self.variable_colnames]

        if any(mh_list):
            binary_mh_list = [int(b) for b in mh_list]
            binary_mh_jnp = jnp.array(binary_mh_list).reshape(-1, 1)
            coef_subset = coefs * binary_mh_jnp
            c = jnp.nansum(coef_subset**2)
        else:
            c = 0
        return(0.5*(1/size)*l2reg*c)

    # Compute the loss function
    def loss_fn(self, coefs, variables, counts, l2reg=0):
        """
        Compute the loss function for optimization.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            variables (ndarray): Data matrix of variables.
            counts (ndarray): Counts of choices.
            l2reg (float): L2 regularization strength.

        Returns:
            float: Total loss, including cross-entropy and L2 regularization.
        """
        probs = self.get_prob(variables, coefs)
        if probs.ndim == 1:
            probs = probs.reshape((probs.shape[0], 1))
        size = probs.shape[0]*probs.shape[1]
        loss = self.cross_entropy(probs, counts) + self.l2regularization(coefs, size, l2reg)
        return loss

    def fit(self, variable_matrix, counts_matrix, l2reg, maxiter, tol, step, initial_coefs):
        """
        Fit the conditional logistic regression model using gradient descent.

        Args:
            variable_matrix (ndarray): Data matrix of variables.
            counts_matrix (ndarray): Counts of choices.
            l2reg (float): L2 regularization strength.
            maxiter (int): Maximum number of optimization iterations.
            tolerance (float): tolerance for optimization.
            initial_coefs (ndarray): Initial coefficients for model training.

        Returns:
            OptimizationResult: Result of the optimization process.
        """
        assert counts_matrix is not None, "counts column is missing"

        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.GradientDescent(fun=self.loss_fn, maxiter=maxiter, tol=tol, stepsize=step, implicit_diff=True, verbose=True)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         variables=variable_matrix,
                         counts=counts_matrix,
                         l2reg=l2reg)
        return(res)

    def cv_loss(self, fold_count, l2reg):
        assert self.counts_matrix is not None, "counts column is needed"

        kf = GroupKFold(n_splits=fold_count)
        scores = []

        for train_index, val_index in kf.split(X=self.variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
            train_data, val_data = self.variable_matrix[train_index], self.variable_matrix[val_index]
            train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]
            train_counts = self.reset_weighted_observations(train_counts)
            val_counts = self.reset_weighted_observations(val_counts)

            # Train the model on the training data
            model = self.fit(train_data, train_counts, l2reg, self.maxiter, self.tolerance, self.step, self.initial_coefs)

            # Compute the loss on the validation data
            loss = self.loss_fn(model.params,
                                val_data,
                                val_counts)

            # Store the loss as a score (lower score is better)
            scores.append(float(loss))
        return(scores)

    def grid_search_cv(self, l2kfold, l2reg_values=[10**i for i in range(-2, 10)] + [0]):
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

        res = self.fit(self.variable_matrix,
                       self.counts_matrix,
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
        return self

    def save_model(self, file_path):
        assert self.coefs is not None, "need to train model before saving"
        self.training_df = None
        self.variable_matrix = None
        self.counts_matrix = None
        self.nonrepeat_grp_matrix = None
        with open(file_path, 'wb') as file:
            # Serialize and save the object to the file
            dill.dump(self, file)


class ConditionalLogisticRegressionPredictor(DataTransformer):
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
    def __init__(self, model, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params):
        super().__init__(None, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, params)
        self.model = model
        if not isinstance(model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.input_variable_colnames = variable_colnames
        self.input_group_colname = group_colname
        self.input_choice_colname = choice_colname

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
        variable_matrix, counts_matrix, nonrepeat_grp_matrix, new_df = self.get_matrices(new_df, pretrain=False, return_df=True)

        if not variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_prob(variable_matrix, self.model.coefs)
        # transform probs to a dataframe
        choice_cols = self.get_mapping_dict(new_df, self.choice_colname)
        group_cols = self.get_mapping_dict(new_df, self.group_colname)
        if len(group_cols) == 1:
            probs = probs.reshape((len(group_cols), len(choice_cols)))
        prob_df = pd.DataFrame(probs)
        prob_df.columns = list(choice_cols.keys())
        prob_df[self.group_colname] = list(group_cols.keys())
        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname],
                            var_name=self.choice_colname,
                            value_name='predicted_prob')
        # merge predicted probabilities with original df
        merged_df = pd.merge(new_df, melted_df,
                             on=[self.group_colname, self.choice_colname],
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
        variable_matrix, counts_matrix, nonrepeat_grp_matrix = self.get_matrices(new_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix)
        return(float(loss))


class ConditionalLogisticRegressionEvaluator(DataTransformer):
    def __init__(self, model_path, params, training_df = None, validation_df = None):
        self.model = self.load_model(model_path)
        self.model.training_df = training_df
        self.validation_df = validation_df
        if not isinstance(self.model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        super().__init__(self.model.training_df, self.model.input_variable_colnames, self.model.count_colname, self.model.input_group_colname, self.model.repeat_obs_colname, self.model.input_choice_colname, params)
        self.log_loss = None
        self.expected_log_loss = None

    def load_model(self, file_path):
        with open(file_path, 'rb') as file:
            model = dill.load(file)
        assert model.coefs is not None, "model is not trained"
        return(model)

    def calculate_log_loss(self):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the loss on the training data
        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix)
        return(loss)

    def calculate_expected_log_loss(self, fold_count=20):
        assert self.model.training_df is not None, 'No input training dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix = self.get_matrices(self.model.training_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        # Compute the expected loss on the training data
        self.model.variable_matrix = variable_matrix
        self.model.counts_matrix = counts_matrix
        self.model.nonrepeat_grp_matrix = nonrepeat_grp_matrix

        expected = self.model.cv_loss(fold_count, self.model.l2reg)
        e_loss = sum((1/fold_count) * np.array(expected))
        return(e_loss)

    def calculate_validation_log_loss(self):
        assert self.validation_df is not None, 'No input validation dataframe provided'
        variable_matrix, counts_matrix, nonrepeat_grp_matrix = self.get_matrices(self.validation_df, pretrain=False)

        assert counts_matrix is not None, "counts column is needed"

        loss = self.model.loss_fn(self.model.coefs,
                                  variable_matrix,
                                  counts_matrix)
        return(loss)

    def compile_evaluation_results_df(self, left_nuc_count, right_nuc_count, motif_type, gene_weight_type, upper_trim_bound, lower_trim_bound, insertion_bound, model_type, base_count_5end_length, calculate_validation_loss = False):
        result = {'motif_length_5_end':[left_nuc_count],
                  'motif_length_3_end':[right_nuc_count],
                  'motif_type':[motif_type],
                  'gene_weight_type':[gene_weight_type],
                  'upper_bound':[upper_trim_bound],
                  'lower_bound':[lower_trim_bound],
                  'insertion_bound':[insertion_bound],
                  'model_type':[model_type],
                  'base_count_5end_length':[base_count_5end_length],
                  'model_parameter_count':[len(self.model.coefs)]}

        results_df = pd.DataFrame(result)

        if calculate_validation_loss is False:
            self.log_loss = self.calculate_log_loss()
            self.expected_log_loss = self.calculate_expected_log_loss()

            # add loss result
            l = results_df.copy()
            l['loss_type'] = 'Log loss on training data'
            l['log_loss'] = self.log_loss

            # add expected loss result
            e = results_df.copy()
            e['loss_type'] = 'Expected log loss across training data'
            e['log_loss'] = self.expected_log_loss

            final = pd.concat([l, e], axis = 0)
        else:
            val_loss = self.calculate_validation_log_loss()
            final = result.copy()
            final['loss type'] = 'Log loss on validation data'
            final['log_loss'] = val_loss
        return(final)


