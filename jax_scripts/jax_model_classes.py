import pandas as pd
import numpy as np
import jax
import jax.numpy as jnp
import jaxopt
import patsy
from patsy.contrasts import Sum
from sklearn.model_selection import GroupKFold

# TODO move on to more complex case
# TODO make compatible with up and downstream R code

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
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname):
        self.training_df = training_df
        self.variable_colnames = variable_colnames
        self.count_colname = count_colname
        self.group_colname = group_colname
        self.repeat_obs_colname = repeat_obs_colname
        self.choice_colname = choice_colname

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

    def get_contrast_matrix(self, training_df, col):
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
            contrast = Sum().code_without_intercept(unique).matrix

            # Create an identity matrix with the number of categories
            identity_matrix = np.eye(len(unique))

            # Apply the contrast matrix to the identity matrix
            encoded_values = np.dot(identity_matrix, contrast)

            # Create a new DataFrame with the encoded values
            encoded_data = pd.DataFrame(encoded_values, columns=[f"{col}_{i}" for i in unique[:-1]])
            encoded_data[col] = unique
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

    def transform_categorical_vars(self, training_df, col):
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
            contrast_df = self.get_contrast_matrix(training_df, col)
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
            training_df[new_col_name] = training_df.apply(lambda row: '_'.join(row[getattr(self, col)].astype(str)), axis=1)
            setattr(self, col, new_col_name)
        return training_df


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
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname)
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
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

    def preprocess_data(self):
        """
        Preprocesses the training data, transforming columns as necessary.

        Returns:
            pd.DataFrame: The preprocessed training DataFrame.
        """
        # Transform group and choice column lists into strings
        self.training_df = self.expand_multivariable(self.training_df, "original_group_colname", "group_colname")
        self.training_df = self.transform_categorical_response_vars(self.training_df, "original_group_colname", "group_colname")
        self.training_df = self.expand_multivariable(self.training_df, "original_choice_colname", "choice_colname")
        nonrepeat_groups = self.group_colname

        if self.repeat_obs_colname != None:
            self.group_colname = [self.group_colname] + [self.repeat_obs_colname]
            self.original_group_colname = self.group_colname
            self.training_df = self.expand_multivariable(self.training_df, "original_group_colname", "group_colname")

        # Transform group and choice columns into integer type
        self.training_df = self.transform_categorical_response_vars(self.training_df, "original_group_colname", "group_colname")
        self.training_df = self.transform_categorical_response_vars(self.training_df, "original_choice_colname", "choice_colname")

        # Transform categorical variable columns into contrast columns
        for col in self.variable_colnames:
            self.training_df = self.transform_categorical_vars(self.training_df, col)
        return(self.training_df)

    def get_matrices(self):
        """
        Prepares data matrices for modeling.

        Returns:
            tuple: A tuple containing three data matrices - variables, counts, and non-repeat groups.
        """
        self.training_df = self.preprocess_data()

        # Create three-dimensional matrix
        groups = pd.unique(self.training_df[self.group_colname])
        choices = pd.unique(self.training_df[self.choice_colname])

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
            temp1 = self.training_df[self.training_df[self.group_colname] == i]
            nonrepeat_groups_mat[i,] = pd.unique(temp1[nonrepeat_groups])
            for j in choices:
                temp2 = temp1[temp1[self.choice_colname] == j]
                counts_mat[i, j, 0] = float(temp2[self.count_colname].iloc[0])
                for k in var_mapping.keys():
                    mat[i, j, k] = float(temp2[var_mapping[k]].iloc[0])
        return mat, counts_mat, nonrepeat_groups_mat

    def get_coefficients(self, coefs=None):
        if coefs is None:
            coefs = self.coefs
        coefs = coefs.flatten().tolist()
        d = {n: c for n, c in zip(self.variable_colnames, coefs)}
        for col in list(set(self.original_variable_colnames) - set(self.variable_colnames)):
            contrast_vars, missing_var = self.get_dropped_contrast_var(self.training_df, col)
            d[missing_var] = sum(d[var] for var in contrast_vars)
        return d


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
        stepsize (float): Step size for optimization.
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

        fit(variable_matrix, counts_matrix, l2reg, maxiter, stepsize, initial_coefs):
            Fits the conditional logistic regression model using gradient descent.

        grid_search_cv(l2kfold, l2reg_values=[10**i for i in range(-2, 8)] + [0]):
            Performs grid search cross-validation for hyperparameter tuning.

        get_l2reg():
            Determines the optimal L2 regularization strength.

        train_model(l2=False, maxiter=None, stepsize=None):
            Trains the conditional logistic regression model with optional L2 regularization.

    Inherits Attributes and Methods from DataTransformer class.
    """
    def __init__(self, training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname, l2kfold=10):
        super().__init__(training_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname)
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix = self.get_matrices()
        self.initial_coefs = self.get_random_coefs()
        self.coefs = None
        self.training_info = None
        self.maxiter = 100000
        self.stepsize = 0.000001
        self.l2reg = 0
        self.l2kfold = None

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
        probs = jax.nn.softmax(reshape)
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
        loss = -jnp.sum(jnp.log(probs) * counts_reshape)
        return loss

    def l2regularization(self, coefs, size, l2reg):
        """
        Compute L2 regularization term.

        Args:
            coefs (ndarray): Coefficients for the logistic regression model.
            size (int): Size of the training data
            l2reg (float): L2 regularization strength.

        Returns:
            float: L2 regularization term.
        """
        c = jnp.sum(coefs**2)
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

    def fit(self, variable_matrix, counts_matrix, l2reg, maxiter, stepsize, initial_coefs):
        """
        Fit the conditional logistic regression model using gradient descent.

        Args:
            variable_matrix (ndarray): Data matrix of variables.
            counts_matrix (ndarray): Counts of choices.
            l2reg (float): L2 regularization strength.
            maxiter (int): Maximum number of optimization iterations.
            stepsize (float): Step size for optimization.
            initial_coefs (ndarray): Initial coefficients for model training.

        Returns:
            OptimizationResult: Result of the optimization process.
        """
        # Create a jaxopt GradientDescent optimizer
        solver = jaxopt.GradientDescent(fun=self.loss_fn, maxiter=maxiter, stepsize=stepsize)

        # Run gradient descent
        res = solver.run(initial_coefs,
                         variables=variable_matrix,
                         counts=counts_matrix,
                         l2reg=l2reg)
        return(res)

    def grid_search_cv(self, l2kfold, l2reg_values=[10**i for i in range(-2, 8)] + [0]):
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
            kf = GroupKFold(n_splits=l2kfold)
            scores = []

            for train_index, val_index in kf.split(X=self.variable_matrix, y=self.counts_matrix, groups=self.nonrepeat_grp_matrix):
                train_data, val_data = self.variable_matrix[train_index], self.variable_matrix[val_index]
                train_counts, val_counts = self.counts_matrix[train_index], self.counts_matrix[val_index]

                # Train the model on the training data
                model = self.fit(train_data, train_counts, l2reg, self.maxiter, self.stepsize, self.initial_coefs)

                # Compute the loss on the validation data
                loss = self.loss_fn(model.params,
                                    val_data,
                                    val_counts)

                # Store the loss as a score (lower score is better)
                scores.append(float(loss))

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
        min_obs = grid[grid.mean_CV_loss == min(grid.mean_CV_loss)]
        return(float(min_obs.l2reg.iloc[0]))

    def train_model(self, l2=False, maxiter=None, stepsize=None):
        """
        Train the conditional logistic regression model with optional L2 regularization.

        Args:
            l2 (bool): Whether to use L2 regularization.
            maxiter (int, optional): Maximum number of optimization iterations. If not provided, the default value is used.
            stepsize (float, optional): Step size for optimization. If not provided, the default value is used.

        Returns:
            self: The trained ConditionalLogisticRegressor instance.
        """
        if l2 is False:
            self.l2reg = 0
        else:
            self.l2reg = self.get_l2reg()

        if maxiter is None:
            maxiter = self.maxiter
        else:
            self.maxiter = maxiter

        if stepsize is None:
            stepsize = self.stepsize
        else:
            self.stepsize = stepsize

        res = self.fit(self.variable_matrix, self.counts_matrix, self.l2reg, maxiter, stepsize, self.initial_coefs)

        self.coefs = res.params
        self.training_info = res
        self.maxiter = maxiter
        self.stepsize = stepsize
        return self


class ConditionalLogisticRegressionPredictor(DataTransformer):
    """
    A class for making predictions using a trained conditional logistic regression model.

    Args:
        model (ConditionalLogisticRegressor): A trained conditional logistic regression model.
        new_df (pd.DataFrame): The DataFrame with data for prediction.
        variable_colnames (list): List of variable column names.
        count_colname (str): The column name for the count variable.
        group_colname (list): List of column names representing group identifiers.
        choice_colname (list): List of column names representing choice identifiers.

    Attributes:
        model (ConditionalLogisticRegressor): The trained conditional logistic regression model.
        new_df (pd.DataFrame): The DataFrame with data for prediction.
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
    def __init__(self, model, new_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname):
        super().__init__(new_df, variable_colnames, count_colname, group_colname, repeat_obs_colname, choice_colname)
        self.model = model
        if not isinstance(model, ConditionalLogisticRegressor):
            raise TypeError("'model' must be a ConditionalLogisticRegressor object")
        self.new_df = new_df
        self.original_variable_colnames = variable_colnames
        self.original_group_colname = group_colname
        self.original_choice_colname = choice_colname
        self.variable_matrix, self.counts_matrix, self.nonrepeat_grp_matrix = self.get_matrices()

    # get probability for input parameters given coefficients
    def predict(self):
        """
        Make predictions using the trained model.

        Returns:
            pd.DataFrame: A DataFrame containing predicted probabilities merged with the original data.
        Raises:
            ValueError: If the number of variable columns in the input DataFrame doesn't match the model's coefficients.
        """
        if not self.variable_matrix.shape[-1] == self.model.coefs.shape[0]:
            raise ValueError("Input dataframe variable column count doesn't match the trained model coefficient count")
        # get predicted probabilities
        probs = self.model.get_prob(self.variable_matrix, self.model.coefs)
        # transform probs to a dataframe
        prob_df = pd.DataFrame(probs)
        choice_cols = self.get_mapping_dict(self.new_df, self.choice_colname)
        group_cols = self.get_mapping_dict(self.new_df, self.group_colname)
        prob_df.columns = list(choice_cols.keys())
        prob_df[self.group_colname] = list(group_cols.keys())
        melted_df = pd.melt(prob_df,
                            id_vars=[self.group_colname],
                            var_name=self.choice_colname,
                            value_name='predicted_prob')
        # merge predicted probabilities with original df
        merged_df = pd.merge(self.new_df, melted_df,
                             on=[self.group_colname, self.choice_colname],
                             how='inner')
        return(merged_df)

    def compute_loss(self):
        """
        Compute the loss for the given data.

        Returns:
            float: The computed loss.
        """
        loss = self.model.loss_fn(self.model.coefs,
                             self.variable_matrix,
                             self.counts_matrix)
        return(float(loss))

    def get_coefficients(self):
        """
        Get the model coefficients as a dictionary.

        Returns:
            dict: A dictionary containing variable names as keys and their corresponding coefficients as values.
        """
        self.coefs = self.model.coefs
        return(super().get_coefficients())
